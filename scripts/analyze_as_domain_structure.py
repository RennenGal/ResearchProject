#!/usr/bin/env python3
"""
Structural impact of AS events on TIM-barrel domain architecture.

Classifies each canonical TIM-barrel element in every AS isoform two ways:

  Combined: each (beta/alpha) repeat treated as one unit [beta_start, alpha_end]
  Separate: beta-strand [beta_start, beta_end] and alpha-helix [alpha_start, alpha_end]
            classified independently

State per element/motif:
  intact  – no overlap with any VSP replacement region
  partial – partially overlapping
  removed – fully contained within a VSP replacement region

Also reports mean pLDDT per region (pre-VSP / VSP / post-VSP) for
single-VSP isoforms with AlphaFold structures (printed to stdout only).

Output:
  figures/as_domain_disruption.png          – combined (β/α as one motif)
  figures/as_domain_disruption_separate.png – β and α classified independently
  data/domain_disruption_summary.tsv        – per-isoform disruption table

Usage:
    python scripts/analyze_as_domain_structure.py
"""

import argparse
import json
import sqlite3
import sys
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats as scipy_stats

sys.path.insert(0, str(Path(__file__).parent.parent))
from protein_data_collector.config import get_config

PDB_ISO_DIR = Path("data/alphafold_isoforms")


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_canonical(conn):
    rows = conn.execute("""
        SELECT vc.uniprot_id, vc.domain_start, vc.domain_end, vc.motif_annotations,
               i.sequence_length
        FROM   view_canonical vc
        JOIN   isoforms i ON i.uniprot_id = vc.uniprot_id AND i.is_canonical = 1
    """).fetchall()
    out = {}
    for uid, ds, de, mj, slen in rows:
        out[uid] = dict(ds=ds, de=de, motifs=json.loads(mj), slen=slen)
    return out


def load_isoforms(conn):
    rows = conn.execute("""
        SELECT nc.isoform_id, nc.uniprot_id, nc.vsp_domain_events,
               i.sequence_length
        FROM   view_noncanonical nc
        JOIN   isoforms i ON i.isoform_id = nc.isoform_id
        WHERE  nc.vsp_domain_events != '[]'
    """).fetchall()
    out = []
    for iso_id, uid, vsp_json, iso_len in rows:
        out.append(dict(isoform_id=iso_id, uniprot_id=uid,
                        vsps=json.loads(vsp_json), iso_len=iso_len))
    return out


# ---------------------------------------------------------------------------
# Motif disruption classification
# ---------------------------------------------------------------------------

def classify_element(elem_start, elem_end, vsps):
    """Return 'removed', 'partial', or 'intact' for a single span vs. all VSPs."""
    state = "intact"
    for vsp in vsps:
        v_s, v_e = vsp["can_start"], vsp["can_end"]
        if v_e < elem_start or v_s > elem_end:
            continue
        if v_s <= elem_start and v_e >= elem_end:
            return "removed"
        state = "partial"
    return state


def classify_isoform_combined(iso, canonicals):
    """Return list of states, one per whole (β/α) motif [beta_start, alpha_end]."""
    uid = iso["uniprot_id"]
    if uid not in canonicals:
        return None
    vsps = iso["vsps"]
    return [classify_element(m["beta_start"], m["alpha_end"], vsps)
            for m in canonicals[uid]["motifs"]]


def classify_isoform_separate(iso, canonicals):
    """Return list of (beta_state, alpha_state) tuples, one per canonical motif."""
    uid = iso["uniprot_id"]
    if uid not in canonicals:
        return None
    vsps = iso["vsps"]
    return [
        (classify_element(m["beta_start"], m["beta_end"],   vsps),
         classify_element(m["alpha_start"], m["alpha_end"], vsps))
        for m in canonicals[uid]["motifs"]
    ]


# ---------------------------------------------------------------------------
# pLDDT extraction
# ---------------------------------------------------------------------------

def parse_pdb_bfactors(pdb_path):
    """Return dict: residue_number → mean pLDDT (averaged over all atoms)."""
    res_vals = {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            try:
                resnum = int(line[22:26])
                bfact  = float(line[60:66])
            except ValueError:
                continue
            res_vals.setdefault(resnum, []).append(bfact)
    return {r: float(np.mean(v)) for r, v in res_vals.items()}


def region_mean_plddt(bfactors, start, end):
    """Mean pLDDT for residues [start, end] inclusive. NaN if no residues."""
    vals = [bfactors[r] for r in range(start, end + 1) if r in bfactors]
    return float(np.mean(vals)) if vals else float("nan")


# ---------------------------------------------------------------------------
# pLDDT analysis — single-VSP isoforms
# ---------------------------------------------------------------------------

def analyze_plddt(isoforms, canonicals):
    """
    For single-VSP isoforms with a PDB file, extract mean pLDDT for:
      pre  : isoform residues 1 .. can_start-1
      vsp  : isoform residues can_start .. can_start + repl_len - 1
      post : isoform residues can_start + repl_len .. iso_len

    replacement_length = iso_len - can_len + (can_end - can_start + 1)
    """
    records = []
    for iso in isoforms:
        if len(iso["vsps"]) != 1:
            continue
        uid = iso["uniprot_id"]
        if uid not in canonicals:
            continue
        pdb_path = PDB_ISO_DIR / f"{iso['isoform_id']}.pdb"
        if not pdb_path.exists():
            continue

        vsp     = iso["vsps"][0]
        can_s   = vsp["can_start"]
        can_e   = vsp["can_end"]
        can_len = canonicals[uid]["slen"]
        iso_len = iso["iso_len"]

        repl_len = iso_len - can_len + (can_e - can_s + 1)
        if repl_len < 0:
            continue

        bf = parse_pdb_bfactors(pdb_path)

        pre_end    = can_s - 1
        vsp_end    = can_s + repl_len - 1
        post_start = vsp_end + 1

        pre  = region_mean_plddt(bf, 1,          pre_end)    if pre_end >= 1          else float("nan")
        vsp_ = region_mean_plddt(bf, can_s,      vsp_end)    if repl_len > 0          else float("nan")
        post = region_mean_plddt(bf, post_start, iso_len)    if post_start <= iso_len  else float("nan")

        records.append(dict(
            isoform_id=iso["isoform_id"], uniprot_id=uid,
            can_start=can_s, can_end=can_e, repl_len=repl_len,
            pre_plddt=pre, vsp_plddt=vsp_, post_plddt=post,
        ))
    return records


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def _bar_panel(ax, pos_rates_list, pos_labels, title, colors=None):
    """Shared helper for the per-position disruption bar panel."""
    if colors is None:
        colors = ["#4C72B0", "#DD8452", "#55A868"]
    max_pos = len(pos_rates_list[0])
    x = np.arange(1, max_pos + 1)
    n_series = len(pos_rates_list)
    w = 0.7 / n_series
    offsets = np.linspace(-(0.7 - w) / 2, (0.7 - w) / 2, n_series)
    for rates, lbl, col, offset in zip(pos_rates_list, pos_labels, colors, offsets):
        bars = ax.bar(x + offset, rates, width=w, color=col,
                      alpha=0.85, edgecolor="white", zorder=3, label=lbl)
        for bar, r in zip(bars, rates):
            if r >= 5:
                ax.text(bar.get_x() + bar.get_width() / 2, r + 1.0,
                        f"{r:.0f}", ha="center", va="bottom", fontsize=7)
    ax.set_xlabel("Barrel repeat position", fontsize=10)
    ax.set_ylabel("Isoforms with disruption (%)", fontsize=10)
    ax.set_xticks(list(x))
    ax.set_ylim(0, 110)
    ax.set_title(title, fontsize=9)
    if n_series > 1:
        ax.legend(fontsize=9, frameon=False)
    ax.yaxis.grid(True, linestyle=":", alpha=0.4, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_disruption_combined(combined_per_isoform, out):
    """Two panels: intact-motif histogram + per-position rate (whole beta/alpha motif)."""
    n = len(combined_per_isoform)
    intact_counts = [sum(1 for s in st if s == "intact") for st in combined_per_isoform]

    max_pos = max(len(st) for st in combined_per_isoform)
    pos_rate = []
    for i in range(max_pos):
        relevant = [st for st in combined_per_isoform if len(st) > i]
        pos_rate.append(100 * sum(1 for st in relevant if st[i] != "intact") / len(relevant))

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    ax = axes[0]
    max_count = max(intact_counts)
    ax.hist(intact_counts, bins=np.arange(-0.5, max_count + 1.5, 1),
            color="#4C72B0", alpha=0.85, edgecolor="white", zorder=3)
    ax.set_xlabel("Intact (beta/alpha) motifs retained", fontsize=10)
    ax.set_ylabel("Number of isoforms", fontsize=10)
    ax.set_xticks(range(max_count + 1))
    mean_val = float(np.mean(intact_counts))
    ax.axvline(mean_val, color="#C44E52", lw=1.5, ls="--", label=f"Mean = {mean_val:.1f}")
    ax.legend(fontsize=9, frameon=False)
    ax.set_title(f"Intact (beta/alpha) motifs per isoform\nn = {n} isoforms", fontsize=9)
    ax.yaxis.grid(True, linestyle=":", alpha=0.4, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    _bar_panel(axes[1], [pos_rate], ["motif"],
               "Per-position disruption rate\n(beta/alpha as one motif, partial + removed)")

    fig.tight_layout()
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")


def plot_disruption_separate(separate_per_isoform, out):
    """Two panels: intact-element histogram + grouped per-position rates (beta vs alpha)."""
    n = len(separate_per_isoform)
    intact_counts = [
        sum(1 for b, _ in st if b == "intact") + sum(1 for _, a in st if a == "intact")
        for st in separate_per_isoform
    ]

    max_pos = max(len(st) for st in separate_per_isoform)
    beta_rates, alpha_rates = [], []
    for i in range(max_pos):
        relevant = [st for st in separate_per_isoform if len(st) > i]
        n_rel = len(relevant)
        beta_rates.append(100 * sum(1 for st in relevant if st[i][0] != "intact") / n_rel)
        alpha_rates.append(100 * sum(1 for st in relevant if st[i][1] != "intact") / n_rel)

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))

    ax = axes[0]
    max_elem = max(2 * len(st) for st in separate_per_isoform)
    ax.hist(intact_counts, bins=np.arange(-0.5, max_elem + 1.5, 1),
            color="#4C72B0", alpha=0.85, edgecolor="white", zorder=3)
    ax.set_xlabel("Intact elements retained (beta-strand + alpha-helix)", fontsize=10)
    ax.set_ylabel("Number of isoforms", fontsize=10)
    ax.set_xticks(range(0, max_elem + 1, 2))
    mean_val = float(np.mean(intact_counts))
    ax.axvline(mean_val, color="#C44E52", lw=1.5, ls="--", label=f"Mean = {mean_val:.1f}")
    ax.legend(fontsize=9, frameon=False)
    ax.set_title(f"Intact elements per isoform (beta + alpha)\nn = {n} isoforms", fontsize=9)
    ax.yaxis.grid(True, linestyle=":", alpha=0.4, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    _bar_panel(axes[1], [beta_rates, alpha_rates], ["beta-strand", "alpha-helix"],
               "Per-position disruption rate (beta and alpha separated)\n(partial + removed)")

    fig.tight_layout()
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# TSV summary
# ---------------------------------------------------------------------------

def write_summary(isoforms, combined_map, separate_map, out):
    """Write per-isoform TSV with combined and separate disruption states."""
    max_motifs = max((len(v) for v in combined_map.values()), default=8)
    motif_cols = [f"motif{i+1}" for i in range(max_motifs)]
    elem_cols  = []
    for i in range(max_motifs):
        elem_cols += [f"beta{i+1}", f"alpha{i+1}"]
    header = "\t".join(
        ["isoform_id", "uniprot_id", "n_vsps", "n_motifs",
         "n_intact_motif", "n_partial_motif", "n_removed_motif",
         "n_intact_beta", "n_intact_alpha"]
        + motif_cols + elem_cols
    )
    lines = [header]
    for iso in isoforms:
        cm = combined_map.get(iso["isoform_id"])
        sm = separate_map.get(iso["isoform_id"])
        if cm is None or sm is None:
            continue
        flat_sep = []
        for b, a in sm:
            flat_sep += [b, a]
        pad_m = [""] * (max_motifs - len(cm))
        pad_e = [""] * (2 * max_motifs - len(flat_sep))
        lines.append("\t".join([
            iso["isoform_id"], iso["uniprot_id"],
            str(len(iso["vsps"])), str(len(cm)),
            str(sum(1 for s in cm if s == "intact")),
            str(sum(1 for s in cm if s == "partial")),
            str(sum(1 for s in cm if s == "removed")),
            str(sum(1 for b, _ in sm if b == "intact")),
            str(sum(1 for _, a in sm if a == "intact")),
        ] + cm + pad_m + flat_sep + pad_e))
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    Path(out).write_text("\n".join(lines) + "\n")
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db",              default=None)
    parser.add_argument("--out-combined",    default="figures/as_domain_disruption.png")
    parser.add_argument("--out-separate",    default="figures/as_domain_disruption_separate.png")
    parser.add_argument("--out-tsv",         default="data/domain_disruption_summary.tsv")
    args = parser.parse_args()

    db_path    = args.db or get_config().db_path
    conn       = sqlite3.connect(db_path)
    canonicals = load_canonical(conn)
    isoforms   = load_isoforms(conn)
    conn.close()

    print(f"Canonical proteins:       {len(canonicals)}")
    print(f"Isoforms with VSP events: {len(isoforms)}")

    # --- Classify ---
    combined_map, combined_list = {}, []
    separate_map, separate_list = {}, []
    for iso in isoforms:
        cm = classify_isoform_combined(iso, canonicals)
        sm = classify_isoform_separate(iso, canonicals)
        if cm is None:
            continue
        combined_map[iso["isoform_id"]] = cm
        separate_map[iso["isoform_id"]] = sm
        combined_list.append(cm)
        separate_list.append(sm)

    print(f"Classified:               {len(combined_list)}\n")

    # --- Combined stats ---
    intact_combined = [sum(1 for s in st if s == "intact") for st in combined_list]
    print(f"[Combined motif] Mean intact: {np.mean(intact_combined):.2f}  "
          f"Median: {np.median(intact_combined):.1f}")
    print("Intact-motif distribution:")
    for k, v in sorted(Counter(intact_combined).items()):
        print(f"  {k}: {v} isoforms ({100*v/len(intact_combined):.1f}%)")

    max_pos = max(len(st) for st in combined_list)
    print("\nPer-position disruption rate (combined):")
    for i in range(max_pos):
        relevant  = [st for st in combined_list if len(st) > i]
        disrupted = sum(1 for st in relevant if st[i] != "intact")
        print(f"  Position {i+1}: {disrupted}/{len(relevant)} "
              f"({100*disrupted/len(relevant):.1f}%)")

    # --- Separate stats ---
    intact_beta  = [sum(1 for b, _ in st if b == "intact") for st in separate_list]
    intact_alpha = [sum(1 for _, a in st if a == "intact") for st in separate_list]
    print(f"\n[Separate] Mean intact beta-strands: {np.mean(intact_beta):.2f}  "
          f"alpha-helices: {np.mean(intact_alpha):.2f}")
    print("\nPer-position disruption rate (beta / alpha):")
    for i in range(max_pos):
        relevant = [st for st in separate_list if len(st) > i]
        n_rel    = len(relevant)
        b_dis = sum(1 for st in relevant if st[i][0] != "intact")
        a_dis = sum(1 for st in relevant if st[i][1] != "intact")
        print(f"  Position {i+1}: beta={b_dis}/{n_rel} ({100*b_dis/n_rel:.1f}%)  "
              f"alpha={a_dis}/{n_rel} ({100*a_dis/n_rel:.1f}%)")

    plot_disruption_combined(combined_list, args.out_combined)
    plot_disruption_separate(separate_list, args.out_separate)
    write_summary(isoforms, combined_map, separate_map, args.out_tsv)

    # --- pLDDT (summary only, no figure) ---
    print("\n--- pLDDT analysis (single-VSP isoforms) ---")
    plddt_records = analyze_plddt(isoforms, canonicals)
    print(f"Single-VSP isoforms with PDB: {len(plddt_records)}")
    if plddt_records:
        pre_v  = [r["pre_plddt"]  for r in plddt_records if not np.isnan(r["pre_plddt"])]
        vsp_v  = [r["vsp_plddt"]  for r in plddt_records if not np.isnan(r["vsp_plddt"])]
        post_v = [r["post_plddt"] for r in plddt_records if not np.isnan(r["post_plddt"])]
        print(f"  Pre-VSP  mean pLDDT: {np.mean(pre_v):.1f}  (n={len(pre_v)})")
        print(f"  VSP reg. mean pLDDT: {np.mean(vsp_v):.1f}  (n={len(vsp_v)})")
        print(f"  Post-VSP mean pLDDT: {np.mean(post_v):.1f}  (n={len(post_v)})")
        for a, b, lbl in [(pre_v, vsp_v, "pre vs VSP"),
                          (vsp_v, post_v, "VSP vs post"),
                          (pre_v, post_v, "pre vs post")]:
            _, p = scipy_stats.mannwhitneyu(a, b, alternative="two-sided")
            stars = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
            print(f"  {lbl}: p={p:.4f} ({stars})")


if __name__ == "__main__":
    main()
