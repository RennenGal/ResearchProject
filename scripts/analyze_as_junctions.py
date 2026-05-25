#!/usr/bin/env python3
"""
AS junction placement analysis (§8 of Statistical-Analysis.md).

  8A. Element-level distribution of AS-affected canonical junctions vs baseline
  8B. Cross-protein positional distribution (KS test vs canonical baseline)
  8C. Within-protein hotspot reuse across multiple AS isoforms

Output figures:
  figures/as_element_enrichment.png
  figures/as_positional.png
  figures/as_hotspot.png

Usage:
    python scripts/analyze_as_junctions.py
    python scripts/analyze_as_junctions.py --B 5000
"""

import argparse
import json
import sqlite3
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent.parent))
from protein_data_collector.config import get_config
from junction_utils import load_canonical_junctions, load_isoform_junctions


# ---------------------------------------------------------------------------
# Category definitions (identical to analyze_junction_enrichment.py)
# ---------------------------------------------------------------------------

CATS5    = ["beta", "alpha", "inter", "loop", "flanking"]
LABELS5  = {"beta": "β-strand", "alpha": "α-helix", "inter": "Inter-motif",
            "loop": "Loop (β→α)", "flanking": "Flanking"}
ALABELS5 = {"beta": "beta-strand", "alpha": "alpha-helix", "inter": "Inter-motif",
            "loop": "Loop (b->a)", "flanking": "Flanking"}
COLS5    = {"beta": "#4C72B0", "alpha": "#DD8452", "inter": "#55A868",
            "loop": "#8c8c8c", "flanking": "#C44E52"}

MIN_MATCH = 15   # minimum AA window for suffix match to find resync
MAX_SLIDE = 5    # ±residues to slide can_end when searching for resync


# ---------------------------------------------------------------------------
# Element assignment — τ_5, identical to analyze_junction_enrichment.py
# ---------------------------------------------------------------------------

def _tau5(pos, motifs):
    """5-category label for pos; first-match, lower motif index wins."""
    if not motifs or pos < motifs[0]["beta_start"]:
        return "flanking"
    for i, m in enumerate(motifs):
        if m["beta_start"] <= pos <= m["beta_end"]:
            return "beta"
        if m["beta_end"] < pos < m["alpha_start"]:
            return "loop"
        if m["alpha_start"] <= pos <= m["alpha_end"]:
            return "alpha"
        if i + 1 < len(motifs) and m["alpha_end"] < pos < motifs[i + 1]["beta_start"]:
            return "inter"
    return "flanking"


# ---------------------------------------------------------------------------
# Resync search (identical to analyze_as_splice_junctions.py)
# ---------------------------------------------------------------------------

def find_resync(can_seq, can_end, iso_seq):
    offsets = [0] + [s for d in range(1, MAX_SLIDE + 1) for s in (d, -d)]
    for delta in offsets:
        pos = can_end + delta
        if pos < 0 or pos >= len(can_seq):
            continue
        window = min(MIN_MATCH, len(can_seq) - pos)
        if window < MIN_MATCH // 2:
            continue
        target  = can_seq[pos : pos + window]
        iso_idx = iso_seq.find(target)
        if iso_idx >= 0:
            return pos, iso_idx, delta
    return None


# ---------------------------------------------------------------------------
# Load canonical proteins (same filter as other analysis scripts)
# ---------------------------------------------------------------------------

def load_proteins(conn):
    """Return dict uid -> protein record. Excludes proteins with no Ensembl transcript."""
    enst_jcts = load_canonical_junctions(conn)
    rows = conn.execute("""
        SELECT uniprot_id, gene_name, domain_start, domain_end,
               motif_annotations
        FROM   view_canonical
    """).fetchall()
    proteins = {}
    for uid, gene, ds, de, ma in rows:
        if uid not in enst_jcts:
            continue                    # no Ensembl transcript — excluded
        motifs    = json.loads(ma)
        junctions = enst_jcts[uid]
        proteins[uid] = dict(uid=uid, gene=gene, ds=ds, de=de,
                             motifs=motifs, n_motifs=len(motifs),
                             junctions=junctions)
    return proteins


# ---------------------------------------------------------------------------
# Load AS isoforms
# ---------------------------------------------------------------------------

def load_as_isoforms(conn, proteins):
    """
    For each isoform with a matched Ensembl transcript, identify lost canonical
    junctions by direct sequence comparison (divergence point D_seq → resync
    R_can), consistent with analyze_as_splice_junctions.py.

    VSP can_end coordinates are used only as the starting point for the suffix
    match that locates R_can; they do not define which junctions are AS-affected.
    Only isoforms where at least one VSP resyncs successfully contribute.
    """
    iso_jcts = load_isoform_junctions(conn)   # isoform_id -> junction list (isoform coords)

    can_seqs = {r[0]: r[1] for r in conn.execute(
        "SELECT uniprot_id, sequence FROM canonical_analysis WHERE sequence IS NOT NULL"
    ).fetchall()}

    rows = conn.execute("""
        SELECT i.isoform_id, i.sequence, nc.canonical_id, nc.vsp_domain_events
        FROM   isoforms i
        JOIN   view_noncanonical nc ON nc.isoform_id = i.isoform_id
        WHERE  i.is_canonical = 0
          AND  i.sequence IS NOT NULL
          AND  nc.vsp_domain_events IS NOT NULL
    """).fetchall()

    isoforms = []
    for iso_id, iso_seq, can_id, vde in rows:
        if iso_id not in iso_jcts:
            continue                    # no Ensembl transcript — excluded
        if can_id not in proteins:
            continue
        if can_id not in can_seqs:
            continue

        can_seq = can_seqs[can_id]
        p = proteins[can_id]

        # Divergence point D_seq (1-indexed canonical position)
        diverge = None
        for i, (ca, ia) in enumerate(zip(can_seq, iso_seq), start=1):
            if ca != ia:
                diverge = i
                break
        if diverge is None:
            continue                    # identical or pure prefix — no AS event

        # Collect lost canonical junctions across all VSPs (union, deduplicated)
        vsps = json.loads(vde)
        as_junctions = set()
        for vsp in sorted(vsps, key=lambda v: v.get("can_start", 0)):
            can_end = vsp.get("can_end")
            if can_end is None:
                continue
            resync = find_resync(can_seq, can_end, iso_seq)
            if resync is None:
                continue
            R_can, _R_iso, _slide = resync
            if diverge >= R_can:
                continue
            for j in p["junctions"]:
                if diverge <= j < R_can:
                    as_junctions.add(j)

        as_junctions = sorted(as_junctions)
        if not as_junctions:
            continue
        isoforms.append(dict(uid=can_id, isoform_id=iso_id,
                             as_junctions=as_junctions,
                             n_as=len(as_junctions)))
    return isoforms


# ---------------------------------------------------------------------------
# 8A — Element-level distribution
# ---------------------------------------------------------------------------

def compute_8a(proteins, isoforms):
    """
    Returns N_t_AS, f_t_AS, rho_t_AS, baseline f_t, canonical N, AS N.
    """
    N_total = 0
    N_t = defaultdict(int)
    for p in proteins.values():
        for j in p["junctions"]:
            N_t[_tau5(j, p["motifs"])] += 1
            N_total += 1
    f_t = {t: N_t[t] / N_total for t in CATS5}

    N_AS_total = 0
    N_t_AS = defaultdict(int)
    for iso in isoforms:
        p = proteins[iso["uid"]]
        for j in iso["as_junctions"]:
            N_t_AS[_tau5(j, p["motifs"])] += 1
            N_AS_total += 1

    f_t_AS  = {t: N_t_AS[t] / N_AS_total for t in CATS5} if N_AS_total else {}
    rho_t_AS = {t: (f_t_AS[t] / f_t[t]) if f_t.get(t, 0) > 0 else float("nan")
                for t in CATS5}
    return N_t_AS, f_t_AS, rho_t_AS, f_t, N_total, N_AS_total


def chi_square_pvalues_9a(N_t_AS, f_t, N_AS_total):
    """Chi-square Pearson z-score test for §9A element enrichment."""
    pvals_raw = {}
    z_scores  = {}
    for t in CATS5:
        E_t = N_AS_total * f_t[t]
        if E_t > 0:
            z = (N_t_AS[t] - E_t) / np.sqrt(E_t)
        else:
            z = float("nan")
        z_scores[t]  = float(z)
        pvals_raw[t] = float(2 * stats.norm.sf(abs(z))) if not np.isnan(z) else 1.0
    pvals_bh = bh_correct(pvals_raw, CATS5)
    return pvals_raw, pvals_bh, z_scores


def plot_8a(rho_t_AS, f_t, N_AS_total, pvals_raw, pvals_bh, n_iso, N_AS, out):
    fig, ax = plt.subplots(figsize=(8, 4))
    x        = np.arange(len(CATS5))
    rho_vals = [rho_t_AS.get(t, 1.0) for t in CATS5]

    # BH-consistent CI: z_r / sqrt(E_t), lower bar clipped so CI doesn't go below 0
    alpha = 0.05
    K = len(CATS5)
    sorted_cats = sorted(CATS5, key=lambda t: pvals_raw[t])
    ranks = {t: i + 1 for i, t in enumerate(sorted_cats)}
    z_bh = {t: stats.norm.ppf(1 - alpha * ranks[t] / (2 * K)) for t in CATS5}
    hws = []
    for t in CATS5:
        E_t = N_AS_total * f_t[t]
        hws.append(z_bh[t] / np.sqrt(E_t) if E_t > 0 else 0.0)
    lo_errs = [min(rho, hw) for rho, hw in zip(rho_vals, hws)]
    hi_errs = hws

    ax.bar(x, rho_vals, width=0.65, color=[COLS5[t] for t in CATS5],
           alpha=0.85, zorder=3)
    ax.errorbar(x, rho_vals, yerr=[lo_errs, hi_errs],
                fmt="none", color="black", capsize=5, lw=1.2, zorder=4)
    ax.axhline(1.0, color="black", lw=0.8, ls="--", zorder=2)

    for i, t in enumerate(CATS5):
        p   = pvals_bh.get(t, float("nan"))
        sig = ("**" if p < 0.01 else ("*" if p < 0.05 else "")) if not np.isnan(p) else ""
        if sig:
            top = rho_vals[i] + hi_errs[i] + 0.04
            ax.text(x[i], top, sig, ha="center", va="bottom",
                    fontsize=11, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([LABELS5[t] for t in CATS5], fontsize=9)
    ax.set_ylabel(r"$\rho_t^{AS}$ (AS / canonical fraction)", fontsize=10)
    ax.set_title(
        f"AS junction element enrichment relative to canonical baseline\n"
        f"{n_iso} isoforms, $N^{{AS}}$ = {N_AS}  |  "
        "error bars = BH-adjusted CI;  ** BH $\\chi^2$ p < 0.01,  * BH $\\chi^2$ p < 0.05",
        fontsize=9,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# 8B — Cross-protein positional distribution
# ---------------------------------------------------------------------------

def normalized_positions_as(proteins, isoforms):
    xs = []
    for iso in isoforms:
        p = proteins[iso["uid"]]
        ds, de = p["ds"], p["de"]
        dlen = de - ds
        for j in iso["as_junctions"]:
            xs.append((j - ds) / dlen)
    return np.array(xs)


def normalized_positions_canonical(proteins):
    xs = []
    for p in proteins.values():
        ds, de = p["ds"], p["de"]
        dlen = de - ds
        for j in p["junctions"]:
            xs.append((j - ds) / dlen)
    return np.array(xs)




def mean_element_spans(proteins):
    spans = defaultdict(list)
    for p in proteins.values():
        ds, de = p["ds"], p["de"]
        dlen = de - ds
        if dlen <= 0:
            continue
        for i, m in enumerate(p["motifs"]):
            n = m["motif"]
            def r(x, ds=ds, dlen=dlen): return (x - ds) / dlen
            spans[f"beta_{n}"].append( (r(m["beta_start"]),      r(m["beta_end"])) )
            spans[f"loop_{n}"].append( (r(m["beta_end"] + 1),    r(m["alpha_start"] - 1)) )
            spans[f"alpha_{n}"].append((r(m["alpha_start"]),      r(m["alpha_end"])) )
            if i + 1 < len(p["motifs"]):
                nxt = p["motifs"][i + 1]
                spans[f"inter_{n}_{n+1}"].append(
                    (r(m["alpha_end"] + 1), r(nxt["beta_start"] - 1)))
    return {k: (np.mean([s for s, e in v]), np.mean([e for s, e in v]))
            for k, v in spans.items() if v}


def plot_8b(x_as, x_const, D_obs, ks_p, spans, n_iso, out):
    bw   = "scott"
    grid = np.linspace(0, 1, 500)
    BG   = {"beta":  ("#4C72B0", 0.18), "alpha": ("#DD8452", 0.18),
            "loop":  ("#eeeeee", 0.80), "inter": ("#eeeeee", 0.80)}

    kde_as    = stats.gaussian_kde(x_as,    bw_method=bw)(grid) if len(x_as) > 1 else np.ones(500)
    kde_const = stats.gaussian_kde(x_const, bw_method=bw)(grid)

    fig, ax = plt.subplots(figsize=(10, 4.5))
    for elem, (s, e) in spans.items():
        if e <= s:
            continue
        prefix = elem.split("_")[0]
        col, alpha = BG.get(prefix, ("#eeeeee", 0.5))
        ax.axvspan(s, e, color=col, alpha=alpha, zorder=0)

    ax.plot(grid, kde_const, color="#4C72B0", lw=1.2, ls="--", zorder=3,
            label=f"Canonical junctions (N = {len(x_const)})")
    ax.plot(grid, kde_as,    color="black",   lw=1.8, zorder=4,
            label=f"AS-affected junctions (N = {len(x_as)})")

    bg_patches = [
        mpatches.Patch(facecolor="#4C72B0", alpha=0.40, label="β-strand region (mean)"),
        mpatches.Patch(facecolor="#DD8452", alpha=0.40, label="α-helix region (mean)"),
        mpatches.Patch(facecolor="#cccccc", alpha=0.80, label="Loop / inter-motif (mean)"),
    ]
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles + bg_patches, fontsize=8, frameon=False,
              loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3)

    ax.set_xlabel("Normalised domain position", fontsize=10)
    ax.set_ylabel("Density", fontsize=10)
    ax.set_xlim(0, 1)
    ax.set_title(
        f"Positional distribution: AS-affected vs canonical junctions\n"
        f"{n_iso} isoforms  |  KS D = {D_obs:.4f},  p = {ks_p:.4f}  (analytical)",
        fontsize=9,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# 8C — Within-protein hotspot reuse
# ---------------------------------------------------------------------------

def compute_8c(proteins, isoforms):
    """
    Returns H_bar, multi-isoform protein dict, per-protein stats,
    total hotspot count, total unique junction count.
    """
    by_protein = defaultdict(list)
    for iso in isoforms:
        by_protein[iso["uid"]].append(iso["as_junctions"])

    multi = {uid: jsets for uid, jsets in by_protein.items() if len(jsets) >= 2}

    total_hotspot = 0
    total_unique  = 0
    per_protein   = {}
    for uid, jsets in multi.items():
        usage = defaultdict(int)
        for jlist in jsets:
            for j in jlist:
                usage[j] += 1
        union_j   = set(usage)
        hotspot_j = {j for j, u in usage.items() if u >= 2}
        H_p = len(hotspot_j) / len(union_j) if union_j else 0.0
        per_protein[uid] = dict(H_p=H_p, hotspot=len(hotspot_j),
                                unique=len(union_j), k=len(jsets))
        total_hotspot += len(hotspot_j)
        total_unique  += len(union_j)

    H_bar = total_hotspot / total_unique if total_unique > 0 else 0.0
    return H_bar, multi, per_protein, total_hotspot, total_unique


def hotspot_z_test(proteins, multi):
    """
    Compare observed total hotspot junctions against the analytically expected
    count under Bernoulli independence (with-replacement approximation).

    For each protein p with isoforms drawing m_1,...,m_k junctions from pool n:
      P(junction j shared by >=2 isoforms) = 1 - P(in 0) - P(in exactly 1)
    Expected hotspot count = n * P(hotspot).

    Returns (observed, expected, z, p_one_sided).
    """
    total_obs = 0
    total_exp = 0.0
    for uid, jsets in multi.items():
        n = len(proteins[uid]["junctions"])
        if n == 0:
            continue
        m_list = [len(set(js)) for js in jsets]
        p_in   = [m / n for m in m_list]
        p_none = float(np.prod([1 - p for p in p_in]))
        p_one  = sum(
            p_in[i] * float(np.prod([1 - p_in[j] for j in range(len(p_in)) if j != i]))
            for i in range(len(p_in))
        )
        p_hotspot = max(0.0, 1 - p_none - p_one)

        usage = defaultdict(int)
        for jlist in jsets:
            for j in jlist:
                usage[j] += 1
        total_obs += sum(1 for u in usage.values() if u >= 2)
        total_exp += n * p_hotspot

    z = (total_obs - total_exp) / np.sqrt(total_exp) if total_exp > 0 else 0.0
    p = float(stats.norm.sf(z))
    return total_obs, total_exp, float(z), p


def plot_8c(proteins, multi, total_obs, total_exp, p_val, n_multi, total_unique, out):
    obs_per_prot = []
    exp_per_prot = []
    for uid, jsets in multi.items():
        n = len(proteins[uid]["junctions"])
        if n == 0:
            continue
        m_list = [len(set(js)) for js in jsets]
        p_in   = [m / n for m in m_list]
        p_none = float(np.prod([1 - p for p in p_in]))
        p_one  = sum(
            p_in[i] * float(np.prod([1 - p_in[j] for j in range(len(p_in)) if j != i]))
            for i in range(len(p_in))
        )
        p_hotspot = max(0.0, 1 - p_none - p_one)

        usage = defaultdict(int)
        for jlist in jsets:
            for j in jlist:
                usage[j] += 1
        obs_per_prot.append(sum(1 for u in usage.values() if u >= 2))
        exp_per_prot.append(n * p_hotspot)

    lim = max(max(obs_per_prot, default=1), max(exp_per_prot, default=1)) + 0.5
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(exp_per_prot, obs_per_prot, color="#4C72B0",
               alpha=0.75, edgecolors="white", s=50, zorder=3)
    ax.plot([0, lim], [0, lim], color="black", lw=0.9,
            ls="--", zorder=2, label="Expected = Observed")
    ax.set_xlabel("Expected hotspot junctions (Bernoulli null)", fontsize=10)
    ax.set_ylabel("Observed hotspot junctions", fontsize=10)
    ax.set_title(
        f"Within-protein AS hotspot reuse  ({n_multi} proteins, k ≥ 2)\n"
        f"Total: {total_obs} observed / {total_exp:.1f} expected  |  "
        f"z-test p = {p_val:.4f}",
        fontsize=9,
    )
    ax.legend(fontsize=8, frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# BH correction
# ---------------------------------------------------------------------------

def bh_correct(pvals_dict, keys):
    pvals = np.array([pvals_dict[k] for k in keys])
    order = np.argsort(pvals)
    m = len(pvals)
    adj = np.zeros(m)
    for rank, idx in enumerate(order, 1):
        adj[idx] = pvals[idx] * m / rank
    for j in range(m - 2, -1, -1):
        adj[order[j]] = min(adj[order[j]], adj[order[j + 1]])
    adj = np.minimum(adj, 1.0)
    return {k: float(adj[i]) for i, k in enumerate(keys)}


# ---------------------------------------------------------------------------
# Update §8 in Statistical-Analysis.md
# ---------------------------------------------------------------------------

def _replace_block(text, start_tag, end_tag, new_content):
    si = text.find(start_tag)
    if si == -1:
        return text
    content_start = si + len(start_tag)
    if end_tag is None:
        return text[:content_start] + "\n" + new_content + "\n"
    ei = text.find(end_tag, content_start)
    if ei == -1:
        return text[:content_start] + "\n" + new_content + "\n"
    return text[:content_start] + "\n" + new_content + "\n" + text[ei:]


def update_section9(md_path, content_9a, content_9b, content_9c):
    text = Path(md_path).read_text(encoding="utf-8")
    text = _replace_block(text, "#### 9A Results\n", "\n---\n\n### 9B.", content_9a)
    text = _replace_block(text, "#### 9B Results\n", "\n---\n\n### 9C.", content_9b)
    text = _replace_block(text, "#### 9C Results\n", "\n---\n\n### 9D.", content_9c)
    Path(md_path).write_text(text, encoding="utf-8")
    print(f"Updated §9 in {md_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db",        default=None)
    parser.add_argument("--B",         type=int, default=2000)
    parser.add_argument("--seed",      type=int, default=42)
    parser.add_argument("--out-8a",    default="figures/as_element_enrichment.png")
    parser.add_argument("--out-8b",    default="figures/as_positional.png")
    parser.add_argument("--out-8c",    default="figures/as_hotspot.png")
    parser.add_argument("--md",        default="Statistical-Analysis.md")
    args = parser.parse_args()

    db_path  = args.db or get_config().db_path
    conn     = sqlite3.connect(db_path)
    proteins = load_proteins(conn)
    isoforms = load_as_isoforms(conn, proteins)
    conn.close()

    n_prot  = len(proteins)
    n_iso   = len(isoforms)
    N_AS    = sum(iso["n_as"] for iso in isoforms)
    n_genes = len(set(iso["uid"] for iso in isoforms))
    print(f"Loaded {n_prot} canonical proteins.")
    print(f"AS isoforms with >=1 domain junction in VSP span: {n_iso} "
          f"({n_genes} distinct canonical proteins), {N_AS} AS-affected junction instances.")

    if n_iso == 0:
        print("No AS isoforms found -- nothing to analyse.")
        return

    # ── 8A. Element distribution ─────────────────────────────────────────────
    print("\n8A. Element enrichment (chi-square) ...")
    N_t_AS, f_t_AS, rho_t_AS, f_t, N_total, N_AS_total = compute_8a(proteins, isoforms)
    for t in CATS5:
        print(f"   {ALABELS5[t]}: N_t_AS={N_t_AS[t]}, "
              f"f_t_AS={f_t_AS.get(t,0):.3f}, f_t={f_t[t]:.3f}, "
              f"rho={rho_t_AS.get(t,float('nan')):.3f}")

    pvals_8a, pvals_8a_bh, z_scores_8a = chi_square_pvalues_9a(N_t_AS, f_t, N_AS_total)
    for t in CATS5:
        print(f"   {ALABELS5[t]}: z={z_scores_8a[t]:.3f}, raw p={pvals_8a[t]:.4f}, "
              f"BH p={pvals_8a_bh[t]:.4f}")

    plot_8a(rho_t_AS, f_t, N_AS_total, pvals_8a, pvals_8a_bh, n_iso, N_AS_total, args.out_8a)

    # ── 8B. Positional distribution ──────────────────────────────────────────
    print("\n8B. Positional distribution ...")
    x_as    = normalized_positions_as(proteins, isoforms)
    x_const = normalized_positions_canonical(proteins)
    D_obs, p_8b = stats.ks_2samp(x_as, x_const)
    print(f"   KS D = {D_obs:.4f}, p = {p_8b:.4f}")

    spans = mean_element_spans(proteins)
    plot_8b(x_as, x_const, D_obs, p_8b, spans, n_iso, args.out_8b)

    # ── 8C. Hotspot reuse ────────────────────────────────────────────────────
    print(f"\n8C. Hotspot reuse (B = {args.B}) ...")
    H_bar, multi, per_protein, total_hotspot, total_unique = compute_8c(proteins, isoforms)
    n_multi = len(multi)
    print(f"   Multi-isoform proteins (k>=2): {n_multi}")
    print(f"   H_bar = {H_bar:.4f}  ({total_hotspot} hotspot / {total_unique} unique junctions)")

    obs_8c, exp_8c, z_8c, p_8c = hotspot_z_test(proteins, multi)
    print(f"   Expected hotspots = {exp_8c:.2f}  |  z = {z_8c:.3f}  |  p = {p_8c:.4f}")

    plot_8c(proteins, multi, obs_8c, exp_8c, p_8c, n_multi, total_unique, args.out_8c)

    # ── Build markdown content ────────────────────────────────────────────────
    def sig(p):
        return "**" if p < 0.01 else ("*" if p < 0.05 else "ns")

    rows_8a = [
        f"| {LABELS5[t]} | {N_t_AS[t]} | {f_t_AS.get(t,0):.3f} | {f_t[t]:.3f} | "
        f"{rho_t_AS.get(t, float('nan')):.3f} | "
        f"{z_scores_8a[t]:.3f} | "
        f"{pvals_8a[t]:.4f} | {pvals_8a_bh[t]:.4f} | {sig(pvals_8a_bh[t])} |"
        for t in CATS5
    ]

    content_8a_md = (
        f"Dataset: $|\\mathcal{{A}}|$ = {n_iso} isoforms across {n_genes} canonical proteins; "
        f"$N^{{AS}}$ = {N_AS_total} AS-affected junction instances.\n\n"
        "| Element | $N_t^{AS}$ | $f_t^{AS}$ | $f_t$ | $\\rho_t^{AS}$ | "
        "$\\chi^2$ $z$ | Raw $p$ | BH $p$ | Sig |\n"
        "|---|---|---|---|---|---|---|---|---|\n"
        + "\n".join(rows_8a)
        + "\n\nPearson $\\chi^2$ $z$-scores are two-sided; BH correction across all 5 element types; "
        "error bars in figure show 95% Poisson CI.\n\n"
        "![AS junction element enrichment](figures/as_element_enrichment.png)"
    )

    reject_8b = p_8b < 0.05
    content_8b_md = (
        f"$N^{{AS}}$ = {len(x_as)} normalised positions from {n_iso} isoforms; "
        f"canonical baseline $N$ = {len(x_const)}.\n\n"
        "| Statistic | Value |\n"
        "|---|---|\n"
        f"| KS statistic $D$ | {D_obs:.4f} |\n"
        f"| Analytical $p$ | {p_8b:.4f} |\n\n"
        f"$H_0^{{AS}}$ **{'rejected' if reject_8b else 'not rejected'}** "
        f"(analytical $p = {p_8b:.4f}$, $\\alpha = 0.05$). "
        + (
            "The AS-affected junction positional distribution is significantly different "
            "from the canonical baseline."
            if reject_8b else
            "The AS-affected junction positional distribution is not significantly different "
            "from the canonical baseline; AS-affected canonical junctions do not cluster at "
            "different domain positions than the full canonical pool."
        )
        + "\n\n![AS junction positional distribution](figures/as_positional.png)"
    )

    reject_8c = p_8c < 0.05
    content_8c_md = (
        f"Multi-isoform proteins ($k_p \\ge 2$): $|\\mathcal{{P}}^{{(2)}}|$ = {n_multi}.\n\n"
        "| Statistic | Value |\n"
        "|---|---|\n"
        f"| Proteins with $k_p \\ge 2$ | {n_multi} |\n"
        f"| Total unique AS-affected junctions | {total_unique} |\n"
        f"| Hotspot junctions ($u_p(j) \\ge 2$) | {total_hotspot} |\n"
        f"| Observed $\\bar{{H}}$ | {H_bar:.4f} |\n"
        f"| Expected hotspots (Bernoulli null) | {exp_8c:.2f} |\n"
        f"| $z$ | {z_8c:.3f} |\n"
        f"| $p$ (one-sided) | {p_8c:.4f} |\n\n"
        f"$H_0^{{AS}}$ **{'rejected' if reject_8c else 'not rejected'}** "
        f"($z = {z_8c:.3f}$, $p = {p_8c:.4f}$, $\\alpha = 0.05$). "
        + (
            "The observed hotspot count significantly exceeds the Bernoulli expectation, "
            "indicating that multiple AS isoforms of the same protein non-randomly share canonical "
            "junctions — consistent with within-protein AS hotspots."
            if reject_8c else
            "The observed hotspot count is consistent with the Bernoulli null; "
            "within-protein AS hotspot reuse is not significantly above chance."
        )
        + "\n\n![Within-protein hotspot reuse](figures/as_hotspot.png)"
    )

    md_path = Path(args.md)
    if md_path.exists():
        update_section9(str(md_path), content_8a_md, content_8b_md, content_8c_md)
    else:
        print(f"  [note] {md_path} not found — markdown not updated.")


if __name__ == "__main__":
    main()
