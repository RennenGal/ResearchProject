#!/usr/bin/env python3
"""
AS splice junction enrichment in TIM-barrel structural elements.

For each VSP event in each non-canonical isoform, identifies the actual
canonical exon junctions that flank the alternatively spliced region:

  Entry junction: last canonical junction strictly before VSP can_start.
  Exit junction:  first canonical junction at or after the resync point,
                  where the resync is found by suffix-matching
                  canonical[can_end : can_end+MIN_MATCH] in the isoform
                  sequence (±MAX_SLIDE residues to absorb annotation noise).

Isoforms are matched to Ensembl transcripts by protein sequence; unmatched
ones are retried by stored ENST ID.  Truncations (isoform ends before
can_end) are flagged and receive no exit junction.  Multiple VSP events in
the same isoform are processed independently; if the shared segment between
two events is shorter than MIN_MATCH the second entry is flagged.

Domain-internal entry and exit junctions are pooled and tested for
enrichment against a length-weighted null using the same chi-square
framework as analyze_junction_enrichment.py.

Output:
  figures/as_splice_junctions.png

Usage:
    python scripts/analyze_as_splice_junctions.py
"""

import argparse
import json
import sqlite3
import sys

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent.parent))
from protein_data_collector.config import get_config


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CATS5   = ["beta", "alpha", "inter", "loop", "flanking"]
LABELS5 = {"beta": "β-strand", "alpha": "α-helix", "inter": "Inter-motif",
           "loop": "Loop (β→α)", "flanking": "Flanking"}
ALABELS = {"beta": "beta-strand", "alpha": "alpha-helix", "inter": "Inter-motif",
           "loop": "Loop (b->a)", "flanking": "Flanking"}
COLS5   = {"beta": "#4C72B0", "alpha": "#DD8452", "inter": "#55A868",
           "loop": "#8c8c8c", "flanking": "#C44E52"}

MIN_MATCH  = 15   # minimum AA window for suffix match to be trusted
MAX_SLIDE  = 5    # max residues to slide can_end when looking for resync


# ---------------------------------------------------------------------------
# Structural element classification
# ---------------------------------------------------------------------------

def _tau5(pos, motifs):
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


def build_type_array(ds, de, motifs):
    cat_idx = {c: i for i, c in enumerate(CATS5)}
    arr = np.zeros(de - ds, dtype=np.int8)
    for idx in range(de - ds):
        arr[idx] = cat_idx[_tau5(ds + idx, motifs)]
    return arr


def compute_residue_null(canonicals):
    """
    Background distribution of structural element types across ALL domain
    residues (length-weighted, same null as Analysis 1 and Ochoa-Leyva).
    """
    counts = {c: 0 for c in CATS5}
    total  = 0
    for can in canonicals.values():
        ds, de = can["ds"], can["de"]
        if de <= ds:
            continue
        arr = build_type_array(ds, de, can["motifs"])
        for v in arr:
            counts[CATS5[v]] += 1
            total += 1
    if total == 0:
        return {c: 1 / len(CATS5) for c in CATS5}
    return {c: counts[c] / total for c in CATS5}


def compute_junction_null(canonicals):
    """
    Background distribution of structural element types at ALL domain-internal
    canonical junction positions (used as null for AS-junction enrichment).
    """
    counts = {c: 0 for c in CATS5}
    total  = 0
    for can in canonicals.values():
        ds, de = can["ds"], can["de"]
        if de <= ds:
            continue
        arr = build_type_array(ds, de, can["motifs"])
        for j in can["junctions"]:
            if ds <= j < de:
                counts[CATS5[arr[j - ds]]] += 1
                total += 1
    if total == 0:
        return {c: 1 / len(CATS5) for c in CATS5}
    return {c: counts[c] / total for c in CATS5}


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_canonicals(conn):
    rows = conn.execute("""
        SELECT i.uniprot_id, i.sequence,
               vc.domain_start, vc.domain_end,
               vc.exon_annotations, vc.motif_annotations
        FROM   isoforms i
        JOIN   view_canonical vc ON vc.uniprot_id = i.uniprot_id
        WHERE  i.is_canonical = 1 AND i.sequence IS NOT NULL
    """).fetchall()
    out = {}
    for uid, seq, ds, de, ea, ma in rows:
        exons     = json.loads(ea)
        motifs    = json.loads(ma)
        junctions = sorted(e["end"] for e in exons[:-1])
        out[uid]  = dict(seq=seq, ds=ds, de=de,
                         junctions=junctions, motifs=motifs)
    return out


def load_enst_lookup(conn):
    """Build two lookups: by protein sequence and by bare ENST ID."""
    by_seq = {}
    by_id  = {}
    for enst_id, seq, ea in conn.execute(
        "SELECT enst_id, sequence, exon_annotations "
        "FROM   ensembl_transcripts WHERE exon_annotations IS NOT NULL"
    ):
        base   = enst_id.split('.')[0]
        parsed = json.loads(ea)
        by_seq[seq]  = (enst_id, parsed)
        by_id[base]  = (enst_id, parsed)
    return by_seq, by_id


def load_isoforms(conn):
    rows = conn.execute("""
        SELECT i.isoform_id, i.uniprot_id, i.sequence,
               i.ensembl_transcript_id,
               nc.vsp_domain_events
        FROM   isoforms i
        JOIN   view_noncanonical nc ON nc.isoform_id = i.isoform_id
        WHERE  i.is_canonical = 0 AND i.sequence IS NOT NULL
    """).fetchall()
    return [
        dict(iso_id=r[0], uid=r[1], seq=r[2],
             enst_hint=r[3], vsps=json.loads(r[4]))
        for r in rows
    ]


# ---------------------------------------------------------------------------
# Transcript matching
# ---------------------------------------------------------------------------

def match_transcript(iso_seq, enst_hint, by_seq, by_id):
    if iso_seq in by_seq:
        return by_seq[iso_seq]
    if enst_hint:
        base = enst_hint.split('.')[0]
        if base in by_id:
            return by_id[base]
    return None


# ---------------------------------------------------------------------------
# Suffix match for exit junction
# ---------------------------------------------------------------------------

def find_resync(can_seq, can_end, iso_seq):
    """
    Confirm that canonical and isoform sequences rejoin at can_end by
    searching for can_seq[can_end : can_end+MIN_MATCH] in iso_seq.
    Tries can_end exactly, then ±1 … ±MAX_SLIDE.

    Returns (resync_pos, match_len, slide_offset) or None.
    """
    offsets = [0] + [s for d in range(1, MAX_SLIDE + 1) for s in (d, -d)]
    for delta in offsets:
        pos = can_end + delta
        if pos < 0 or pos >= len(can_seq):
            continue
        window = min(MIN_MATCH, len(can_seq) - pos)
        if window < MIN_MATCH // 2:
            continue
        target = can_seq[pos : pos + window]
        if iso_seq.find(target) >= 0:
            return pos, window, delta
    return None


# ---------------------------------------------------------------------------
# Core: find entry and exit junctions per VSP event
# ---------------------------------------------------------------------------

def find_as_junctions(canonicals, isoforms, by_seq, by_id):
    """
    Returns one record per (isoform, VSP-event) pair:
      iso_id, uid, feat_id,
      can_start, can_end,
      entry, entry_in_domain,
      exit,  exit_in_domain,
      is_truncation, resync_slide, flags
    """
    records = []

    for iso in isoforms:
        uid, iso_id, iso_seq = iso["uid"], iso["iso_id"], iso["seq"]

        if uid not in canonicals:
            continue
        if match_transcript(iso_seq, iso["enst_hint"], by_seq, by_id) is None:
            continue

        can       = canonicals[uid]
        can_seq   = can["seq"]
        can_junc  = can["junctions"]
        ds, de    = can["ds"], can["de"]

        # Compute per-isoform sequence divergence point once.
        # diverge = first AA (1-indexed) where sequences differ.
        # If isoform is a pure prefix (no mismatch before isoform ends),
        # set diverge = len(iso_seq) + 1 as a truncation sentinel.
        diverge = None
        for i, (ca, ia) in enumerate(zip(can_seq, iso_seq), start=1):
            if ca != ia:
                diverge = i
                break
        if diverge is None:
            diverge = len(iso_seq) + 1
        isoform_is_prefix = (diverge == len(iso_seq) + 1)

        # Process VSPs in genomic order
        for vsp in sorted(iso["vsps"], key=lambda v: v.get("can_start", 0)):
            can_start = vsp.get("can_start")
            can_end   = vsp.get("can_end")
            feat_id   = vsp.get("feature_id", "")
            if can_start is None or can_end is None:
                continue

            flags = []

            # ── Entry: last canonical junction before can_start ───────────
            before = [j for j in can_junc if j < can_start]
            entry  = max(before) if before else None

            # ── Exit via suffix match ─────────────────────────────────────
            # Always attempt the suffix match first, regardless of isoform
            # length vs can_end. Comparing lengths mixes coordinate systems
            # (isoform vs canonical) and misclassifies exon-skip isoforms.
            # Only classify as truncation if the match fails AND the isoform
            # is a pure prefix of the canonical (no alternative sequence).
            exit_j       = None
            resync_slide = None
            is_trunc     = False

            resync = find_resync(can_seq, can_end, iso_seq)
            if resync is not None:
                resync_pos, match_len, slide = resync
                resync_slide = slide
                if slide != 0:
                    flags.append(f"resync_slid{slide:+d}")
                if match_len < MIN_MATCH:
                    flags.append(f"short_match_{match_len}aa")
                after  = [j for j in can_junc if j >= resync_pos]
                exit_j = min(after) if after else None
            elif isoform_is_prefix and len(iso_seq) <= can_end:
                is_trunc = True
                flags.append("truncation")
            else:
                flags.append("no_resync")

            records.append(dict(
                iso_id           = iso_id,
                uid              = uid,
                feat_id          = feat_id,
                can_start        = can_start,
                can_end          = can_end,
                entry            = entry,
                entry_in_domain  = (entry  is not None and ds <= entry  < de),
                exit             = exit_j,
                exit_in_domain   = (exit_j is not None and ds <= exit_j < de),
                is_truncation    = is_trunc,
                resync_slide     = resync_slide,
                flags            = flags,
            ))

    return records


# ---------------------------------------------------------------------------
# Enrichment analysis
# ---------------------------------------------------------------------------

def compute_enrichment(junction_positions, canonicals, pi_null):
    """
    junction_positions: list of (uid, aa_position) tuples (domain-internal).
    pi_null: background element-type fractions (from compute_junction_null).
    Returns dict with N, N_t, pi_t, f_t, rho_t.
    """
    N_t = {c: 0 for c in CATS5}
    N   = 0

    for uid, pos in junction_positions:
        if uid not in canonicals:
            continue
        can    = canonicals[uid]
        ds, de = can["ds"], can["de"]
        if de <= ds:
            continue
        arr = build_type_array(ds, de, can["motifs"])
        N_t[CATS5[arr[pos - ds]]] += 1
        N += 1

    if N == 0:
        return None

    f_t   = {c: N_t[c] / N for c in CATS5}
    rho_t = {c: f_t[c] / pi_null[c] if pi_null[c] > 0 else 0.0 for c in CATS5}
    return dict(N=N, N_t=N_t, pi_t=pi_null, f_t=f_t, rho_t=rho_t)


def chi_square_pvalues(result):
    N, pvals_raw = result["N"], {}
    for c in CATS5:
        E_t = N * result["pi_t"][c]
        if E_t > 0:
            z = (result["N_t"][c] - E_t) / np.sqrt(E_t)
            pvals_raw[c] = float(2 * stats.norm.sf(abs(z)))
        else:
            pvals_raw[c] = 1.0
    return pvals_raw, bh_correction(pvals_raw)


def bh_correction(pvals_dict):
    cats  = list(pvals_dict.keys())
    pvals = np.array([pvals_dict[c] for c in cats])
    n     = len(pvals)
    order = np.argsort(pvals)
    adj   = np.zeros(n)
    for rank, idx in enumerate(order, 1):
        adj[idx] = pvals[idx] * n / rank
    for j in range(n - 2, -1, -1):
        adj[order[j]] = min(adj[order[j]], adj[order[j + 1]])
    return {c: float(min(adj[i], 1.0)) for i, c in enumerate(cats)}


# ---------------------------------------------------------------------------
# Console output
# ---------------------------------------------------------------------------

def sig_stars(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"


def print_enrichment(label, result, pvals_raw, pvals_bh, chi2, p_global, dof):
    N = result["N"]
    print(f"\n{'='*72}")
    print(f"  {label}  (N = {N})")
    print(f"  Global chi2({dof}) = {chi2:.2f},  p = {p_global:.4g}")
    print(f"{'='*72}")
    print(f"  {'Element':<16}  {'N_t':>5}  {'f_t':>7}  {'pi_t':>7}  "
          f"{'rho':>6}  {'p_raw':>7}  {'p_BH':>7}  Sig")
    print("  " + "-"*16 + "  " + "  ".join(["-"*5] + ["-"*7]*5) + "  ---")
    for c in CATS5:
        print(f"  {ALABELS[c]:<16}  {result['N_t'][c]:>5}  "
              f"{result['f_t'][c]:>7.4f}  {result['pi_t'][c]:>7.4f}  "
              f"{result['rho_t'][c]:>6.3f}  "
              f"{pvals_raw[c]:>7.4f}  {pvals_bh[c]:>7.4f}  {sig_stars(pvals_bh[c])}")
    print(f"{'='*72}")


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def poisson_ci(k, E_t, alpha=0.05):
    """Exact 95% Poisson CI for rate rho = k / E_t (chi-squared method)."""
    from scipy.stats import chi2 as _chi2
    lo = _chi2.ppf(alpha / 2,     2 * k)       / (2 * E_t) if k > 0 else 0.0
    hi = _chi2.ppf(1 - alpha / 2, 2 * (k + 1)) / (2 * E_t)
    return lo, hi


def plot_enrichment_panel(ax, result, pvals_raw, title):
    """Draw one enrichment panel onto ax using raw p-values for stars."""
    x    = np.arange(len(CATS5))
    rhos = [result["rho_t"][c] for c in CATS5]
    N    = result["N"]

    lo_errs, hi_errs = [], []
    for c in CATS5:
        k   = result["N_t"][c]
        E_t = N * result["pi_t"][c]
        lo, hi = poisson_ci(k, E_t) if E_t > 0 else (0.0, 0.0)
        rho = result["rho_t"][c]
        lo_errs.append(max(rho - lo, 0.0))
        hi_errs.append(max(hi - rho, 0.0))

    ax.bar(x, rhos, color=[COLS5[c] for c in CATS5],
           alpha=0.85, width=0.6, zorder=3)
    ax.axhline(1.0, color="black", lw=0.9, ls="--", zorder=2)
    ax.errorbar(x, rhos, yerr=[lo_errs, hi_errs], fmt="none",
                color="#888888", capsize=5, lw=1.2, zorder=4)

    for i, c in enumerate(CATS5):
        s = sig_stars(pvals_raw[c])
        if s != "ns":
            ax.text(i, rhos[i] + hi_errs[i] + 0.03, s,
                    ha="center", va="bottom", fontsize=10, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([LABELS5[c] for c in CATS5],
                       rotation=20, ha="right", fontsize=9)
    ax.set_ylabel("Enrichment ratio $\\rho_t$", fontsize=10)
    ax.set_title(title, fontsize=9)
    top = max(r + e for r, e in zip(rhos, hi_errs))
    ax.set_ylim(0, max(top + 0.3, 1.85))
    ax.yaxis.grid(True, linestyle=":", alpha=0.4, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_enrichment(res_entry, pv_raw_entry, res_exit, pv_raw_exit, out,
                    null_label="all canonical junctions"):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)

    plot_enrichment_panel(
        axes[0], res_entry, pv_raw_entry,
        f"Entry junctions  ($N$ = {res_entry['N']})",
    )
    plot_enrichment_panel(
        axes[1], res_exit, pv_raw_exit,
        f"Exit junctions  ($N$ = {res_exit['N']})",
    )
    axes[1].set_ylabel("")

    fig.suptitle(
        f"AS splice junction enrichment vs. {null_label}\n"
        "error bars = 95% CI;  *, **, *** = raw $p$ < 0.05 / 0.01 / 0.001  "
        "(BH-corrected $p$ in table)",
        fontsize=9,
    )
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db",  default=None)
    parser.add_argument("--out",         default="figures/as_splice_junctions.png")
    parser.add_argument("--out-residue", default="figures/as_splice_junctions_residue_null.png")
    args = parser.parse_args()

    db_path = args.db or get_config().db_path
    conn    = sqlite3.connect(db_path)

    print("Loading data ...")
    canonicals      = load_canonicals(conn)
    by_seq, by_id   = load_enst_lookup(conn)
    isoforms        = load_isoforms(conn)
    conn.close()

    print(f"  Canonical proteins:  {len(canonicals)}")
    print(f"  Non-canonical isoforms: {len(isoforms)}")

    pi_null    = compute_junction_null(canonicals)
    pi_residue = compute_residue_null(canonicals)
    print("\nJunction null (all canonical domain-internal junctions):")
    for c in CATS5:
        print(f"  {ALABELS[c]:<16}  {pi_null[c]:.4f}")
    print("\nResidue null (all canonical domain residues):")
    for c in CATS5:
        print(f"  {ALABELS[c]:<16}  {pi_residue[c]:.4f}")

    print("\nFinding AS splice junctions ...")
    records = find_as_junctions(canonicals, isoforms, by_seq, by_id)

    # ── Diagnostics ──────────────────────────────────────────────────────────
    n_events   = len(records)
    n_trunc    = sum(1 for r in records if r["is_truncation"])
    n_no_rsync = sum(1 for r in records if "no_resync" in r["flags"])
    n_slid     = sum(1 for r in records
                     if r["resync_slide"] is not None and r["resync_slide"] != 0)
    print(f"\n  VSP events processed:          {n_events}")
    print(f"    Truncations (no exit):        {n_trunc}")
    print(f"    No resync found:              {n_no_rsync}")
    print(f"    Resync with annotation slide: {n_slid}")

    entry_dom = [(r["uid"], r["entry"]) for r in records if r["entry_in_domain"]]
    exit_dom  = [(r["uid"], r["exit"])  for r in records if r["exit_in_domain"]]
    pooled    = entry_dom + exit_dom

    print(f"\n  Entry junctions in domain: {len(entry_dom)}")
    print(f"  Exit  junctions in domain: {len(exit_dom)}")
    print(f"  Pooled total:              {len(pooled)}")

    # ── Enrichment helper ────────────────────────────────────────────────────
    def run(positions, label, null):
        res = compute_enrichment(positions, canonicals, null)
        if res is None:
            print(f"\n  {label}: no domain-internal junctions.")
            return None
        pv_raw, pv_bh = chi_square_pvalues(res)
        obs = np.array([res["N_t"][c] for c in CATS5], dtype=float)
        exp = np.array([res["N"] * res["pi_t"][c] for c in CATS5], dtype=float)
        chi2, p_gl = stats.chisquare(obs, f_exp=exp)
        print_enrichment(label, res, pv_raw, pv_bh, chi2, p_gl, len(CATS5) - 1)
        return res, pv_raw, pv_bh, chi2, p_gl

    # ── Junction-null comparison ─────────────────────────────────────────────
    print("\n" + "="*72)
    print("  JUNCTION NULL")
    print("="*72)
    en_jn = run(entry_dom, "Entry junctions only",   pi_null)
    ex_jn = run(exit_dom,  "Exit  junctions only",   pi_null)
    run(pooled,     "Pooled (entry + exit)",          pi_null)

    if en_jn and ex_jn:
        res_en, pv_raw_en, _, _, _ = en_jn
        res_ex, pv_raw_ex, _, _, _ = ex_jn
        plot_enrichment(res_en, pv_raw_en, res_ex, pv_raw_ex, args.out)

    # ── Residue-null comparison (Ochoa-Leyva equivalent) ────────────────────
    print("\n" + "="*72)
    print("  RESIDUE NULL  (comparable to Ochoa-Leyva)")
    print("="*72)
    en_rn = run(entry_dom, "Entry junctions only [residue null]", pi_residue)
    ex_rn = run(exit_dom,  "Exit  junctions only [residue null]", pi_residue)
    run(pooled,     "Pooled (entry + exit) [residue null]",       pi_residue)

    if en_rn and ex_rn:
        res_en_r, pv_raw_en_r, _, _, _ = en_rn
        res_ex_r, pv_raw_ex_r, _, _, _ = ex_rn
        plot_enrichment(res_en_r, pv_raw_en_r, res_ex_r, pv_raw_ex_r,
                        args.out_residue,
                        null_label="all domain residues")


if __name__ == "__main__":
    main()
