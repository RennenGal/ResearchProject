#!/usr/bin/env python3
"""
Transcript-derived AS boundary enrichment in TIM-barrel structural elements.

For each isoform, the boundaries of the alternatively spliced region are
located by direct protein sequence comparison between the canonical and
isoform sequences:

  Divergence start D_seq : first residue where canonical and isoform
                           sequences differ (= start of the AS region in
                           canonical coordinates).

  Resync end R_can       : first canonical residue where the sequences
                           rejoin, found by suffix-matching 15 residues of
                           canonical[can_end : can_end+15] in the isoform
                           sequence (±5 residue slide).  VSP can_end is
                           used only as a starting hint for this search.

Both boundaries are domain-clipped and classified by structural element
(β-strand / α-helix / Loop / Inter-motif / Flanking) using the τ₅
classifier and compared against the length-weighted residue null π_t^0.

This analysis is the transcript-level counterpart of Analysis 4 (VSP
boundary enrichment, analyze_vsp_boundaries.py), which uses VSP can_start
and can_end annotations instead of sequence-derived boundaries.  Results
from the two analyses are directly comparable.

Proteins without a matched Ensembl transcript are excluded (see Results.md).

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
from junction_utils import load_canonical_junctions


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

MIN_MATCH = 15   # minimum AA window for suffix match to find resync
MAX_SLIDE = 5    # ±residues to slide can_end when searching for resync


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


def compute_residue_null(canonicals):
    counts = {c: 0 for c in CATS5}
    total  = 0
    for can in canonicals.values():
        ds, de = can["ds"], can["de"]
        if de <= ds:
            continue
        for idx in range(de - ds):
            counts[_tau5(ds + idx, can["motifs"])] += 1
            total += 1
    if total == 0:
        return {c: 1 / len(CATS5) for c in CATS5}
    return {c: counts[c] / total for c in CATS5}


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_canonicals(conn):
    """
    Load canonical proteins with sequences and motif annotations.
    Ensembl match is required only for the exclusion filter (13 proteins);
    sequences are needed for the divergence/resync search.
    """
    enst_matched = set(load_canonical_junctions(conn).keys())

    rows = conn.execute("""
        SELECT vc.uniprot_id, vc.domain_start, vc.domain_end,
               vc.motif_annotations, ca.sequence
        FROM   view_canonical vc
        JOIN   canonical_analysis ca ON ca.uniprot_id = vc.uniprot_id
        WHERE  vc.motif_annotations IS NOT NULL
          AND  vc.domain_start IS NOT NULL
          AND  vc.domain_end   IS NOT NULL
          AND  ca.sequence     IS NOT NULL
    """).fetchall()

    out = {}
    for uid, ds, de, ma, seq in rows:
        if uid not in enst_matched:
            continue
        out[uid] = dict(uid=uid, ds=ds, de=de, seq=seq, motifs=json.loads(ma))
    return out


def load_isoforms(conn):
    """
    Load non-canonical isoforms with protein sequences and VSP events.
    No Ensembl isoform-transcript match is required — the sequence
    comparison uses only the stored protein sequence.
    """
    rows = conn.execute("""
        SELECT i.isoform_id, i.sequence, nc.canonical_id, nc.vsp_domain_events
        FROM   isoforms i
        JOIN   view_noncanonical nc ON nc.isoform_id = i.isoform_id
        WHERE  i.is_canonical = 0
          AND  i.sequence IS NOT NULL
          AND  nc.vsp_domain_events IS NOT NULL
    """).fetchall()

    return [dict(iso_id=r[0], seq=r[1], can_id=r[2], vsps=json.loads(r[3]))
            for r in rows]


# ---------------------------------------------------------------------------
# Resync search
# ---------------------------------------------------------------------------

def find_resync(can_seq, can_end, iso_seq):
    """
    Search for canonical[can_end : can_end+MIN_MATCH] in iso_seq,
    sliding can_end by ±MAX_SLIDE.

    Returns (canonical_resync_pos, isoform_resync_pos, slide_offset) or None.
    canonical_resync_pos is the first shared position after the AS region.
    """
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
# Core: find transcript-derived AS boundaries
# ---------------------------------------------------------------------------

def find_as_boundaries(canonicals, isoforms):
    """
    For each isoform+VSP, compute the domain-clipped start and end of the
    alternatively spliced region:

      start = max(D_seq,   domain_start)
      end   = min(R_can-1, domain_end-1)   when resync is found (replacement)
            = min(can_end, domain_end-1)   when resync fails (truncating isoform)

    Truncating isoforms (where the isoform sequence does not rejoin the
    canonical) contribute can_end as the end boundary, identical to Analysis 4.
    Only records where the AS region overlaps the domain are included.
    Returns list of (canonical_uid, start_pos, end_pos).
    """
    records = []
    for iso in isoforms:
        can_id = iso["can_id"]
        if can_id not in canonicals:
            continue
        can     = canonicals[can_id]
        can_seq = can["seq"]
        iso_seq = iso["seq"]
        ds, de  = can["ds"], can["de"]

        # Divergence point D_seq (1-indexed canonical position)
        diverge = None
        for i, (ca, ia) in enumerate(zip(can_seq, iso_seq), start=1):
            if ca != ia:
                diverge = i
                break
        if diverge is None:
            continue

        for vsp in sorted(iso["vsps"], key=lambda v: v.get("can_start", 0)):
            can_end = vsp.get("can_end")
            if can_end is None:
                continue

            resync = find_resync(can_seq, can_end, iso_seq)
            if resync is not None:
                R_can, _R_iso, _slide = resync
                end = R_can - 1          # last diverged position (sequence-derived)
            else:
                end = can_end            # truncating isoform: use VSP annotation

            if diverge >= end:
                continue

            # Domain-clip both boundaries
            s_t = max(diverge, ds)
            e_t = min(end, de - 1)

            if s_t > e_t:               # AS region has no domain overlap
                continue

            records.append((can_id, s_t, e_t))
    return records


# ---------------------------------------------------------------------------
# Enrichment analysis
# ---------------------------------------------------------------------------

def compute_enrichment(positions, canonicals, pi_null):
    """
    positions : list of (canonical_uid, residue_position) tuples.
    Returns dict with N, N_t, pi_t, f_t, rho_t or None if N == 0.
    """
    N_t = {c: 0 for c in CATS5}
    N   = 0
    for uid, pos in positions:
        if uid not in canonicals:
            continue
        can = canonicals[uid]
        ds, de = can["ds"], can["de"]
        if not (ds <= pos < de):
            continue
        N_t[_tau5(pos, can["motifs"])] += 1
        N += 1
    if N == 0:
        return None
    f_t   = {c: N_t[c] / N for c in CATS5}
    rho_t = {c: f_t[c] / pi_null[c] if pi_null[c] > 0 else 0.0 for c in CATS5}
    return dict(N=N, N_t=N_t, pi_t=pi_null, f_t=f_t, rho_t=rho_t)


def chi_square_pvalues(result):
    pvals_raw = {}
    for c in CATS5:
        E_t = result["N"] * result["pi_t"][c]
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

def bh_ci_halfwidths(pvals_raw, pi_null, N, alpha=0.05):
    """CI half-widths consistent with BH-corrected stars: z_BH(t) / sqrt(N * pi_t)."""
    m     = len(CATS5)
    order = sorted(CATS5, key=lambda c: pvals_raw[c])
    ranks = {c: r for r, c in enumerate(order, 1)}
    errs  = {}
    for c in CATS5:
        alpha_eff = alpha * ranks[c] / m
        z_c = float(stats.norm.ppf(1 - alpha_eff / 2))
        E_t = N * pi_null[c]
        errs[c] = z_c / np.sqrt(E_t) if E_t > 0 else 0.0
    return errs


def plot_results(res_start, pv_raw_start, pv_bh_start, res_end, pv_raw_end, pv_bh_end, out):
    """
    Grouped bar chart: start and end enrichment side-by-side per structural element,
    matching the style of Analysis 4 (plot_combined_enrichment in analyze_vsp_boundaries.py).
    """
    w  = 0.35
    x  = np.arange(len(CATS5))
    xs = x - w / 2
    xe = x + w / 2

    rho_s = [res_start["rho_t"][c] for c in CATS5]
    rho_e = [res_end["rho_t"][c]   for c in CATS5]

    hw_s    = bh_ci_halfwidths(pv_raw_start, res_start["pi_t"], res_start["N"])
    hw_e    = bh_ci_halfwidths(pv_raw_end,   res_end["pi_t"],   res_end["N"])
    err_s_lo = [min(rho_s[i], hw_s[c]) for i, c in enumerate(CATS5)]
    err_s_hi = [hw_s[c] for c in CATS5]
    err_e_lo = [min(rho_e[i], hw_e[c]) for i, c in enumerate(CATS5)]
    err_e_hi = [hw_e[c] for c in CATS5]

    fig, ax = plt.subplots(figsize=(9, 4.5))

    ax.bar(xs, rho_s, w, color=[COLS5[c] for c in CATS5],
           alpha=0.90, zorder=3, label="AS start ($D_{\\rm seq}$)")
    ax.bar(xe, rho_e, w, color=[COLS5[c] for c in CATS5],
           alpha=0.45, zorder=3, hatch="//", label="AS end ($R_{\\rm can}-1$)")

    ax.errorbar(xs, rho_s, yerr=[err_s_lo, err_s_hi],
                fmt="none", color="black", capsize=4, lw=1.0, zorder=4)
    ax.errorbar(xe, rho_e, yerr=[err_e_lo, err_e_hi],
                fmt="none", color="black", capsize=4, lw=1.0, zorder=4)

    ax.axhline(1.0, color="black", lw=0.8, ls="--", zorder=2)

    for i, c in enumerate(CATS5):
        sl = sig_stars(pv_bh_start[c])
        el = sig_stars(pv_bh_end[c])
        if sl != "ns":
            ax.text(xs[i], rho_s[i] + err_s_hi[i] + 0.04, sl,
                    ha="center", va="bottom", fontsize=9, fontweight="bold")
        if el != "ns":
            ax.text(xe[i], rho_e[i] + err_e_hi[i] + 0.04, el,
                    ha="center", va="bottom", fontsize=9, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([LABELS5[c] for c in CATS5], fontsize=9)
    ax.set_ylabel(r"$\rho_t$ (observed / length-weighted null)", fontsize=10)
    ax.set_title(
        f"Transcript-derived AS boundary enrichment in TIM-barrel structural elements\n"
        f"$N$ = {res_start['N']} boundaries  |  "
        "error bars = BH-adjusted CI;  *, **, *** = BH-adjusted $p$",
        fontsize=9,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(fontsize=9, frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, -0.14), ncol=2)

    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")


def plot_pooled_enrichment(N_t_s, N_t_e, N_total, pi_null, pvals_bh, pvals_raw, out):
    """Single bar chart treating start and end positions as one pooled set."""
    N_t   = {c: N_t_s[c] + N_t_e[c] for c in CATS5}
    f_t   = {c: N_t[c] / N_total for c in CATS5}
    rho_t = {c: f_t[c] / pi_null[c] if pi_null[c] > 0 else 0.0 for c in CATS5}

    x        = np.arange(len(CATS5))
    rho_vals = [rho_t[c] for c in CATS5]
    hw       = bh_ci_halfwidths(pvals_raw, pi_null, N_total)
    err_lo   = [min(rho_vals[i], hw[c]) for i, c in enumerate(CATS5)]
    err_hi   = [hw[c] for c in CATS5]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(x, rho_vals, width=0.6, color=[COLS5[c] for c in CATS5],
           alpha=0.85, zorder=3)
    ax.errorbar(x, rho_vals, yerr=[err_lo, err_hi],
                fmt="none", color="black", capsize=5, lw=1.2, zorder=4)
    ax.axhline(1.0, color="black", lw=0.8, ls="--", zorder=2)

    for i, c in enumerate(CATS5):
        sl = sig_stars(pvals_bh[c])
        if sl != "ns":
            ax.text(x[i], rho_vals[i] + err_hi[i] + 0.04, sl,
                    ha="center", va="bottom", fontsize=10, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([LABELS5[c] for c in CATS5], fontsize=9)
    ax.set_ylabel(r"$\rho_t$ (observed / length-weighted null)", fontsize=10)
    ax.set_title(
        "Transcript-derived AS boundary enrichment in TIM-barrel structural elements\n"
        f"(start + end pooled, $N$ = {N_total}  |  error bars = BH-adjusted CI)",
        fontsize=9,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")
    return N_t, f_t, rho_t


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db",         default=None)
    parser.add_argument("--out",        default="figures/as_splice_junctions.png")
    parser.add_argument("--out-pooled", default="figures/as_splice_junctions_pooled.png")
    args = parser.parse_args()

    db_path = args.db or get_config().db_path
    conn    = sqlite3.connect(db_path)

    print("Loading data ...")
    canonicals = load_canonicals(conn)
    isoforms   = load_isoforms(conn)
    conn.close()

    print(f"  Canonical proteins (Ensembl-matched): {len(canonicals)}")
    print(f"  Non-canonical isoforms: {len(isoforms)}")

    pi_null = compute_residue_null(canonicals)
    print("\nResidue null (domain-length-weighted):")
    for c in CATS5:
        print(f"  {ALABELS[c]:<16}  {pi_null[c]:.4f}")

    print("\nFinding transcript-derived AS boundaries ...")
    boundaries = find_as_boundaries(canonicals, isoforms)

    starts = [(uid, s) for uid, s, e in boundaries]
    ends   = [(uid, e) for uid, s, e in boundaries]
    print(f"\n  AS boundary pairs found: {len(boundaries)}")
    print(f"  Start positions (domain-clipped): {len(starts)}")
    print(f"  End positions   (domain-clipped): {len(ends)}")

    def run(positions, label):
        res = compute_enrichment(positions, canonicals, pi_null)
        if res is None:
            print(f"\n  {label}: no domain-internal positions.")
            return None
        pv_raw, pv_bh = chi_square_pvalues(res)
        obs  = np.array([res["N_t"][c] for c in CATS5], dtype=float)
        exp  = np.array([res["N"] * res["pi_t"][c] for c in CATS5], dtype=float)
        chi2, p_gl = stats.chisquare(obs, f_exp=exp)
        print_enrichment(label, res, pv_raw, pv_bh, chi2, p_gl, len(CATS5) - 1)
        return res, pv_raw, pv_bh

    r_start = run(starts, "Transcript start positions (D_seq)")
    r_end   = run(ends,   "Transcript end positions (R_can - 1)")

    if r_start and r_end:
        res_s, pv_raw_s, pv_bh_s = r_start
        res_e, pv_raw_e, pv_bh_e = r_end
        plot_results(res_s, pv_raw_s, pv_bh_s, res_e, pv_raw_e, pv_bh_e, args.out)

        # Pooled: start + end as one set
        N_pooled = 2 * len(boundaries)
        pv_pool_raw = {}
        for c in CATS5:
            N_t_pool = res_s["N_t"][c] + res_e["N_t"][c]
            E_t = N_pooled * pi_null[c]
            z   = (N_t_pool - E_t) / np.sqrt(E_t) if pi_null[c] > 0 else 0.0
            pv_pool_raw[c] = float(2 * stats.norm.sf(abs(z)))
        pv_pool_bh = bh_correction(pv_pool_raw)
        N_t_pool_d, _, rho_pool = plot_pooled_enrichment(
            res_s["N_t"], res_e["N_t"], N_pooled, pi_null,
            pv_pool_bh, pv_pool_raw, args.out_pooled,
        )
        print("\nPooled boundary enrichment:")
        for c in CATS5:
            print(f"   {ALABELS[c]}: N={N_t_pool_d[c]}, rho={rho_pool[c]:.3f}, "
                  f"BH p={pv_pool_bh[c]:.4f}")


if __name__ == "__main__":
    main()
