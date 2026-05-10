#!/usr/bin/env python3
"""
Statistical analysis of exon junction placement in TIM-barrel structural elements.

Implements the formal framework from Statistical-Analysis.md:
  - Eligible junction positions E_p = {d_p^s, ..., d_p^e - 1}
  - Junction-count-weighted null expectation pi_t^0
  - Chi-square goodness-of-fit test (analytical approximation)
  - Within-protein permutation test (preferred, B replicates)
  - Benjamini–Hochberg FDR correction on per-element p-values
  - Bar-chart and permutation-distribution figures

Output figures:
  figures/enrichment_bars.png        — enrichment ratios (m=3 and m=5)
  figures/enrichment_permutation.png — permutation null distributions (m=5)

Usage:
    python scripts/analyze_junction_enrichment.py
    python scripts/analyze_junction_enrichment.py --full-only   # K_p=8 sensitivity check
    python scripts/analyze_junction_enrichment.py --B 10000
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
# Category definitions
# ---------------------------------------------------------------------------

CATS5   = ["beta", "alpha", "inter", "loop", "flanking"]
CATS3   = ["beta", "alpha", "other"]
LABELS5 = {"beta": "β-strand", "alpha": "α-helix", "inter": "Inter-motif",
           "loop": "Loop (β→α)", "flanking": "Flanking"}
LABELS3 = {"beta": "β-strand", "alpha": "α-helix", "other": "Other"}
# ASCII equivalents for console printing on Windows
ALABELS5 = {"beta": "beta-strand", "alpha": "alpha-helix", "inter": "Inter-motif",
            "loop": "Loop (b->a)", "flanking": "Flanking"}
ALABELS3 = {"beta": "beta-strand", "alpha": "alpha-helix", "other": "Other"}
COLS5   = {"beta": "#4C72B0", "alpha": "#DD8452", "inter": "#55A868",
           "loop": "#8c8c8c", "flanking": "#C44E52"}
COLS3   = {"beta": "#4C72B0", "alpha": "#DD8452", "other": "#8c8c8c"}


# ---------------------------------------------------------------------------
# Per-protein type array
# ---------------------------------------------------------------------------

def _tau5(pos, motifs):
    """Return the 5-category label for position pos (first-match, lower motif wins)."""
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


def _tau3(pos, motifs):
    t = _tau5(pos, motifs)
    return t if t in ("beta", "alpha") else "other"


def build_type_array(ds, de, motifs, model, cat_idx):
    """
    Return an int8 array of length |E_p| = de - ds.
    Entry i gives the category index of position ds+i under model m.
    Uses per-position classification so overlapping motif annotations
    are handled correctly (first-match / lower-numbered motif wins).
    """
    ep_len = de - ds
    arr    = np.zeros(ep_len, dtype=np.int8)
    fn     = _tau3 if model == 3 else _tau5
    for idx in range(ep_len):
        arr[idx] = cat_idx[fn(ds + idx, motifs)]
    return arr


# ---------------------------------------------------------------------------
# Load proteins
# ---------------------------------------------------------------------------

def load_proteins(conn, full_only=False):
    rows = conn.execute("""
        SELECT uniprot_id, gene_name, domain_start, domain_end,
               exon_annotations, motif_annotations
        FROM   canonical_analysis
        WHERE  exon_annotations  IS NOT NULL
          AND  motif_annotations IS NOT NULL
          AND  domain_start      IS NOT NULL
          AND  domain_end        IS NOT NULL
    """).fetchall()

    proteins = []
    for uid, gene, ds, de, ea, ma in rows:
        motifs = json.loads(ma)
        if full_only and len(motifs) != 8:
            continue
        exons     = json.loads(ea)
        junctions = [e["end"] for e in exons[:-1] if ds <= e["end"] < de]
        if not junctions:
            continue
        proteins.append({"uid": uid, "gene": gene, "ds": ds, "de": de,
                         "motifs": motifs, "n_motifs": len(motifs),
                         "junctions": junctions, "n_p": len(junctions)})
    return proteins


# ---------------------------------------------------------------------------
# Core analysis — enrichment ratios
# ---------------------------------------------------------------------------

def run_analysis(proteins, model):
    cats    = CATS3 if model == 3 else CATS5
    cat_idx = {c: i for i, c in enumerate(cats)}
    N_t     = {c: 0   for c in cats}
    wq      = {c: 0.0 for c in cats}   # Σ_p n_p * q_pt
    N       = 0

    for p in proteins:
        ds, de = p["ds"], p["de"]
        ep_len = de - ds                       # |E_p|
        arr    = build_type_array(ds, de, p["motifs"], model, cat_idx)

        # Lambda_t^p — eligible positions per element
        for i, c in enumerate(cats):
            wq[c] += p["n_p"] * float(np.sum(arr == i)) / ep_len

        # Observed junction counts (look up in pre-built array for consistency)
        for j in p["junctions"]:
            t = cats[arr[j - ds]]
            N_t[t] += 1
            N += 1

    pi_t  = {c: wq[c] / N                        for c in cats}
    f_t   = {c: N_t[c] / N                        for c in cats}
    rho_t = {c: f_t[c] / pi_t[c] if pi_t[c] > 0 else 0.0 for c in cats}
    return {"N": N, "N_t": N_t, "pi_t": pi_t, "f_t": f_t,
            "rho_t": rho_t, "cats": cats, "cat_idx": cat_idx}


# ---------------------------------------------------------------------------
# Chi-square test (analytical approximation)
# ---------------------------------------------------------------------------

def chi_square_test(result):
    cats     = result["cats"]
    N        = result["N"]
    observed = np.array([result["N_t"][c]       for c in cats], dtype=float)
    expected = np.array([N * result["pi_t"][c]  for c in cats], dtype=float)
    chi2, p  = stats.chisquare(observed, f_exp=expected)
    return chi2, p, len(cats) - 1


# ---------------------------------------------------------------------------
# Within-protein permutation test
# ---------------------------------------------------------------------------

def permutation_test(proteins, model, result, B=10000, seed=42):
    cats    = result["cats"]
    cat_idx = result["cat_idx"]
    n_cats  = len(cats)
    N       = result["N"]
    pi_t    = result["pi_t"]
    rng     = np.random.default_rng(seed)

    # Pre-build type arrays once (reused across all B replicates)
    type_arrays = [
        (p["n_p"], build_type_array(p["ds"], p["de"], p["motifs"], model, cat_idx))
        for p in proteins
    ]

    rho_b = np.zeros((B, n_cats))
    for b in range(B):
        counts_b = np.zeros(n_cats, dtype=np.int64)
        for n_p, arr in type_arrays:
            drawn     = rng.choice(len(arr), size=n_p, replace=False)
            counts_b += np.bincount(arr[drawn].astype(int), minlength=n_cats)
        f_b = counts_b / N
        for i, c in enumerate(cats):
            rho_b[b, i] = f_b[i] / pi_t[c] if pi_t[c] > 0 else 0.0

    return rho_b


def compute_pvalues(result, rho_b):
    """Two-sided p-values using |rho - 1| distance, with +1 finite-sample correction."""
    cats  = result["cats"]
    B     = len(rho_b)
    return {
        c: (1 + int(np.sum(np.abs(rho_b[:, i] - 1) >= abs(result["rho_t"][c] - 1)))) / (B + 1)
        for i, c in enumerate(cats)
    }


def bh_correction(pvals_dict):
    """Benjamini–Hochberg FDR correction."""
    cats  = list(pvals_dict.keys())
    pvals = np.array([pvals_dict[c] for c in cats])
    n     = len(pvals)
    order = np.argsort(pvals)
    adj   = np.zeros(n)
    for rank, idx in enumerate(order, 1):
        adj[idx] = pvals[idx] * n / rank
    for j in range(n - 2, -1, -1):
        adj[order[j]] = min(adj[order[j]], adj[order[j + 1]])
    adj = np.minimum(adj, 1.0)
    return {c: float(adj[i]) for i, c in enumerate(cats)}


# ---------------------------------------------------------------------------
# Console output
# ---------------------------------------------------------------------------

def sig_stars(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"


def print_results(label, n_prot, result, chi2, chi2_p, dof, pvals, pvals_bh, rho_b):
    cats   = result["cats"]
    labels = ALABELS3 if len(cats) == 3 else ALABELS5
    N      = result["N"]
    print(f"\n{'='*78}")
    print(f"  {label}")
    print(f"  Proteins: {n_prot}  |  Total junctions (N): {N}  "
          f"|  Mean n_p: {N/n_prot:.1f}")
    print(f"{'='*78}")
    print(f"  {'Category':<16}  {'N_t':>5}  {'f_t':>7}  {'pi_t0':>7}  "
          f"{'rho_t':>6}  {'95% CI permutation':>20}  {'p_raw':>7}  {'p_BH':>7}  Sig")
    print("  " + "-"*16 + ("  " + "-"*5) + ("  " + "-"*7)*3
          + "  " + "-"*20 + "  " + "-"*7 + "  " + "-"*7 + "  ---")
    for i, c in enumerate(cats):
        lo = np.percentile(rho_b[:, i], 2.5)
        hi = np.percentile(rho_b[:, i], 97.5)
        print(f"  {labels[c]:<16}  {result['N_t'][c]:>5}  "
              f"{result['f_t'][c]:>7.4f}  {result['pi_t'][c]:>7.4f}  "
              f"{result['rho_t'][c]:>6.3f}  "
              f"  [{lo:.3f}, {hi:.3f}]        "
              f"{pvals[c]:>7.4f}  {pvals_bh[c]:>7.4f}  {sig_stars(pvals_bh[c])}")
    print(f"\n  Chi-square (approx): chi2({dof}) = {chi2:.2f}, p = {chi2_p:.4g}")
    print(f"{'='*78}")


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def plot_enrichment_bars(r3, r5, rb3, rb5, p3_bh, p5_bh,
                         out="figures/enrichment_bars.png"):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    for ax, result, rho_b, pvals_bh, cols, labels, title in [
        (axes[0], r3, rb3, p3_bh, COLS3, LABELS3, "Simplified model (m = 3)"),
        (axes[1], r5, rb5, p5_bh, COLS5, LABELS5, "Full model (m = 5)"),
    ]:
        cats    = result["cats"]
        x       = np.arange(len(cats))
        rhos    = [result["rho_t"][c] for c in cats]
        null_lo = [np.percentile(rho_b[:, i], 2.5)  for i in range(len(cats))]
        null_hi = [np.percentile(rho_b[:, i], 97.5) for i in range(len(cats))]

        ax.bar(x, rhos, color=[cols[c] for c in cats],
               alpha=0.85, width=0.6, zorder=3)
        ax.axhline(1.0, color="black", linewidth=0.9, linestyle="--", zorder=2)

        # Null 95% CI as grey error bars centred at rho=1
        ax.errorbar(x, [1.0] * len(cats),
                    yerr=[[1.0 - lo for lo in null_lo],
                          [hi - 1.0 for hi in null_hi]],
                    fmt="none", color="#888888", capsize=5,
                    linewidth=1.2, zorder=4, label="Null 95% CI")

        for i, c in enumerate(cats):
            star = sig_stars(pvals_bh[c])
            if star != "ns":
                ax.text(i, max(rhos[i], null_hi[i]) + 0.03, star,
                        ha="center", va="bottom", fontsize=10, fontweight="bold")

        ax.set_xticks(x)
        ax.set_xticklabels([labels[c] for c in cats],
                           rotation=20, ha="right", fontsize=9)
        ax.set_ylabel("Enrichment ratio $\\rho_t$", fontsize=10)
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_ylim(0, max(max(rhos) + 0.25, max(null_hi) + 0.15, 1.85))
        ax.yaxis.grid(True, linestyle=":", alpha=0.4, zorder=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.suptitle(
        "Exon junction enrichment in TIM-barrel structural elements\n"
        "(error bars = 95% permutation null interval; *, **, *** = BH-adjusted p < 0.05/0.01/0.001)",
        fontsize=10, y=1.02,
    )
    fig.tight_layout()
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")


def plot_permutation_distributions(result, rho_b, pvals_raw, pvals_bh,
                                   out="figures/enrichment_permutation.png"):
    cats   = result["cats"]
    labels = LABELS3 if len(cats) == 3 else LABELS5
    cols   = COLS3   if len(cats) == 3 else COLS5
    ncols  = 3
    nrows  = (len(cats) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = np.array(axes).flatten()

    for i, c in enumerate(cats):
        ax  = axes[i]
        obs = result["rho_t"][c]
        ax.hist(rho_b[:, i], bins=60, color=cols[c], alpha=0.65,
                edgecolor="white", linewidth=0.3, density=True, zorder=2)
        ax.axvline(obs, color="black", linewidth=2.0, zorder=4)
        ax.axvline(1.0, color="gray",  linewidth=0.8, linestyle="--", zorder=3)
        tx = 0.03 if obs > 1.0 else 0.97
        ha = "left" if obs > 1.0 else "right"
        ax.text(tx, 0.95,
                f"$\\rho$ = {obs:.3f}\n"
                f"$p$ = {pvals_raw[c]:.4f}\n"
                f"$p_{{BH}}$ = {pvals_bh[c]:.4f}   {sig_stars(pvals_bh[c])}",
                transform=ax.transAxes, ha=ha, va="top", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
        ax.set_title(labels[c], fontsize=10, fontweight="bold")
        ax.set_xlabel("$\\rho_t$ (permuted)", fontsize=8)
        ax.set_ylabel("Density", fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(
        f"Permutation null distributions  (B = {len(rho_b):,} replicates)\n"
        "Vertical line = observed $\\rho_t$",
        fontsize=11, y=1.01,
    )
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# Update Statistical-Analysis.md §5
# ---------------------------------------------------------------------------

def _md_table(result, rho_b, pvals_raw, pvals_bh):
    cats   = result["cats"]
    labels = LABELS3 if len(cats) == 3 else LABELS5
    rows   = ["| Category | $N_t$ | $f_t$ | $\\pi_t^0$ | $\\rho_t$ | Null 95% interval | $p$ (raw) | $p$ (BH) |",
              "|---|---|---|---|---|---|---|---|"]
    for i, c in enumerate(cats):
        lo   = np.percentile(rho_b[:, i], 2.5)
        hi   = np.percentile(rho_b[:, i], 97.5)
        sig  = sig_stars(pvals_bh[c])
        rows.append(
            f"| {labels[c]} | {result['N_t'][c]} "
            f"| {result['f_t'][c]:.4f} | {result['pi_t'][c]:.4f} "
            f"| **{result['rho_t'][c]:.3f}** {sig} "
            f"| [{lo:.3f}, {hi:.3f}] "
            f"| {pvals_raw[c]:.4f} | {pvals_bh[c]:.4f} |"
        )
    return "\n".join(rows)


def update_results_section(md_path, content):
    text      = Path(md_path).read_text(encoding="utf-8")
    start_tag = "## 5. Results\n"
    end_tag   = "\n## 6."
    si = text.find(start_tag)
    ei = text.find(end_tag, si)
    if si == -1 or ei == -1:
        print(f"  [warn] Could not locate §5 Results in {md_path} — skipping update.")
        return
    new_text = text[: si + len(start_tag)] + "\n" + content + "\n" + text[ei:]
    Path(md_path).write_text(new_text, encoding="utf-8")
    print(f"Updated §5 Results in {md_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db",        default=None)
    parser.add_argument("--B",         type=int, default=10000,
                        help="Permutation replicates (default 10 000)")
    parser.add_argument("--seed",      type=int, default=42)
    parser.add_argument("--full-only", action="store_true",
                        help="Restrict to proteins with K_p = 8 (sensitivity check)")
    parser.add_argument("--out-bars",  default="figures/enrichment_bars.png")
    parser.add_argument("--out-perm",  default="figures/enrichment_permutation.png")
    parser.add_argument("--md",        default="Statistical-Analysis.md")
    args = parser.parse_args()

    db_path  = args.db or get_config().db_path
    conn     = sqlite3.connect(db_path)
    proteins = load_proteins(conn, full_only=args.full_only)
    conn.close()
    print(f"Loaded {len(proteins)} proteins "
          f"({'K_p = 8 only' if args.full_only else 'full + partial'}).")

    # ── Simplified model (m = 3) ──────────────────────────────────────────
    print("\nSimplified model (m = 3) …")
    r3          = run_analysis(proteins, model=3)
    chi2_3, p3, dof3 = chi_square_test(r3)
    print(f"  Permutation test (B = {args.B}) …")
    rb3         = permutation_test(proteins, 3, r3, B=args.B, seed=args.seed)
    pv3         = compute_pvalues(r3, rb3)
    pv3_bh      = bh_correction(pv3)
    print_results("Simplified model (m = 3)", len(proteins),
                  r3, chi2_3, p3, dof3, pv3, pv3_bh, rb3)

    # ── Full model (m = 5) ────────────────────────────────────────────────
    print("\nFull model (m = 5) …")
    r5          = run_analysis(proteins, model=5)
    chi2_5, p5, dof5 = chi_square_test(r5)
    print(f"  Permutation test (B = {args.B}) …")
    rb5         = permutation_test(proteins, 5, r5, B=args.B, seed=args.seed)
    pv5         = compute_pvalues(r5, rb5)
    pv5_bh      = bh_correction(pv5)
    print_results("Full model (m = 5)", len(proteins),
                  r5, chi2_5, p5, dof5, pv5, pv5_bh, rb5)

    # ── Figures ───────────────────────────────────────────────────────────
    plot_enrichment_bars(r3, r5, rb3, rb5, pv3_bh, pv5_bh, out=args.out_bars)
    plot_permutation_distributions(r5, rb5, pv5, pv5_bh, out=args.out_perm)

    # ── Update Statistical-Analysis.md §5 ─────────────────────────────────
    N       = r3["N"]
    n_prot  = len(proteins)
    content = f"""\
### Dataset

| | |
|---|---|
| Proteins analysed | {n_prot} |
| Total domain junctions ($N$) | {N} |
| Mean junctions per protein ($\\bar{{n}}_p$) | {N/n_prot:.1f} |
| Proteins with full annotation ($K_p = 8$) | {sum(1 for p in proteins if p['n_motifs'] == 8)} |
| Proteins with partial annotation ($K_p < 8$) | {sum(1 for p in proteins if p['n_motifs'] < 8)} |

### Simplified model ($m = 3$): β-strand, α-helix, other

Chi-square global test (approximation): χ²({dof3}) = {chi2_3:.2f}, p = {p3:.4g}

{_md_table(r3, rb3, pv3, pv3_bh)}

### Full model ($m = 5$): β-strand, α-helix, inter-motif linker, loop (β→α), flanking

Chi-square global test (approximation): χ²({dof5}) = {chi2_5:.2f}, p = {p5:.4g}

{_md_table(r5, rb5, pv5, pv5_bh)}

Significance codes (BH-adjusted permutation p-values):
\\* p < 0.05 \\*\\* p < 0.01 \\*\\*\\* p < 0.001

Figures: `figures/enrichment_bars.png` (enrichment ratios with 95% permutation CI),
`figures/enrichment_permutation.png` (permutation null distributions, full model).
"""
    md_path = Path(args.md)
    if md_path.exists():
        update_results_section(str(md_path), content)
    else:
        print(f"  [note] {md_path} not found — markdown not updated.")


if __name__ == "__main__":
    main()
