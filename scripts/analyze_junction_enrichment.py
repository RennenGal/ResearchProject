#!/usr/bin/env python3
"""
Statistical analysis of exon junction placement in TIM-barrel structural elements.

Implements the formal framework from Statistical-Analysis.md:
  - Eligible junction positions E_p = {d_p^s, ..., d_p^e - 1}
  - Junction-count-weighted null expectation pi_t^0
  - Chi-square goodness-of-fit test (analytical approximation)
  - Per-element Pearson z-score p-values with Benjamini-Hochberg FDR correction
  - Bar-chart figure

Output figures:
  figures/enrichment_bars.png        — enrichment ratios (m=3 and m=5)

Usage:
    python scripts/analyze_junction_enrichment.py
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

def load_proteins(conn):
    rows = conn.execute("""
        SELECT uniprot_id, gene_name, domain_start, domain_end,
               exon_annotations, motif_annotations
        FROM   view_canonical
    """).fetchall()

    proteins = []
    for uid, gene, ds, de, ea, ma in rows:
        motifs = json.loads(ma)
        exons     = json.loads(ea)
        junctions = [e["end"] for e in exons[:-1] if ds <= e["end"] < de]
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
    N_t     = {c: 0   for c in cats}   # observed junction count per element type
    wq      = {c: 0.0 for c in cats}   # Σ_p n_p * q_pt  (weighted null accumulator → π_t^0)
    N       = 0                         # total domain-internal junctions across all proteins

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
# Per-element chi-square p-values (Pearson z-score)
# ---------------------------------------------------------------------------

def chi_square_pvalues(result):
    cats = result["cats"]
    N    = result["N"]
    pvals_raw = {}
    for c in cats:
        O_t = result["N_t"][c]
        E_t = N * result["pi_t"][c]
        if E_t > 0:
            z_t = (O_t - E_t) / np.sqrt(E_t)
            p   = 2 * stats.norm.sf(abs(z_t))
        else:
            p = 1.0
        pvals_raw[c] = float(p)
    pvals_bh = bh_correction(pvals_raw)
    return pvals_raw, pvals_bh


def bh_correction(pvals_dict):
    """Benjamini-Hochberg FDR correction."""
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


def print_results(label, n_prot, result, chi2, chi2_p, dof, pvals, pvals_bh):
    cats   = result["cats"]
    labels = ALABELS3 if len(cats) == 3 else ALABELS5
    N      = result["N"]
    print(f"\n{'='*78}")
    print(f"  {label}")
    print(f"  Proteins: {n_prot}  |  Total junctions (N): {N}  "
          f"|  Mean n_p: {N/n_prot:.1f}")
    print(f"{'='*78}")
    print(f"  {'Category':<16}  {'N_t':>5}  {'f_t':>7}  {'pi_t0':>7}  "
          f"{'rho_t':>6}  {'p_raw':>7}  {'p_BH':>7}  Sig")
    print("  " + "-"*16 + ("  " + "-"*5) + ("  " + "-"*7)*3
          + "  " + "-"*7 + "  " + "-"*7 + "  ---")
    for i, c in enumerate(cats):
        print(f"  {labels[c]:<16}  {result['N_t'][c]:>5}  "
              f"{result['f_t'][c]:>7.4f}  {result['pi_t'][c]:>7.4f}  "
              f"{result['rho_t'][c]:>6.3f}  "
              f"{pvals[c]:>7.4f}  {pvals_bh[c]:>7.4f}  {sig_stars(pvals_bh[c])}")
    print(f"\n  Chi-square (approx): chi2({dof}) = {chi2:.2f}, p = {chi2_p:.4g}")
    print(f"{'='*78}")


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def plot_enrichment_bars(result, pvals_bh, out="figures/enrichment_bars.png"):
    cats = result["cats"]
    N    = result["N"]
    x    = np.arange(len(cats))
    rhos = [result["rho_t"][c] for c in cats]
    errs = [1.96 * np.sqrt(result["rho_t"][c] / (N * result["pi_t"][c]))
            if result["pi_t"][c] > 0 else 0.0
            for c in cats]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(x, rhos, color=[COLS5[c] for c in cats],
           alpha=0.85, width=0.6, zorder=3)
    ax.axhline(1.0, color="black", linewidth=0.9, linestyle="--", zorder=2)

    # 95% Poisson CI error bars centred at rho_t
    ax.errorbar(x, rhos, yerr=errs,
                fmt="none", color="#888888", capsize=5,
                linewidth=1.2, zorder=4, label="95% CI (Poisson)")

    for i, c in enumerate(cats):
        star = sig_stars(pvals_bh[c])
        if star != "ns":
            ax.text(i, rhos[i] + errs[i] + 0.03, star,
                    ha="center", va="bottom", fontsize=10, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([LABELS5[c] for c in cats], rotation=20, ha="right", fontsize=9)
    ax.set_ylabel("Enrichment ratio $\\rho_t$", fontsize=10)
    ax.set_title(
        "Exon junction enrichment in TIM-barrel structural elements\n"
        "(error bars = 95% Poisson CI; *, **, *** = BH-adjusted p < 0.05/0.01/0.001)",
        fontsize=9,
    )
    ax.set_ylim(0, max(max(r + e for r, e in zip(rhos, errs)) + 0.25, 1.85))
    ax.yaxis.grid(True, linestyle=":", alpha=0.4, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# Update Statistical-Analysis.md §5
# ---------------------------------------------------------------------------

def _md_table(result, pvals_raw, pvals_bh):
    cats   = result["cats"]
    N      = result["N"]
    labels = LABELS3 if len(cats) == 3 else LABELS5
    rows   = ["| Category | $N_t$ | $f_t$ | $\\pi_t^0$ | $\\rho_t$ | $\\chi^2$ | $p$ (raw) | $p$ (BH) |",
              "|---|---|---|---|---|---|---|---|"]
    for i, c in enumerate(cats):
        O_t  = result["N_t"][c]
        E_t  = N * result["pi_t"][c]
        chi2 = (O_t - E_t) ** 2 / E_t if E_t > 0 else 0.0
        sig  = sig_stars(pvals_bh[c])
        rows.append(
            f"| {labels[c]} | {result['N_t'][c]} "
            f"| {result['f_t'][c]:.4f} | {result['pi_t'][c]:.4f} "
            f"| **{result['rho_t'][c]:.3f}** {sig} "
            f"| {chi2:.2f} "
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
    parser.add_argument("--out-bars",  default="figures/enrichment_bars.png")
    parser.add_argument("--md",        default="Statistical-Analysis.md")
    args = parser.parse_args()

    db_path  = args.db or get_config().db_path
    conn     = sqlite3.connect(db_path)
    proteins = load_proteins(conn)
    conn.close()
    print(f"Loaded {len(proteins)} proteins.")

    # -- Full model (m = 5) ------------------------------------------------
    print("\nFull model (m = 5) ...")
    r5               = run_analysis(proteins, model=5)
    chi2_5, p5, dof5 = chi_square_test(r5)
    pv5, pv5_bh      = chi_square_pvalues(r5)
    print_results("Full model (m = 5)", len(proteins),
                  r5, chi2_5, p5, dof5, pv5, pv5_bh)

    # -- Figures -----------------------------------------------------------
    plot_enrichment_bars(r5, pv5_bh, out=args.out_bars)

    # -- Update Statistical-Analysis.md §5 ---------------------------------
    N      = r5["N"]
    n_prot = len(proteins)
    content = f"""\
### Dataset

| | |
|---|---|
| Proteins analysed | {n_prot} |
| Total domain junctions ($N$) | {N} |
| Mean junctions per protein ($\\bar{{n}}_p$) | {N/n_prot:.1f} |
| Proteins with full annotation ($K_p = 8$) | {sum(1 for p in proteins if p['n_motifs'] == 8)} |
| Proteins with partial annotation ($K_p < 8$) | {sum(1 for p in proteins if p['n_motifs'] < 8)} |

### Full model ($m = 5$): β-strand, α-helix, inter-motif linker, loop (β→α), flanking

Chi-square global test (approximation): χ²({dof5}) = {chi2_5:.2f}, p = {p5:.4g}

{_md_table(r5, pv5, pv5_bh)}

Significance codes (BH-adjusted chi-square p-values):
\\* p < 0.05 \\*\\* p < 0.01 \\*\\*\\* p < 0.001

Figures: `figures/enrichment_bars.png` (enrichment ratios with 95% Poisson CI).
"""
    md_path = Path(args.md)
    if md_path.exists():
        update_results_section(str(md_path), content)
    else:
        print(f"  [note] {md_path} not found — markdown not updated.")


if __name__ == "__main__":
    main()
