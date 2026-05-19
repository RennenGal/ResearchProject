#!/usr/bin/env python3
"""
Motif-specific element enrichment analysis (§8A replacement).

For each exon junction, assigns a (element_type, motif_number) label
τ_motif(j, p) ∈ {(β,1),(loop,1),(α,1),(β,2),...,(α,8)} instead of the
five-category τ5 used in §6.  This tests whether junctions recur at the
same named structural position (e.g. loop 3, α-helix 5) across proteins
of different lengths — a question that normalised-domain position (old §8A)
cannot address because it conflates structural position with domain length.

The null follows §5: for each category the expected count under the
length-weighted null is E_{(t,k)} = N * π_{(t,k)}^0.  A z-score
z_{(t,k)} = (N_{(t,k)} - E_{(t,k)}) / sqrt(E_{(t,k)}) is converted to a
two-sided p-value via the standard normal, and BH correction is applied
across all 31 primary categories.

Primary categories: (β,k), (loop,k), (α,k) for k = 1..8  plus
(inter,k) for k = 1..7 — 31 categories total.
Flanking is tracked but excluded (positions before/after the annotated
repeat region are not motif-numbered).
BH correction is applied across all 31 primary categories.

Output:
  figures/motif_enrichment_heatmap.png  — 3×8 heatmap of ρ values
  Statistical-Analysis.md               — §8A section updated in-place

Usage:
    python scripts/analyze_motif_enrichment.py
"""

import argparse
import json
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

sys.path.insert(0, str(Path(__file__).parent.parent))
from protein_data_collector.config import get_config


# ---------------------------------------------------------------------------
# Category encoding
# ---------------------------------------------------------------------------

ELEM_TYPES  = ["beta", "loop", "alpha"]
N_MOTIFS    = 8
INTER_BASE  = 24   # inter_{k,k+1} → INTER_BASE + (k-1), k=1..7
FLANKING    = 31
N_TOTAL_CATS = 32

# Primary categories: 24 (β/loop/α × 1..8) + 7 inter-motif (1..7) = 31
_BETA_LOOP_ALPHA = [(t, k) for k in range(1, N_MOTIFS + 1) for t in ELEM_TYPES]
_INTER_CATS      = [("inter", k) for k in range(1, N_MOTIFS)]   # inter_{k,k+1}
PRIMARY_CATS     = _BETA_LOOP_ALPHA + _INTER_CATS

# Integer index maps
CAT_IDX = {(t, k): 3 * (k - 1) + i
           for k in range(1, N_MOTIFS + 1) for i, t in enumerate(ELEM_TYPES)}
CAT_IDX.update({("inter", k): INTER_BASE + (k - 1) for k in range(1, N_MOTIFS)})
IDX_CAT = {v: k for k, v in CAT_IDX.items()}

ELEM_LABELS  = {"beta": "β-strand", "loop": "Loop (β→α)",
                "alpha": "α-helix",  "inter": "Inter-motif"}
ELEM_ALABELS = {"beta": "beta-strand", "loop": "Loop (b->a)",
                "alpha": "alpha-helix", "inter": "Inter-motif"}
ELEM_COLORS  = {"beta": "#4C72B0", "loop": "#8c8c8c",
                "alpha": "#DD8452", "inter": "#55A868"}


# ---------------------------------------------------------------------------
# τ_motif: assign (element_type, motif_number) to a domain position
# ---------------------------------------------------------------------------

def _tau_motif_idx(pos, motifs):
    """
    Return integer category index for pos given the protein's motif list.
    Categories 0..23: primary (beta/loop/alpha × motifs 1..8)
    Categories 24..30: inter-motif linkers (between motif k and k+1)
    Category 31: flanking (before first motif or after last)
    """
    if not motifs or pos < motifs[0]["beta_start"]:
        return FLANKING
    for i, m in enumerate(motifs):
        k = m["motif"]
        if m["beta_start"] <= pos <= m["beta_end"]:
            return CAT_IDX.get(("beta", k), FLANKING)
        if m["beta_end"] < pos < m["alpha_start"]:
            return CAT_IDX.get(("loop", k), FLANKING)
        if m["alpha_start"] <= pos <= m["alpha_end"]:
            return CAT_IDX.get(("alpha", k), FLANKING)
        if i + 1 < len(motifs):
            nxt = motifs[i + 1]
            if m["alpha_end"] < pos < nxt["beta_start"]:
                inter_idx = INTER_BASE + (k - 1)
                return inter_idx if inter_idx < FLANKING else FLANKING
    return FLANKING


def build_type_array(ds, de, motifs):
    """
    Int8 array of length (de - ds) where entry i = τ_motif_idx(ds + i, motifs).
    Built once per protein; reused across all B permutation replicates.
    """
    ep_len = de - ds
    arr = np.zeros(ep_len, dtype=np.int8)
    for idx in range(ep_len):
        arr[idx] = _tau_motif_idx(ds + idx, motifs)
    return arr


# ---------------------------------------------------------------------------
# Load proteins (identical filter to §6)
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
        proteins.append(dict(uid=uid, gene=gene, ds=ds, de=de,
                             motifs=motifs, n_motifs=len(motifs),
                             junctions=junctions, n_p=len(junctions)))
    return proteins


# ---------------------------------------------------------------------------
# Core enrichment computation
# ---------------------------------------------------------------------------

def run_analysis(proteins):
    """
    Returns observed counts, null expectations, enrichment ratios, and
    per-category protein counts for the 24 primary categories.
    """
    N_t   = defaultdict(int)     # observed junction counts
    wq    = defaultdict(float)   # Σ_p n_p * q_p(t,k)
    n_prot_cat = defaultdict(int)  # proteins contributing to each category
    N     = 0

    for p in proteins:
        ds, de = p["ds"], p["de"]
        ep_len = de - ds
        if ep_len <= 0:
            continue
        arr = build_type_array(ds, de, p["motifs"])

        # Null expectation: length-weighted contribution of each category
        for cat_i in range(N_TOTAL_CATS):
            cnt = int(np.sum(arr == cat_i))
            if cnt > 0:
                key = IDX_CAT.get(cat_i)
                if key is not None:  # primary category
                    wq[key]         += p["n_p"] * cnt / ep_len
                    n_prot_cat[key] += 1

        # Observed junction counts
        for j in p["junctions"]:
            cat_i = int(arr[j - ds])
            key   = IDX_CAT.get(cat_i)
            if key is not None:
                N_t[key] += 1
            N += 1

    pi_t  = {cat: wq[cat] / N              for cat in PRIMARY_CATS}
    f_t   = {cat: N_t[cat] / N             for cat in PRIMARY_CATS}
    rho_t = {cat: (f_t[cat] / pi_t[cat])
             if pi_t[cat] > 0 else float("nan")
             for cat in PRIMARY_CATS}

    return dict(N=N, N_t=N_t, pi_t=pi_t, f_t=f_t,
                rho_t=rho_t, n_prot_cat=n_prot_cat)


# ---------------------------------------------------------------------------
# Chi-square p-values (replaces permutation test)
# ---------------------------------------------------------------------------

def chi_square_pvalues(result):
    """
    For each of the 31 primary categories compute a z-score against the
    length-weighted null, convert to a two-sided p-value via the standard
    normal, and BH-correct across all categories.

    Returns (pvals_raw, pvals_bh, z_scores) — all dicts keyed by category.
    """
    N    = result["N"]
    N_t  = result["N_t"]
    pi_t = result["pi_t"]

    z_scores  = {}
    pvals_raw = {}
    for cat in PRIMARY_CATS:
        pi = pi_t[cat]
        if pi <= 0:
            z_scores[cat]  = float("nan")
            pvals_raw[cat] = float("nan")
            continue
        E_c = N * pi
        z   = (N_t[cat] - E_c) / np.sqrt(E_c)
        z_scores[cat]  = float(z)
        pvals_raw[cat] = float(2 * norm.sf(abs(z)))

    pvals_bh = bh_correction(pvals_raw)
    return pvals_raw, pvals_bh, z_scores


def bh_correction(pvals):
    cats  = [c for c in PRIMARY_CATS if not np.isnan(pvals.get(c, float("nan")))]
    arr   = np.array([pvals[c] for c in cats])
    n     = len(arr)
    order = np.argsort(arr)
    adj   = np.zeros(n)
    for rank, idx in enumerate(order, 1):
        adj[idx] = arr[idx] * n / rank
    for j in range(n - 2, -1, -1):
        adj[order[j]] = min(adj[order[j]], adj[order[j + 1]])
    adj = np.minimum(adj, 1.0)
    result = {c: float("nan") for c in PRIMARY_CATS}
    for i, c in enumerate(cats):
        result[c] = float(adj[i])
    return result


def sig_stars(p):
    if np.isnan(p):   return ""
    if p < 0.001:     return "***"
    if p < 0.01:      return "**"
    if p < 0.05:      return "*"
    return ""


# ---------------------------------------------------------------------------
# Heatmap figure
# ---------------------------------------------------------------------------

def plot_heatmap(result, pvals_bh, n_prot_full, out):
    """
    4 × 8 heatmap: rows = β/loop/α/inter, columns = motifs 1–8.
    Inter-motif has only 7 positions (k=1..7); column 8 is left blank.
    Color encodes ρ (diverging, centred at 1).  Cell text shows ρ and
    significance stars.
    """
    all_rows  = ELEM_TYPES + ["inter"]
    rho_grid  = np.full((4, 8), np.nan)
    sig_grid  = [[""] * 8 for _ in range(4)]
    nprot_row = np.zeros((4, 8), dtype=int)

    for ri, elem in enumerate(all_rows):
        n_k = 7 if elem == "inter" else 8
        for k in range(1, n_k + 1):
            cat = (elem, k)
            rho_grid[ri, k - 1]  = result["rho_t"].get(cat, float("nan"))
            sig_grid[ri][k - 1]  = sig_stars(pvals_bh.get(cat, float("nan")))
            nprot_row[ri, k - 1] = result["n_prot_cat"].get(cat, 0)

    # Diverging colormap centred at 1.0
    vmax = max(0.3, float(np.nanmax(np.abs(rho_grid - 1)))) + 0.05
    vmin, vmax_plot = 1 - vmax, 1 + vmax

    fig, ax = plt.subplots(figsize=(11, 5.0))
    im = ax.imshow(rho_grid, cmap="RdBu_r",
                   vmin=vmin, vmax=vmax_plot, aspect="auto")

    # Cell annotations
    for ri in range(4):
        for ci in range(8):
            rho = rho_grid[ri, ci]
            star = sig_grid[ri][ci]
            if not np.isnan(rho):
                color = "white" if abs(rho - 1) > vmax * 0.6 else "black"
                label = f"{rho:.2f}"
                if star:
                    label += f"\n{star}"
                ax.text(ci, ri, label, ha="center", va="center",
                        fontsize=8, color=color, fontweight="bold" if star else "normal")

    ax.set_xticks(range(8))
    ax.set_xticklabels(
        [f"Motif {k}\n(n={nprot_row[0, k-1]})" for k in range(1, 9)],
        fontsize=8,
    )
    ax.set_yticks(range(4))
    ax.set_yticklabels([ELEM_LABELS[e] for e in ELEM_TYPES + ["inter"]], fontsize=9)
    ax.set_title(
        "Motif-specific junction enrichment  ($\\rho_{(t,k)}$ = observed / null fraction)\n"
        "* BH $p$ < 0.05   ** $p$ < 0.01   *** $p$ < 0.001",
        fontsize=9,
    )

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Enrichment ratio $\\rho$", fontsize=8)
    cbar.ax.axhline(1.0, color="black", lw=0.8, ls="--")

    fig.tight_layout()
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# Update §8A in Statistical-Analysis.md
# ---------------------------------------------------------------------------

def build_md_content(result, pvals_raw, pvals_bh, z_scores, n_prot):
    N = result["N"]

    # Table rows grouped by element type
    def cat_label(t, k):
        return f"{ELEM_LABELS[t]} {k}"

    header = (
        "| Position | $N_{(t,k)}$ | $f_{(t,k)}$ | $\\pi_{(t,k)}^0$ | "
        "$\\rho_{(t,k)}$ | $\\chi^2$ | Raw $p$ | BH $p$ | Sig |\n"
        "|---|---|---|---|---|---|---|---|---|"
    )

    rows = []
    for elem in ELEM_TYPES + ["inter"]:
        for k in range(1, 9 if elem != "inter" else 8):
            cat = (elem, k)
            rho = result["rho_t"][cat]
            if np.isnan(rho):
                continue
            z   = z_scores.get(cat, float("nan"))
            chi2 = z ** 2 if not np.isnan(z) else float("nan")
            p_r = pvals_raw[cat]
            p_b = pvals_bh[cat]
            n_p = result["n_prot_cat"].get(cat, 0)
            rows.append(
                f"| {cat_label(elem, k)} (n={n_p}) "
                f"| {result['N_t'][cat]} "
                f"| {result['f_t'][cat]:.4f} "
                f"| {result['pi_t'][cat]:.4f} "
                f"| {rho:.3f} "
                f"| {chi2:.3f} "
                f"| {p_r:.4f} "
                f"| {p_b:.4f} "
                f"| {sig_stars(p_b) or 'ns'} |"
            )
        rows.append("|---|---|---|---|---|---|---|---|---|")

    sig_cats = [(cat, result["rho_t"][cat], pvals_bh[cat])
                for cat in PRIMARY_CATS
                if not np.isnan(pvals_bh.get(cat, float("nan")))
                and pvals_bh[cat] < 0.05]
    sig_cats.sort(key=lambda x: x[2])

    if sig_cats:
        sig_lines = "  ".join(
            f"{ELEM_LABELS[t]} {k} ($\\rho={rho:.3f}$, BH $p={p:.4f}$)"
            for (t, k), rho, p in sig_cats
        )
        sig_summary = f"Significant positions (BH $p < 0.05$): {sig_lines}."
    else:
        sig_summary = "No individual motif-element position reaches significance after BH correction."

    return f"""\
Dataset: {n_prot} proteins, {N} domain-internal junctions.

{header}
{chr(10).join(rows)}

Significance codes (BH-adjusted chi-square p-values):
\\* p < 0.05  \\*\\* p < 0.01  \\*\\*\\* p < 0.001

{sig_summary}

![Motif-specific element enrichment heatmap](figures/motif_enrichment_heatmap.png)"""


def update_8a(md_path, new_content):
    text     = Path(md_path).read_text(encoding="utf-8")
    start    = "### 8A."
    end      = "\n---\n\n### 8B."
    si = text.find(start)
    ei = text.find(end, si)
    if si == -1 or ei == -1:
        print(f"  [warn] Could not locate §8A in {md_path}")
        return

    new_8a = (
        "### 8A. Motif-specific element enrichment\n\n"
        "#### Setup\n\n"
        "For each junction $j \\in \\mathcal{J}_p$, assign the **motif-element label**\n\n"
        "$$\\tau_\\text{motif}(j, p) = (t,\\, k)$$\n\n"
        "where $t \\in \\{\\beta, \\text{loop}, \\alpha\\}$ is the structural element type "
        "and $k \\in \\{1, \\ldots, K_p\\}$ is the motif number.  This extends the "
        "five-category $\\tau_5$ of §6 by preserving which repeat unit the junction falls in.  "
        "A junction at loop 3 in a short protein and loop 3 in a long protein receive the same "
        "label $(\\text{loop}, 3)$ regardless of their normalised domain coordinates.\n\n"
        "The 31 **primary categories** are $(t, k)$ for $t \\in \\{\\beta, \\text{loop}, \\alpha\\}$ "
        "and $k = 1, \\ldots, 8$ (24 categories), plus $(\\text{inter}, k)$ for $k = 1, \\ldots, 7$ "
        "(7 categories representing the linker between motif $k$ and motif $k+1$).  "
        "Flanking positions (before the first motif or after the last) are tracked but "
        "excluded from the primary test because they are not assigned to a specific repeat unit.  "
        "For proteins with $K_p < 8$, categories for motifs $k > K_p$ receive no contribution "
        "from that protein.\n\n"
        "The null expectation follows §5:\n\n"
        "$$\\pi_{(t,k)}^0 = \\frac{1}{N} \\sum_{p} n_p\\, q_{p,(t,k)}$$\n\n"
        "where $q_{p,(t,k)} = |E_{p,(t,k)}| / L_p$ is the fraction of domain positions "
        "in element $(t, k)$ for protein $p$, and the sum runs only over proteins that "
        "have motif $k$.  The enrichment ratio and BH-corrected chi-square test are "
        "applied across all 31 primary categories.  For each category, "
        "$z_{(t,k)} = (N_{(t,k)} - E_{(t,k)}) / \\sqrt{E_{(t,k)}}$ where "
        "$E_{(t,k)} = N \\pi_{(t,k)}^0$, and $p = 2\\,\\Phi(-|z|)$.\n\n"
        "#### 8A Results\n\n"
        + new_content
    )

    text = text[:si] + new_8a + text[ei:]
    Path(md_path).write_text(text, encoding="utf-8")
    print(f"Updated §8A in {md_path}")


# ---------------------------------------------------------------------------
# Also update §8 intro paragraph to reflect the new 8A question
# ---------------------------------------------------------------------------

def update_8_intro(md_path):
    text  = Path(md_path).read_text(encoding="utf-8")
    old   = (
        "Section 8 asks whether exon junctions show positional consistency across TIM-barrel domains\n"
        "beyond the element-level enrichment tested in §6. Analysis 8A tests global clustering along\n"
        "normalised domain coordinates. Analysis 8B tests whether, conditional on falling inside a\n"
        "specific structural element type, junctions preferentially occur near the beginning, middle,\n"
        "or end of that element."
    )
    new   = (
        "Section 8 asks whether exon junctions show positional consistency across TIM-barrel domains\n"
        "beyond the element-level enrichment tested in §6. Analysis 8A tests whether junctions recur\n"
        "at the same named structural position — specifically the same motif-element combination\n"
        "$(t, k)$ (e.g.\\ loop 3, $\\alpha$-helix 5) — across proteins of different lengths,\n"
        "using a 24-category motif-specific enrichment test. Analysis 8B tests whether, conditional\n"
        "on falling inside a specific structural element type, junctions preferentially occur near\n"
        "the beginning, middle, or end of that element."
    )
    if old not in text:
        print("  [warn] Could not locate §8 intro paragraph — skipping update.")
        return
    Path(md_path).write_text(text.replace(old, new), encoding="utf-8")
    print(f"Updated §8 intro in {md_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db",        default=None)
    parser.add_argument("--out-heatmap", default="figures/motif_enrichment_heatmap.png")
    parser.add_argument("--md",        default="Statistical-Analysis.md")
    args = parser.parse_args()

    db_path  = args.db or get_config().db_path
    conn     = sqlite3.connect(db_path)
    proteins = load_proteins(conn)
    conn.close()

    n_prot      = len(proteins)
    n_full      = sum(1 for p in proteins if p["n_motifs"] == 8)
    n_partial   = n_prot - n_full
    print(f"Loaded {n_prot} proteins ({n_full} full K_p=8, {n_partial} partial).")

    # Restrict to fully-annotated proteins: partial proteins concentrate junctions
    # into fewer motif slots, skewing enrichment estimates at specific (t,k) positions.
    proteins = [p for p in proteins if p["n_motifs"] == 8]
    print(f"Restricting to K_p=8 proteins: {len(proteins)} proteins for analysis.")

    print("\nComputing motif-specific enrichment ...")
    result = run_analysis(proteins)
    print(f"  Total junctions N = {result['N']}")

    print("\nComputing chi-square p-values ...")
    pvals_raw, pvals_bh, z_scores = chi_square_pvalues(result)

    # Console summary
    print(f"\n{'='*70}")
    print(f"  Motif-specific enrichment  (BH-adjusted chi-square)")
    print(f"{'='*70}")
    print(f"  {'Position':<18}  {'N_t':>5}  {'rho':>6}  {'z':>7}  {'p_raw':>7}  {'p_BH':>7}  Sig")
    for elem in ELEM_TYPES + ["inter"]:
        for k in range(1, 9 if elem != "inter" else 8):
            cat = (elem, k)
            rho = result["rho_t"][cat]
            if np.isnan(rho):
                continue
            label = f"{ELEM_ALABELS[elem][:5]} {k}"
            z = z_scores.get(cat, float("nan"))
            print(f"  {label:<18}  {result['N_t'][cat]:>5}  {rho:>6.3f}"
                  f"  {z:>7.3f}  {pvals_raw[cat]:>7.4f}  {pvals_bh[cat]:>7.4f}  "
                  f"{sig_stars(pvals_bh[cat]) or 'ns'}")
        print()

    plot_heatmap(result, pvals_bh, n_full, args.out_heatmap)

    md_path = Path(args.md)
    if md_path.exists():
        md_content = build_md_content(result, pvals_raw, pvals_bh, z_scores, n_prot)
        update_8a(str(md_path), md_content)
        update_8_intro(str(md_path))
    else:
        print(f"  [note] {md_path} not found — markdown not updated.")


if __name__ == "__main__":
    main()
