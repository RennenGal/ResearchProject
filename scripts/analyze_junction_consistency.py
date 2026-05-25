#!/usr/bin/env python3
"""
Cross-protein junction consistency analysis (§7 of Statistical-Analysis.md).

  A. Global positional clustering   — KDE of normalised domain positions vs. permutation null
  B. Within-element phase           — KDE of within-element phases for α, β, loop

Output figures:
  figures/consistency_global.png
  figures/consistency_phase.png

Usage:
    python scripts/analyze_junction_consistency.py
    python scripts/analyze_junction_consistency.py --B 5000
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
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent.parent))
from protein_data_collector.config import get_config
from junction_utils import load_canonical_junctions


# ---------------------------------------------------------------------------
# Load proteins  (same filter as analyze_junction_enrichment.py)
# ---------------------------------------------------------------------------

def load_proteins(conn):
    enst_jcts = load_canonical_junctions(conn)
    rows = conn.execute("""
        SELECT uniprot_id, gene_name, domain_start, domain_end,
               motif_annotations
        FROM   view_canonical
        WHERE  motif_annotations IS NOT NULL
          AND  domain_start      IS NOT NULL
          AND  domain_end        IS NOT NULL
    """).fetchall()
    proteins = []
    for uid, gene, ds, de, ma in rows:
        if uid not in enst_jcts:
            continue                    # no Ensembl transcript — excluded
        junctions = enst_jcts[uid]
        if not junctions:
            continue
        motifs = json.loads(ma)
        proteins.append(dict(uid=uid, gene=gene, ds=ds, de=de,
                             motifs=motifs, n_motifs=len(motifs),
                             junctions=junctions, n_p=len(junctions)))
    return proteins


# ---------------------------------------------------------------------------
# Mean element spans for background shading
# ---------------------------------------------------------------------------

def mean_element_spans(proteins):
    spans = defaultdict(list)
    for p in proteins:
        ds, de = p["ds"], p["de"]
        dlen   = de - ds
        if dlen <= 0:
            continue
        for i, m in enumerate(p["motifs"]):
            n   = m["motif"]
            r   = lambda x: (x - ds) / dlen
            spans[f"beta_{n}"].append( (r(m["beta_start"]),      r(m["beta_end"])) )
            spans[f"loop_{n}"].append( (r(m["beta_end"] + 1),    r(m["alpha_start"] - 1)) )
            spans[f"alpha_{n}"].append((r(m["alpha_start"]),      r(m["alpha_end"])) )
            if i + 1 < len(p["motifs"]):
                nxt = p["motifs"][i + 1]
                spans[f"inter_{n}_{n+1}"].append(
                    (r(m["alpha_end"] + 1), r(nxt["beta_start"] - 1)))
    return {k: (np.mean([s for s, e in v]), np.mean([e for s, e in v]))
            for k, v in spans.items() if v}


# ---------------------------------------------------------------------------
# A. Global positional clustering
# ---------------------------------------------------------------------------

def normalized_positions(proteins):
    xs = []
    for p in proteins:
        ds, de = p["ds"], p["de"]
        dlen   = de - ds
        for j in p["junctions"]:
            xs.append((j - ds) / dlen)
    return np.array(xs)


def permutation_positions(proteins, B, seed):
    rng = np.random.default_rng(seed)
    N   = sum(p["n_p"] for p in proteins)
    all_xs = np.zeros((B, N))
    for b in range(B):
        idx = 0
        for p in proteins:
            ds, de = p["ds"], p["de"]
            dlen   = de - ds
            drawn  = rng.choice(dlen, size=p["n_p"], replace=False)
            all_xs[b, idx : idx + p["n_p"]] = drawn / dlen
            idx += p["n_p"]
    return all_xs


def plot_global_clustering(x_obs, perm_xs, ks_obs, ks_perm_p, spans, n_prot, out):
    grid     = np.linspace(0, 1, 500)
    bw       = "scott"
    kde_obs  = stats.gaussian_kde(x_obs, bw_method=bw)(grid)
    kde_null = np.array([stats.gaussian_kde(perm_xs[b], bw_method=bw)(grid)
                         for b in range(len(perm_xs))])
    null_lo  = np.percentile(kde_null,  2.5, axis=0)
    null_hi  = np.percentile(kde_null, 97.5, axis=0)
    null_mean= kde_null.mean(axis=0)

    BG = {"beta":  ("#4C72B0", 0.18), "alpha": ("#DD8452", 0.18),
          "loop":  ("#eeeeee", 0.80), "inter": ("#eeeeee", 0.80)}

    fig, ax = plt.subplots(figsize=(10, 4.5))
    for elem, (s, e) in spans.items():
        if e <= s:
            continue
        prefix = elem.split("_")[0]
        col, alpha = BG.get(prefix, ("#eeeeee", 0.5))
        ax.axvspan(s, e, color=col, alpha=alpha, zorder=0)

    ax.fill_between(grid, null_lo, null_hi, color="#888888", alpha=0.25,
                    zorder=2, label="Null 95% interval")
    ax.plot(grid, null_mean, color="#888888", lw=0.9, ls="--", zorder=3)
    ax.plot(grid, kde_obs,   color="black",   lw=1.8, zorder=4, label="Observed KDE")
    ax.axhline(1.0, color="#aaaaaa", lw=0.5, ls=":", zorder=1)

    import matplotlib.patches as mpatches
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
        f"Global junction density — {n_prot} proteins, N = {len(x_obs)} junctions\n"
        f"KS: D = {ks_obs:.4f},  permutation p = {ks_perm_p:.4f}",
        fontsize=10,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# B. Within-element phase
# ---------------------------------------------------------------------------

def _tau5_with_interval(pos, motifs):
    """
    Assign pos to exactly one element using the same first-match priority as
    _tau5() in analyze_junction_enrichment.py, and return (type, s, e) for
    the matched element instance.
    """
    if not motifs or pos < motifs[0]["beta_start"]:
        return "flanking", pos, pos
    for i, m in enumerate(motifs):
        if m["beta_start"] <= pos <= m["beta_end"]:
            return "beta", m["beta_start"], m["beta_end"]
        if m["beta_end"] < pos < m["alpha_start"]:
            return "loop", m["beta_end"] + 1, m["alpha_start"] - 1
        if m["alpha_start"] <= pos <= m["alpha_end"]:
            return "alpha", m["alpha_start"], m["alpha_end"]
        if i + 1 < len(motifs) and m["alpha_end"] < pos < motifs[i + 1]["beta_start"]:
            return "inter", m["alpha_end"] + 1, motifs[i + 1]["beta_start"] - 1
    return "flanking", pos, pos


def collect_phases(proteins, element_type):
    """
    For each junction whose §5-consistent exclusive assignment is element_type,
    return its phase phi = (j - s) / (e - s).  Element instances of length 1
    (e == s) are excluded because phi is undefined.
    """
    phases = []
    for p in proteins:
        for j in p["junctions"]:
            t, s, e = _tau5_with_interval(j, p["motifs"])
            if t == element_type and e > s:
                phases.append((j - s) / (e - s))
    return np.array(phases)


def permutation_phases(proteins, element_type, B, seed):
    """
    For each junction exclusively assigned to element_type, re-draw its
    position uniformly within the same element instance [s, e].
    Preserves which element instances received junctions; tests within-element
    uniformity only, not element-level enrichment.
    """
    rng  = np.random.default_rng(seed)
    lens = []
    for p in proteins:
        for j in p["junctions"]:
            t, s, e = _tau5_with_interval(j, p["motifs"])
            if t == element_type and e > s:
                lens.append(e - s)
    if not lens:
        return np.empty((B, 0))
    lens = np.array(lens)
    perm = np.zeros((B, len(lens)))
    for b in range(B):
        offsets = rng.integers(0, lens + 1, size=len(lens))
        perm[b] = offsets / lens
    return perm


def plot_phase_consistency(proteins, perm_phases_dict, ks_results, out):
    elements = [
        ("alpha", "α-helix",    "#DD8452"),
        ("beta",  "β-strand",   "#4C72B0"),
        ("loop",  "Loop (β→α)", "#8c8c8c"),
    ]
    grid = np.linspace(0, 1, 300)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))

    for ax, (etype, label, col) in zip(axes, elements):
        phases          = collect_phases(proteins, etype)
        n               = len(phases)
        ks_D, ks_p, perm_p = ks_results[etype]

        if n < 5:
            ax.set_title(f"{label}\n(n = {n}, insufficient data)", fontsize=10)
            ax.set_visible(True)
            continue

        kde = stats.gaussian_kde(phases, bw_method="scott")(grid)

        pp = perm_phases_dict[etype]
        if pp.shape[1] > 0:
            null_kdes = np.array([
                stats.gaussian_kde(pp[b], bw_method="scott")(grid)
                for b in range(len(pp))
            ])
            ax.fill_between(grid,
                            np.percentile(null_kdes,  2.5, axis=0),
                            np.percentile(null_kdes, 97.5, axis=0),
                            color="#888888", alpha=0.20, zorder=2,
                            label="Null 95% interval")

        ax.plot(grid, kde, color=col, lw=1.8, zorder=3, label="Observed")
        ax.axhline(1.0, color="black", lw=0.8, ls="--", zorder=4, label="Uniform")
        ax.text(0.97, 0.97,
                f"n = {n}\nKS D = {ks_D:.3f}\np = {ks_p:.4g}\nperm p = {perm_p:.4f}",
                transform=ax.transAxes, ha="right", va="top", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85))
        ax.set_title(label, fontsize=10, fontweight="bold")
        ax.set_xlabel(r"Phase within element $\phi$", fontsize=9)
        ax.set_ylabel("Density", fontsize=9)
        ax.set_xlim(0, 1)
        ax.set_xticks([0, 0.5, 1.0])
        ax.set_xticklabels(["0\n(start)", "0.5", "1\n(end)"], fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if etype == "alpha":
            ax.legend(fontsize=8, frameon=False, loc="upper center",
                      bbox_to_anchor=(0.5, -0.45), ncol=3)

    fig.suptitle(
        r"Within-element phase of exon junctions  ($\phi = 0$: start, $\phi = 1$: end)"
        "\nGrey band = permutation null 95% interval",
        fontsize=10, y=1.02,
    )
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# Update §7 results in Statistical-Analysis.md
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


def update_section8(md_path, content_8a, content_8b):
    text = Path(md_path).read_text(encoding="utf-8")
    text = _replace_block(text,
                          "#### 8A Results\n",
                          "\n---\n\n### 8B.",
                          content_8a)
    text = _replace_block(text,
                          "#### 8B Results\n",
                          "\n---\n\n### 8C.",
                          content_8b)
    Path(md_path).write_text(text, encoding="utf-8")
    print(f"Updated §8 in {md_path}")


def _row(etype, label, ks_results, proteins):
    D, p, pp = ks_results[etype]
    n = len(collect_phases(proteins, etype))
    return f"| {label} | {n} | {D:.4f} | {p:.4g} | {pp:.4f} |"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db",         default=None)
    parser.add_argument("--B",          type=int, default=2000)
    parser.add_argument("--seed",       type=int, default=42)
    parser.add_argument("--out-global", default="figures/consistency_global.png")
    parser.add_argument("--out-phase",  default="figures/consistency_phase.png")
    parser.add_argument("--md",         default="Statistical-Analysis.md")
    args = parser.parse_args()

    db_path  = args.db or get_config().db_path
    conn     = sqlite3.connect(db_path)
    proteins = load_proteins(conn)
    conn.close()
    N = sum(p["n_p"] for p in proteins)
    print(f"Loaded {len(proteins)} proteins, {N} junctions total.")

    # ── A. Global positional clustering ──────────────────────────────────────
    print(f"\nA. Global positional clustering (B = {args.B}) …")
    x_obs  = normalized_positions(proteins)
    ks_res = stats.kstest(x_obs, "uniform")
    print(f"   KS: D = {ks_res.statistic:.4f}, p = {ks_res.pvalue:.4g}")

    print("   Computing permutation null …")
    perm_xs   = permutation_positions(proteins, args.B, args.seed)
    ks_null   = np.array([stats.kstest(perm_xs[b], "uniform").statistic
                          for b in range(args.B)])
    ks_perm_p = (1 + np.sum(ks_null >= ks_res.statistic)) / (args.B + 1)
    print(f"   Permutation p = {ks_perm_p:.4f}")

    spans = mean_element_spans(proteins)
    plot_global_clustering(x_obs, perm_xs, ks_res.statistic, ks_perm_p,
                           spans, len(proteins), args.out_global)

    # ── B. Within-element phase ───────────────────────────────────────────────
    print(f"\nB. Within-element phase (B = {args.B}) …")
    ks_results       = {}
    perm_phases_dict = {}
    for etype in ("alpha", "beta", "loop"):
        ph  = collect_phases(proteins, etype)
        n   = len(ph)
        ks  = stats.kstest(ph, "uniform") if n >= 5 else None
        pp  = permutation_phases(proteins, etype, args.B, args.seed)
        if ks is not None and pp.shape[1] > 0:
            perm_ks = np.array([stats.kstest(pp[b], "uniform").statistic
                                for b in range(args.B)])
            perm_p  = (1 + np.sum(perm_ks >= ks.statistic)) / (args.B + 1)
        else:
            perm_p = float("nan")
        ks_results[etype]       = (
            ks.statistic if ks else float("nan"),
            ks.pvalue    if ks else float("nan"),
            perm_p,
        )
        perm_phases_dict[etype] = pp
        D, p = (ks.statistic, ks.pvalue) if ks else (float("nan"), float("nan"))
        print(f"   {etype}: n = {n}, KS D = {D:.4f}, p = {p:.4g}, perm p = {perm_p:.4f}")

    plot_phase_consistency(proteins, perm_phases_dict, ks_results, args.out_phase)

    # BH correction across the three permutation p-values
    etypes_b = ["alpha", "beta", "loop"]
    perm_pvals = np.array([ks_results[e][2] for e in etypes_b])
    order      = np.argsort(perm_pvals)
    adj        = np.zeros(3)
    for rank, idx in enumerate(order, 1):
        adj[idx] = perm_pvals[idx] * 3 / rank
    for j in range(1, -1, -1):
        adj[order[j]] = min(adj[order[j]], adj[order[j + 1]])
    adj = np.minimum(adj, 1.0)
    bh_pvals = {e: float(adj[i]) for i, e in enumerate(etypes_b)}
    for e in etypes_b:
        print(f"   {e}: perm p = {ks_results[e][2]:.4f}, BH p = {bh_pvals[e]:.4f}")

    # ── Update §7 ─────────────────────────────────────────────────────────────
    content_7a = f"""\
| Statistic | Value |
|---|---|
| Proteins | {len(proteins)} |
| Total junctions ($N$) | {N} |
| KS statistic $D_N$ | {ks_res.statistic:.4f} |
| KS $p$ (analytical, descriptive only) | {ks_res.pvalue:.4g} |
| Permutation $p$ ($B = {args.B}$) | {ks_perm_p:.4f} |

![Global junction density across TIM-barrel domain](figures/consistency_global.png)"""

    def _row_bh(etype, label):
        D, p, pp = ks_results[etype]
        n  = len(collect_phases(proteins, etype))
        bh = bh_pvals[etype]
        return f"| {label} | {n} | {D:.4f} | {p:.4g} | {pp:.4f} | {bh:.4f} |"

    content_7b = f"""\
| Element | $n$ | KS $D$ | KS $p$ (analytical) | Perm $p$ | Perm $p$ (BH) |
|---|---|---|---|---|---|
{_row_bh("alpha", "α-helix")}
{_row_bh("beta",  "β-strand")}
{_row_bh("loop",  "Loop (β→α)")}

Significance codes apply to BH-adjusted permutation p-values.

![Within-element phase distributions](figures/consistency_phase.png)"""

    md_path = Path(args.md)
    if md_path.exists():
        update_section8(str(md_path), content_7a, content_7b)
    else:
        print(f"  [note] {md_path} not found — markdown not updated.")


if __name__ == "__main__":
    main()
