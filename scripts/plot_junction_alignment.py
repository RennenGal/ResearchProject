#!/usr/bin/env python3
"""
Domain junction alignment plot.

For proteins with domain lengths in a given range, draws one horizontal
row per protein (sorted by domain length).  Each exon junction inside the
domain is marked as a coloured dot whose x-position is the junction's
fractional position within the domain (0 = domain_start, 1 = domain_end).
Dot colour reflects the structural element at that junction.

A background stripe shows the mean span of each structural element
(beta / loop / alpha) across all proteins in the panel, computed as the
average normalised start and end of every element across proteins.

Usage
-----
    python scripts/plot_junction_alignment.py
    python scripts/plot_junction_alignment.py --min 200 --max 400 --step 50
    python scripts/plot_junction_alignment.py --min 200 --max 400 --no-subgroups
    python scripts/plot_junction_alignment.py --out figures/junctions.png
"""

import argparse
import json
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))
from protein_data_collector.config import get_config

# ---------------------------------------------------------------------------
# Colours
# ---------------------------------------------------------------------------
COLOURS = {
    "in_motif_beta":  "#4C72B0",   # blue
    "in_motif_loop":  "#8c8c8c",   # grey
    "in_motif_alpha": "#DD8452",   # orange
    "between_motifs": "#55A868",   # green
    "flanking":       "#C44E52",   # red
}
BG_COLOURS = {
    "beta":  "#5B9BD5",   # strong blue
    "loop":  "#f0f0f0",   # near-white
    "alpha": "#E8714A",   # strong orange-red
    "inter": "#f0f0f0",   # same near-white as loop
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def classify(pos, motifs):
    if not motifs:
        return "flanking"
    if pos < motifs[0]["beta_start"]:
        return "flanking"
    for i, m in enumerate(motifs):
        n = m["motif"]
        if m["beta_start"] <= pos <= m["beta_end"]:
            return "in_motif_beta"
        if m["beta_end"] < pos < m["alpha_start"]:
            return "in_motif_loop"
        if m["alpha_start"] <= pos <= m["alpha_end"]:
            return "in_motif_alpha"
        if i + 1 < len(motifs):
            if m["alpha_end"] < pos < motifs[i + 1]["beta_start"]:
                return "between_motifs"
    return "flanking"


def norm(pos, ds, de):
    """Normalise absolute position to [0, 1] within [ds, de]."""
    return (pos - ds) / (de - ds) if de > ds else 0.0


def mean_element_spans(rows):
    """
    For a set of proteins compute mean normalised [start, end] of every
    structural element (beta_N, loop_N, alpha_N, inter_N_{N+1}).
    Returns dict: element_label -> (mean_start, mean_end)
    """
    spans = defaultdict(list)
    for uid, gene, ds, de, ea, ma, dl in rows:
        if not ma:
            continue
        motifs = json.loads(ma)
        dlen = de - ds
        if dlen <= 0:
            continue
        for i, m in enumerate(motifs):
            n = m["motif"]
            spans[f"beta_{n}"].append(
                (norm(m["beta_start"], ds, de), norm(m["beta_end"], ds, de)))
            spans[f"loop_{n}"].append(
                (norm(m["beta_end"] + 1, ds, de), norm(m["alpha_start"] - 1, ds, de)))
            spans[f"alpha_{n}"].append(
                (norm(m["alpha_start"], ds, de), norm(m["alpha_end"], ds, de)))
            if i + 1 < len(motifs):
                nxt = motifs[i + 1]
                spans[f"inter_{n}_{n+1}"].append(
                    (norm(m["alpha_end"] + 1, ds, de),
                     norm(nxt["beta_start"] - 1, ds, de)))

    return {k: (np.mean([s for s, e in v]), np.mean([e for s, e in v]))
            for k, v in spans.items() if v}


def draw_panel(ax, rows, title, show_yticks=True):
    """Draw one subplot panel for a group of proteins."""
    # Sort by domain length so similar sizes are adjacent
    rows = sorted(rows, key=lambda r: r[6])

    # Background element spans
    spans = mean_element_spans(rows)
    bg_order = sorted(spans.keys(),
                      key=lambda k: spans[k][0])  # left to right
    for elem, (s, e) in spans.items():
        if e <= s:
            continue
        if elem.startswith("beta_"):
            colour, alpha = BG_COLOURS["beta"],  0.25
        elif elem.startswith("loop_"):
            colour, alpha = BG_COLOURS["loop"],  0.80
        elif elem.startswith("alpha_"):
            colour, alpha = BG_COLOURS["alpha"], 0.22
        else:
            colour, alpha = BG_COLOURS["inter"], 0.80
        ax.axvspan(s, e, color=colour, alpha=alpha, zorder=0)

    # Motif labels along the top — beta just above, alpha a step higher
    n_motifs = max(
        (int(k.split("_")[1]) for k in spans
         if k.startswith(("beta_", "alpha_"))),
        default=8,
    )
    y_offsets = {"beta": 1.01, "alpha": 1.018}
    for m in range(1, n_motifs + 1):
        for prefix, short in [("beta", "β"), ("alpha", "α")]:
            key = f"{prefix}_{m}"
            if key in spans:
                s, e = spans[key]
                mid = (s + e) / 2
                ax.text(mid, y_offsets[prefix], f"{short}{m}",
                        ha="center", va="bottom", fontsize=5.5,
                        color="#2c5f8a" if prefix == "beta" else "#8b3a1e",
                        fontweight="bold", zorder=5,
                        transform=ax.get_xaxis_transform())

    # One row per protein
    n_proteins = len(rows)
    for row_idx, (uid, gene, ds, de, ea, ma, dl) in enumerate(rows):
        if not ea or not ma:
            continue
        exons  = json.loads(ea)
        motifs = json.loads(ma)
        y      = row_idx

        for exon in exons[:-1]:
            pos = exon["end"]
            if ds <= pos < de:
                x     = norm(pos, ds, de)
                label = classify(pos, motifs)
                ax.scatter(x, y, color=COLOURS[label], s=30, zorder=3,
                           linewidths=0, alpha=0.9)

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-1, n_proteins)
    ax.set_title(title, fontsize=10, fontweight="bold", pad=28)
    ax.set_xlabel("Normalised domain position", fontsize=8)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"], fontsize=7)

    if show_yticks:
        ax.set_ylabel("Proteins (sorted by domain length)", fontsize=8)
        ax.set_yticks([])
    else:
        ax.set_yticks([])

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)



# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db",           default=None)
    parser.add_argument("--min",          type=int, default=200)
    parser.add_argument("--max",          type=int, default=400)
    parser.add_argument("--step",         type=int, default=50)
    parser.add_argument("--no-subgroups", action="store_true",
                        help="Single panel for the whole range")
    parser.add_argument("--out",          default="figures/junction_alignment.png")
    args = parser.parse_args()

    db_path = args.db or get_config().db_path
    conn    = sqlite3.connect(db_path)
    rows    = conn.execute("""
        SELECT uniprot_id, gene_name, domain_start, domain_end,
               exon_annotations, motif_annotations,
               domain_end - domain_start + 1 AS domain_len
        FROM canonical_analysis
        WHERE exon_annotations  IS NOT NULL
          AND motif_annotations IS NOT NULL
          AND domain_start      IS NOT NULL
          AND domain_end        IS NOT NULL
        ORDER BY domain_len
    """).fetchall()
    conn.close()

    in_range = [r for r in rows if args.min <= r[6] < args.max]
    print(f"Proteins {args.min}–{args.max} aa with exon+motif data: {len(in_range)}")

    # Build panels
    if args.no_subgroups:
        panels = [("", in_range)]
    else:
        panels = []
        for lo in range(args.min, args.max, args.step):
            hi     = lo + args.step
            bucket = [r for r in in_range if lo <= r[6] < hi]
            if bucket:
                panels.append((f"{lo}–{hi} aa", bucket))

    n_panels = len(panels)
    if n_panels == 0:
        print("No data found.")
        return

    # Dynamically set figure height based on max proteins in a single panel
    max_rows = max(len(p[1]) for p in panels)
    row_height = 0.10 if max_rows > 100 else 0.22
    fig_h = max(6, max_rows * row_height + 2)
    fig_w = max(6, 4.5 * n_panels)

    fig, axes = plt.subplots(
        1, n_panels,
        figsize=(fig_w, fig_h),
        sharey=False,
    )
    if n_panels == 1:
        axes = [axes]

    for i, (title, panel_rows) in enumerate(panels):
        draw_panel(axes[i], panel_rows, title, show_yticks=(i == 0))

    # Legend
    legend_patches = [
        mpatches.Patch(color=COLOURS["in_motif_beta"],  label="β-strand junction"),
        mpatches.Patch(color=COLOURS["in_motif_loop"],  label="Loop junction"),
        mpatches.Patch(color=COLOURS["in_motif_alpha"], label="α-helix junction"),
        mpatches.Patch(color=COLOURS["between_motifs"], label="Inter-motif junction"),
        mpatches.Patch(color=COLOURS["flanking"],       label="Flanking junction"),
    ]
    bg_patches = [
        mpatches.Patch(facecolor=BG_COLOURS["beta"],  alpha=0.5, label="β-strand region (mean)"),
        mpatches.Patch(facecolor=BG_COLOURS["alpha"], alpha=0.5, label="α-helix region (mean)"),
        mpatches.Patch(facecolor=BG_COLOURS["loop"],  alpha=0.8, label="Loop / inter-motif (mean)"),
    ]
    fig.legend(
        handles=legend_patches + bg_patches,
        loc="lower center",
        ncol=4,
        fontsize=8,
        frameon=False,
        bbox_to_anchor=(0.5, -0.04),
    )

    fig.suptitle(
        "Exon junction positions across TIM barrel domain\n"
        "(each row = one protein, sorted by domain length within panel)",
        fontsize=11, y=1.01,
    )

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"Saved to {args.out}")


if __name__ == "__main__":
    main()
