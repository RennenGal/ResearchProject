#!/usr/bin/env python3
"""
Standalone per-isoform disruption heatmap for HYAL1.

Usage:
    python scripts/plot_hyal1_heatmap.py
"""

import json
import sqlite3
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist

sys.path.insert(0, str(Path(__file__).parent.parent))
from protein_data_collector.config import get_config

GENE      = "HYAL1"
OUT       = Path("figures/alt") / f"{GENE.lower()}_isoform_heatmap.png"
N_POS     = 8


# ---------------------------------------------------------------------------
# Data helpers (mirrors run_alternative_analysis.py)
# ---------------------------------------------------------------------------

def load_canonical(conn):
    rows = conn.execute("""
        SELECT vc.uniprot_id, vc.gene_name, vc.domain_start, vc.domain_end,
               vc.motif_annotations, ca.sequence
        FROM   view_canonical vc
        JOIN   canonical_analysis ca ON ca.uniprot_id  = vc.uniprot_id
                                    AND ca.domain_index = vc.domain_index
        WHERE  vc.domain_index = 1
    """).fetchall()
    out = {}
    for uid, gene, ds, de, mj, seq in rows:
        out[uid] = dict(gene=gene or uid, ds=ds, de=de,
                        motifs=json.loads(mj), seq=seq)
    return out


def load_isoforms(conn):
    rows = conn.execute("""
        SELECT nc.isoform_id, nc.uniprot_id, nc.vsp_domain_events,
               i.sequence_length, i.sequence
        FROM   view_noncanonical nc
        JOIN   isoforms i ON i.isoform_id = nc.isoform_id
        WHERE  nc.vsp_domain_events != '[]'
          AND  i.sequence IS NOT NULL
    """).fetchall()
    return [dict(isoform_id=iso_id, uniprot_id=uid,
                 vsps=json.loads(vj), iso_len=l, seq=s)
            for iso_id, uid, vj, l, s in rows]


def is_core_affected(vsps, motifs):
    if not motifs:
        return False
    core_s = motifs[0]["beta_start"]
    core_e = motifs[-1]["alpha_end"]
    return any(v["can_end"] >= core_s and v["can_start"] <= core_e for v in vsps)


def classify_combined(iso, canonicals):
    uid = iso["uniprot_id"]
    states = []
    for m in canonicals[uid]["motifs"]:
        state = "intact"
        for vsp in iso["vsps"]:
            vs, ve = vsp["can_start"], vsp["can_end"]
            if ve < m["beta_start"] or vs > m["alpha_end"]:
                continue
            if vs <= m["beta_start"] and ve >= m["alpha_end"]:
                state = "removed"
                break
            state = "partial"
        states.append(state)
    return states


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    cfg  = get_config()
    conn = sqlite3.connect(cfg.db_path)

    canonicals = load_canonical(conn)
    all_isos   = load_isoforms(conn)

    # Strict filter: VSP must overlap motif core
    filtered = [iso for iso in all_isos
                if iso["uniprot_id"] in canonicals
                and is_core_affected(iso["vsps"], canonicals[iso["uniprot_id"]]["motifs"])]

    # Keep only HYAL1 isoforms
    hyal1_uid = next((uid for uid, c in canonicals.items() if c["gene"] == GENE), None)
    if hyal1_uid is None:
        print(f"ERROR: {GENE} not found in canonical set")
        return

    isoforms = [iso for iso in filtered if iso["uniprot_id"] == hyal1_uid]
    print(f"{GENE} ({hyal1_uid}): {len(isoforms)} isoforms after strict filter")

    if not isoforms:
        print("No isoforms to plot.")
        return

    # Build matrix
    combined_map = {iso["isoform_id"]: classify_combined(iso, canonicals)
                    for iso in isoforms}
    kp      = len(canonicals[hyal1_uid]["motifs"])
    iso_ids = [iso["isoform_id"] for iso in isoforms]
    n       = len(iso_ids)

    mat = np.full((n, N_POS), np.nan)
    for row, iso_id in enumerate(iso_ids):
        for col in range(min(kp, N_POS)):
            s = combined_map[iso_id][col]
            mat[row, col] = 0.0 if s == "intact" else 1.0

    # Cluster rows
    mat_filled = np.where(np.isnan(mat), 0.0, mat)
    order   = leaves_list(linkage(pdist(mat_filled, "euclidean"), method="ward")) if n > 1 else [0]
    mat     = mat[order]
    iso_ids = [iso_ids[i] for i in order]

    # Y-axis labels: isoform number only (e.g. "HYAL1-3")
    labels = []
    seen   = {}
    for iso_id in iso_ids:
        seen[GENE] = seen.get(GENE, 0) + 1
        labels.append(f"{GENE}-{seen[GENE]}")

    # Plot
    fig_h = max(3, n * 0.45)
    fig, ax = plt.subplots(figsize=(7, fig_h))

    cmap = matplotlib.colors.ListedColormap(["#f7f7f7", "#C44E52"])
    cmap.set_bad(color="#d0d0d0")
    masked = np.ma.masked_invalid(mat)
    im = ax.imshow(masked, cmap=cmap, vmin=0, vmax=1,
                   aspect="auto", interpolation="nearest")

    ax.set_xticks(range(N_POS))
    ax.set_xticklabels([f"$\\beta\\alpha${k}" for k in range(1, N_POS + 1)], fontsize=14)
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels, fontsize=13)
    ax.set_xlabel("Barrel motif position", fontsize=15)
    ax.set_title(f"{GENE} isoform disruption ({n} isoforms)\n"
                 "red = disrupted",
                 fontsize=16)

    ax.set_xticks(np.arange(-0.5, N_POS, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.6)
    ax.tick_params(which="minor", length=0)

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_ticks([0.25, 0.75])
    cbar.set_ticklabels(["intact", "disrupted"], fontsize=12)
    cbar.ax.tick_params(length=0)

    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {OUT}")
    conn.close()


if __name__ == "__main__":
    main()
