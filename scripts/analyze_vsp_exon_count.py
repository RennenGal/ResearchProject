#!/usr/bin/env python3
"""
Distribution of canonical domain exons affected per VSP event.

For each VSP (variant splice protein) domain event, counts the number of
canonical domain exons whose span overlaps with the VSP region [can_start, can_end].
A domain exon is any exon that overlaps the TIM-barrel domain [ds, de].

Output:
  figures/vsp_exon_count.png

Usage:
    python scripts/analyze_vsp_exon_count.py
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

sys.path.insert(0, str(Path(__file__).parent.parent))
from protein_data_collector.config import get_config


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

def load_proteins(conn):
    rows = conn.execute("""
        SELECT uniprot_id, domain_start, domain_end, exon_annotations
        FROM   view_canonical
    """).fetchall()
    proteins = {}
    for uid, ds, de, ea in rows:
        proteins[uid] = dict(ds=ds, de=de, exons=json.loads(ea))
    return proteins


def load_vsps(conn, proteins):
    rows = conn.execute("""
        SELECT uniprot_id, isoform_id, vsp_domain_events
        FROM   view_noncanonical
    """).fetchall()
    records = []
    for uid, iso_id, vsp_json in rows:
        if uid not in proteins:
            continue
        for vsp in json.loads(vsp_json):
            can_start = vsp.get("can_start")
            can_end   = vsp.get("can_end")
            if can_start is None or can_end is None:
                continue
            records.append(dict(uid=uid, isoform_id=iso_id,
                                can_start=can_start, can_end=can_end))
    return records


# ---------------------------------------------------------------------------
# Count domain exons per VSP
# ---------------------------------------------------------------------------

def count_domain_exons(vsp, protein):
    ds, de   = protein["ds"], protein["de"]
    can_start = vsp["can_start"]
    can_end   = vsp["can_end"]
    count = 0
    for exon in protein["exons"]:
        es, ee = exon["start"], exon["end"]
        if ee < ds or es > de:
            continue          # exon completely outside domain
        if es <= can_end and ee >= can_start:
            count += 1        # domain exon overlaps VSP span
    return count


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def plot_histogram(counts, n_vsps, n_isoforms, out):
    max_count = max(counts)
    bins      = np.arange(0, max_count + 2) - 0.5

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(counts, bins=bins, color="#4C72B0", alpha=0.85,
            edgecolor="white", zorder=3)

    ax.set_xlabel("Domain exons overlapping VSP span", fontsize=10)
    ax.set_ylabel("Number of VSP events", fontsize=10)
    ax.set_xticks(range(0, max_count + 1))
    ax.set_title(
        f"Canonical domain exons affected per VSP event\n"
        f"{n_vsps} VSP events across {n_isoforms} isoforms",
        fontsize=9,
    )
    ax.yaxis.grid(True, linestyle=":", alpha=0.4, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Annotate mean
    mean_val = np.mean(counts)
    ax.axvline(mean_val, color="#C44E52", lw=1.5, ls="--",
               label=f"Mean = {mean_val:.2f}")
    ax.legend(fontsize=9, frameon=False)

    fig.tight_layout()
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db",  default=None)
    parser.add_argument("--out", default="figures/vsp_exon_count.png")
    args = parser.parse_args()

    db_path  = args.db or get_config().db_path
    conn     = sqlite3.connect(db_path)
    proteins = load_proteins(conn)
    vsps     = load_vsps(conn, proteins)
    conn.close()

    print(f"Loaded {len(proteins)} canonical proteins.")
    print(f"Loaded {len(vsps)} VSP events.")

    counts = [count_domain_exons(v, proteins[v["uid"]]) for v in vsps]

    dist = Counter(counts)
    print(f"\nDomain exons per VSP event:")
    for k in sorted(dist):
        print(f"  {k} exon(s): {dist[k]} VSP events")
    print(f"\nMean:   {np.mean(counts):.2f}")
    print(f"Median: {np.median(counts):.1f}")

    n_isoforms = len({v["isoform_id"] for v in vsps})
    plot_histogram(counts, len(vsps), n_isoforms, args.out)


if __name__ == "__main__":
    main()
