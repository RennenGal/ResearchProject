#!/usr/bin/env python3
"""
Exon junction enrichment stratified by TIM barrel domain length.

For proteins with domain length 200-400 aa, split into 50 aa windows and
compare whether junction placement (beta / loop / alpha / between / flanking)
is consistent across subgroups.

Usage:
    python scripts/analyze_domain_length_subgroups.py
    python scripts/analyze_domain_length_subgroups.py --db db/protein_data.db
    python scripts/analyze_domain_length_subgroups.py --min 200 --max 400 --step 50
"""

import argparse
import json
import sqlite3
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from protein_data_collector.config import get_config

ORDER = ["between_motifs", "in_motif_beta", "in_motif_loop", "in_motif_alpha", "flanking"]


def classify(pos, motifs):
    if not motifs:
        return "no_motifs"
    if pos < motifs[0]["beta_start"]:
        return "before_barrel"
    for i, m in enumerate(motifs):
        n = m["motif"]
        if m["beta_start"] <= pos <= m["beta_end"]:
            return f"beta_{n}"
        if m["beta_end"] < pos < m["alpha_start"]:
            return f"loop_{n}"
        if m["alpha_start"] <= pos <= m["alpha_end"]:
            return f"alpha_{n}"
        if i + 1 < len(motifs):
            nxt = motifs[i + 1]
            if m["alpha_end"] < pos < nxt["beta_start"]:
                return f"inter_{n}_{n+1}"
    return "after_barrel"


def broad(label):
    if label.startswith("beta_"):  return "in_motif_beta"
    if label.startswith("loop_"):  return "in_motif_loop"
    if label.startswith("alpha_"): return "in_motif_alpha"
    if label.startswith("inter_"): return "between_motifs"
    return "flanking"


def element_lengths(motifs, ds, de):
    lens = {}
    if not motifs:
        return lens
    bb_end = motifs[0]["beta_start"] - 1
    if bb_end >= ds:
        lens["before_barrel"] = bb_end - ds + 1
    for i, m in enumerate(motifs):
        n = m["motif"]
        lens[f"beta_{n}"]  = m["beta_end"]    - m["beta_start"]  + 1
        lens[f"loop_{n}"]  = m["alpha_start"] - m["beta_end"]    - 1
        lens[f"alpha_{n}"] = m["alpha_end"]   - m["alpha_start"] + 1
        if i + 1 < len(motifs):
            gap = motifs[i + 1]["beta_start"] - m["alpha_end"] - 1
            if gap > 0:
                lens[f"inter_{n}_{n+1}"] = gap
    ab_start = motifs[-1]["alpha_end"] + 1
    if ab_start <= de:
        lens["after_barrel"] = de - ab_start + 1
    return lens


def analyze_bucket(bucket_rows):
    cat_counts = Counter()
    cat_res    = defaultdict(int)
    total_j    = 0
    domain_lens = []

    for uid, gene, ds, de, ea, ma, dl in bucket_rows:
        exons  = json.loads(ea)
        motifs = json.loads(ma)
        domain_lens.append(dl)

        for exon in exons[:-1]:
            pos = exon["end"]
            if ds <= pos < de:
                label = classify(pos, motifs)
                cat_counts[broad(label)] += 1
                total_j += 1

        for elem, length in element_lengths(motifs, ds, de).items():
            cat_res[broad(elem)] += length

    total_res   = sum(cat_res.values())
    exp_density = total_j / total_res if total_res > 0 else 0

    return {
        "n":          len(bucket_rows),
        "mean_len":   statistics.mean(domain_lens) if domain_lens else 0,
        "total_j":    total_j,
        "cat_counts": cat_counts,
        "cat_res":    cat_res,
        "exp_density": exp_density,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db",   default=None)
    parser.add_argument("--min",  type=int, default=200)
    parser.add_argument("--max",  type=int, default=400)
    parser.add_argument("--step", type=int, default=50)
    args = parser.parse_args()

    db_path = args.db or get_config().db_path
    conn    = sqlite3.connect(db_path)

    rows = conn.execute("""
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

    in_range = [(uid, gene, ds, de, ea, ma, dl)
                for uid, gene, ds, de, ea, ma, dl in rows
                if args.min <= dl < args.max]

    print(f"\nProteins with domain {args.min}-{args.max} aa + both exon/motif data: {len(in_range)}")

    buckets = range(args.min, args.max, args.step)
    results = []
    for lo in buckets:
        hi = lo + args.step
        bucket = [r for r in in_range if lo <= r[6] < hi]
        if not bucket:
            continue
        res = analyze_bucket(bucket)
        res["range"] = f"{lo}-{hi}"
        results.append(res)

    if not results:
        print("No data.")
        return

    col_w = 16

    # --- Enrichment table ---
    print(f"\n{'='*70}")
    print("  Enrichment by domain length subgroup  (junction density / expected)")
    print(f"{'='*70}")
    header = f"  {'Category':<18}"
    for r in results:
        label = f"{r['range']} ({r['n']}p)"
        header += f"  {label:>{col_w}}"
    print(header)
    print("  " + "-"*18 + ("  " + "-"*col_w) * len(results))
    for cat in ORDER:
        line = f"  {cat:<18}"
        for r in results:
            j       = r["cat_counts"].get(cat, 0)
            res_aa  = r["cat_res"].get(cat, 0)
            density = j / res_aa if res_aa > 0 else 0
            enrich  = density / r["exp_density"] if r["exp_density"] > 0 else 0
            line   += f"  {enrich:>{col_w-1}.2f}x"
        print(line)

    # --- Raw % table ---
    print(f"\n{'='*70}")
    print("  Raw junction distribution (% of in-domain junctions per subgroup)")
    print(f"{'='*70}")
    print(header)
    print("  " + "-"*18 + ("  " + "-"*col_w) * len(results))
    for cat in ORDER:
        line = f"  {cat:<18}"
        for r in results:
            j     = r["cat_counts"].get(cat, 0)
            total = r["total_j"]
            pct   = 100 * j / total if total > 0 else 0
            cell  = f"{j} ({pct:.1f}%)"
            line += f"  {cell:>{col_w}}"
        print(line)

    # --- Subgroup summary ---
    print(f"\n{'='*70}")
    print("  Subgroup summary")
    print(f"{'='*70}")
    print(f"  {'Range':<10}  {'N':>4}  {'Mean len':>9}  {'Junctions':>10}  {'J/protein':>10}")
    print(f"  {'-'*10}  {'-'*4}  {'-'*9}  {'-'*10}  {'-'*10}")
    for r in results:
        jp = r["total_j"] / r["n"] if r["n"] > 0 else 0
        print(f"  {r['range']:<10}  {r['n']:>4}  {r['mean_len']:>9.1f}  {r['total_j']:>10}  {jp:>10.1f}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
