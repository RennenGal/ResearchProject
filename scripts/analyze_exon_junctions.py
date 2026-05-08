#!/usr/bin/env python3
"""
Exon junction analysis for TIM barrel canonical isoforms (UniProt data only).

Prints three key numbers, then classifies every exon junction that falls inside
the TIM barrel domain relative to the 8 beta-alpha motif units.

Key numbers
-----------
  1. Genes encoding a TIM barrel domain (distinct gene names)
  2. Genes where AS disrupts the domain (from tb_affected_isoforms)
  3. Distribution of domain variants per gene

Junction classification
-----------------------
For each protein with both exon_annotations and motif_annotations, every
exon-exon boundary inside [domain_start, domain_end] is placed into one of:

  before_barrel         — upstream of motif 1
  beta_N                — inside the beta-strand of motif N
  loop_N                — inside the beta→alpha loop of motif N (between strand and helix)
  alpha_N               — inside the alpha-helix of motif N
  inter_N_{N+1}         — between the end of motif N and the start of motif N+1
  after_barrel          — downstream of the last annotated motif

Analysis is stratified by whether the protein has a full 8-motif annotation or
a partial one (the hypothesis is testable on either set).

Usage
-----
    python scripts/analyze_exon_junctions.py
    python scripts/analyze_exon_junctions.py --db db/protein_data.db
    python scripts/analyze_exon_junctions.py --full-only   # only 8-motif proteins
"""

import argparse
import json
import sqlite3
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from protein_data_collector.config import get_config


# ---------------------------------------------------------------------------
# Key numbers
# ---------------------------------------------------------------------------

def print_key_numbers(conn: sqlite3.Connection) -> None:
    # 1. Genes encoding TIM barrel
    total_proteins = conn.execute(
        "SELECT COUNT(*) FROM tb_canonical_analysis"
    ).fetchone()[0]
    n_genes = conn.execute(
        "SELECT COUNT(DISTINCT gene_name) FROM tb_canonical_analysis WHERE gene_name IS NOT NULL"
    ).fetchone()[0]
    no_gene = total_proteins - conn.execute(
        "SELECT COUNT(*) FROM tb_canonical_analysis WHERE gene_name IS NOT NULL"
    ).fetchone()[0]

    # 2. Genes with AS in domain
    n_affected_isoforms = conn.execute(
        "SELECT COUNT(*) FROM tb_affected_isoforms"
    ).fetchone()[0]
    n_genes_as = conn.execute("""
        SELECT COUNT(DISTINCT p.gene_name)
        FROM tb_affected_isoforms ai
        JOIN tb_proteins p ON p.uniprot_id = ai.uniprot_id
        WHERE p.gene_name IS NOT NULL
    """).fetchone()[0]
    n_proteins_as = conn.execute(
        "SELECT COUNT(DISTINCT uniprot_id) FROM tb_affected_isoforms"
    ).fetchone()[0]

    # 3. Variants per gene
    variant_counts = conn.execute("""
        SELECT p.gene_name, COUNT(*) AS n
        FROM tb_affected_isoforms ai
        JOIN tb_proteins p ON p.uniprot_id = ai.uniprot_id
        WHERE p.gene_name IS NOT NULL
        GROUP BY p.gene_name
        ORDER BY n DESC
    """).fetchall()
    dist: Counter = Counter(n for _, n in variant_counts)

    print(f"\n{'='*60}")
    print("  TIM barrel - key numbers (UniProt data, Homo sapiens)")
    print(f"{'='*60}")
    print(f"  Canonical proteins in tb_canonical_analysis : {total_proteins}")
    print(f"  Distinct gene names (Swiss-Prot / reviewed) : {n_genes}")
    print(f"  Proteins without a gene name (TrEMBL)       : {no_gene}")
    print()
    print(f"  Proteins with AS disrupting the domain      : {n_proteins_as}")
    print(f"  Genes with AS disrupting the domain         : {n_genes_as}")
    print(f"  Total AS-affected isoforms (variants)       : {n_affected_isoforms}")
    print()
    print("  Distribution of AS variants per gene:")
    print(f"    {'Variants':>8}  {'Genes':>6}")
    for n_var in sorted(dist):
        print(f"    {n_var:>8}  {dist[n_var]:>6}")
    if variant_counts:
        print(f"\n  Genes with most variants:")
        for gene, n in variant_counts[:5]:
            print(f"    {gene}: {n}")
    print(f"{'='*60}")


# ---------------------------------------------------------------------------
# Junction extraction
# ---------------------------------------------------------------------------

def junctions_in_domain(exon_annotations: list[dict],
                         domain_start: int,
                         domain_end: int) -> list[int]:
    """
    Return sorted list of exon junction positions (1-based) that fall inside
    [domain_start, domain_end].

    A junction is the boundary between exon N and exon N+1.  Its position is
    defined as the last residue of exon N (exon_N['end']).  We include junctions
    where domain_start <= junction_pos < domain_end (the terminal boundary of
    the domain itself is not a splice site).
    """
    positions = []
    for i, exon in enumerate(exon_annotations[:-1]):   # skip last exon
        j_pos = exon["end"]
        if domain_start <= j_pos < domain_end:
            positions.append(j_pos)
    return positions


# ---------------------------------------------------------------------------
# Motif classification
# ---------------------------------------------------------------------------

def classify_junction(pos: int, motifs: list[dict]) -> str:
    """
    Classify a junction at *pos* relative to *motifs*.

    Returns a string label such as 'inter_1_2', 'beta_3', 'loop_5', 'alpha_7',
    'before_barrel', or 'after_barrel'.
    """
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
        # Between this motif and the next
        if i + 1 < len(motifs):
            next_m = motifs[i + 1]
            if m["alpha_end"] < pos < next_m["beta_start"]:
                return f"inter_{n}_{n+1}"

    return "after_barrel"


def label_category(label: str) -> str:
    """Map a detailed label to a broad category."""
    if label.startswith("beta_"):
        return "in_motif_beta"
    if label.startswith("loop_"):
        return "in_motif_loop"
    if label.startswith("alpha_"):
        return "in_motif_alpha"
    if label.startswith("inter_"):
        return "between_motifs"
    return "flanking"   # before_barrel or after_barrel


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def run_analysis(conn: sqlite3.Connection, full_only: bool) -> None:
    motif_filter = "AND json_array_length(motif_annotations) = 8" if full_only else ""
    rows = conn.execute(f"""
        SELECT uniprot_id, gene_name, domain_start, domain_end,
               exon_annotations, motif_annotations
        FROM tb_canonical_analysis
        WHERE exon_annotations IS NOT NULL
          AND motif_annotations IS NOT NULL
          {motif_filter}
        ORDER BY uniprot_id
    """).fetchall()

    subset_label = "8-motif proteins only" if full_only else "all proteins with motif annotations"
    print(f"\n{'='*60}")
    print(f"  Exon junction analysis - {subset_label}")
    print(f"{'='*60}")
    print(f"  Proteins analysed: {len(rows)}")

    total_junctions     = 0
    in_domain_junctions = 0

    label_counts:    Counter = Counter()
    category_counts: Counter = Counter()
    inter_gap_counts: Counter = Counter()   # inter_N_{N+1} → count

    per_protein_in_domain: list[int] = []

    for uid, gene, ds, de, exon_json, motif_json in rows:
        exons  = json.loads(exon_json)
        motifs = json.loads(motif_json)

        all_junctions = [exon["end"] for exon in exons[:-1]]
        total_junctions += len(all_junctions)

        in_domain = junctions_in_domain(exons, ds, de)
        in_domain_junctions += len(in_domain)
        per_protein_in_domain.append(len(in_domain))

        for pos in in_domain:
            label = classify_junction(pos, motifs)
            label_counts[label] += 1
            category_counts[label_category(label)] += 1
            if label.startswith("inter_"):
                inter_gap_counts[label] += 1

    if in_domain_junctions == 0:
        print("  No junctions found in domain.")
        return

    # --- broad category summary ---
    print(f"  Total exon junctions (all proteins, any position) : {total_junctions}")
    print(f"  Junctions inside TIM barrel domain                : {in_domain_junctions}")
    if per_protein_in_domain:
        avg = sum(per_protein_in_domain) / len(per_protein_in_domain)
        med = sorted(per_protein_in_domain)[len(per_protein_in_domain) // 2]
        print(f"  Per-protein: mean={avg:.1f}  median={med}")

    print(f"\n  Broad classification of domain junctions:")
    print(f"    {'Category':<20}  {'Count':>6}  {'%':>7}")
    print(f"    {'-'*20}  {'-'*6}  {'-'*7}")
    order = ["between_motifs", "in_motif_beta", "in_motif_loop", "in_motif_alpha", "flanking"]
    for cat in order:
        n = category_counts[cat]
        pct = 100 * n / in_domain_junctions
        print(f"    {cat:<20}  {n:>6}  {pct:>6.1f}%")

    # --- inter-motif gap distribution ---
    print(f"\n  Junctions in each inter-motif gap:")
    print(f"    {'Gap':<14}  {'Count':>6}  {'%':>7}")
    print(f"    {'-'*14}  {'-'*6}  {'-'*7}")
    for gap in [f"inter_{i}_{i+1}" for i in range(1, 9)]:
        n = inter_gap_counts.get(gap, 0)
        pct = 100 * n / in_domain_junctions if in_domain_junctions else 0
        label = gap.replace("inter_", "motif ").replace("_", " > ")
        print(f"    {label:<14}  {n:>6}  {pct:>6.1f}%")

    # --- per-motif detail ---
    print(f"\n  Per-motif detail (beta / loop / alpha junctions):")
    print(f"    {'Region':<16}  {'Count':>6}  {'%':>7}")
    print(f"    {'-'*16}  {'-'*6}  {'-'*7}")
    n_motifs = 8 if full_only else max(
        (int(k.split("_")[1]) for k in label_counts
         if k.startswith(("beta_", "loop_", "alpha_"))),
        default=8
    )
    for m in range(1, n_motifs + 1):
        for prefix, tag in [("beta", "beta"), ("loop", "loop"), ("alpha", "alpha")]:
            key = f"{prefix}_{m}"
            n = label_counts.get(key, 0)
            pct = 100 * n / in_domain_junctions
            print(f"    {tag} motif {m:<8}  {n:>6}  {pct:>6.1f}%")

    print(f"{'='*60}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Exon junction statistics for TIM barrel canonical isoforms"
    )
    parser.add_argument("--db",        default=None)
    parser.add_argument("--full-only", action="store_true",
                        help="Restrict analysis to proteins with full 8-motif annotation")
    args = parser.parse_args()

    db_path = args.db or get_config().db_path
    conn = sqlite3.connect(db_path)

    print_key_numbers(conn)
    run_analysis(conn, full_only=args.full_only)

    conn.close()


if __name__ == "__main__":
    main()
