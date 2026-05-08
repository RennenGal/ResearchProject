#!/usr/bin/env python3
"""
Exon junction analysis for TIM barrel canonical isoforms (UniProt data only).

Prints three key numbers, then classifies every exon junction that falls inside
the TIM barrel domain relative to the 8 beta-alpha motif units.

Key numbers
-----------
  1. Genes encoding a TIM barrel domain (distinct gene names)
  2. Genes where AS disrupts the domain (from affected_isoforms)
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
        "SELECT COUNT(*) FROM canonical_analysis"
    ).fetchone()[0]
    n_genes = conn.execute(
        "SELECT COUNT(DISTINCT gene_name) FROM canonical_analysis WHERE gene_name IS NOT NULL"
    ).fetchone()[0]
    no_gene = total_proteins - conn.execute(
        "SELECT COUNT(*) FROM canonical_analysis WHERE gene_name IS NOT NULL"
    ).fetchone()[0]

    # 2. Genes with AS in domain
    n_affected_isoforms = conn.execute(
        "SELECT COUNT(*) FROM affected_isoforms"
    ).fetchone()[0]
    n_genes_as = conn.execute("""
        SELECT COUNT(DISTINCT p.gene_name)
        FROM affected_isoforms ai
        JOIN proteins p ON p.uniprot_id = ai.uniprot_id
        WHERE p.gene_name IS NOT NULL
    """).fetchone()[0]
    n_proteins_as = conn.execute(
        "SELECT COUNT(DISTINCT uniprot_id) FROM affected_isoforms"
    ).fetchone()[0]

    # 3. Variants per gene
    variant_counts = conn.execute("""
        SELECT p.gene_name, COUNT(*) AS n
        FROM affected_isoforms ai
        JOIN proteins p ON p.uniprot_id = ai.uniprot_id
        WHERE p.gene_name IS NOT NULL
        GROUP BY p.gene_name
        ORDER BY n DESC
    """).fetchall()
    dist: Counter = Counter(n for _, n in variant_counts)

    print(f"\n{'='*60}")
    print("  TIM barrel - key numbers (UniProt data, Homo sapiens)")
    print(f"{'='*60}")
    print(f"  Canonical proteins in canonical_analysis : {total_proteins}")
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

def _element_lengths(motifs: list[dict], domain_start: int, domain_end: int) -> dict[str, int]:
    """
    Return a mapping of element label -> total residues for one protein.
    Labels match classify_junction() output: beta_N, loop_N, alpha_N, inter_N_{N+1},
    before_barrel, after_barrel.
    """
    lengths: dict[str, int] = {}
    if not motifs:
        return lengths

    # before_barrel: domain_start to first beta_start - 1
    bb_end = motifs[0]["beta_start"] - 1
    if bb_end >= domain_start:
        lengths["before_barrel"] = bb_end - domain_start + 1

    for i, m in enumerate(motifs):
        n = m["motif"]
        lengths[f"beta_{n}"]  = m["beta_end"]   - m["beta_start"]  + 1
        lengths[f"loop_{n}"]  = m["alpha_start"] - m["beta_end"]    - 1
        lengths[f"alpha_{n}"] = m["alpha_end"]   - m["alpha_start"] + 1
        if i + 1 < len(motifs):
            next_m = motifs[i + 1]
            gap = next_m["beta_start"] - m["alpha_end"] - 1
            if gap > 0:
                lengths[f"inter_{n}_{n+1}"] = gap

    # after_barrel: last alpha_end + 1 to domain_end
    ab_start = motifs[-1]["alpha_end"] + 1
    if ab_start <= domain_end:
        lengths["after_barrel"] = domain_end - ab_start + 1

    return lengths


def run_analysis(conn: sqlite3.Connection, full_only: bool) -> None:
    motif_filter = "AND json_array_length(motif_annotations) = 8" if full_only else ""
    rows = conn.execute(f"""
        SELECT uniprot_id, gene_name, domain_start, domain_end,
               exon_annotations, motif_annotations
        FROM canonical_analysis
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
    inter_gap_counts: Counter = Counter()

    # For length-normalized density
    element_residues: Counter = Counter()   # label -> total residues across all proteins
    category_residues: defaultdict = defaultdict(int)

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

        for elem, length in _element_lengths(motifs, ds, de).items():
            element_residues[elem] += length
            category_residues[label_category(elem)] += length

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

    print(f"\n  Broad classification of domain junctions (raw counts):")
    print(f"    {'Category':<20}  {'Count':>6}  {'%':>7}")
    print(f"    {'-'*20}  {'-'*6}  {'-'*7}")
    order = ["between_motifs", "in_motif_beta", "in_motif_loop", "in_motif_alpha", "flanking"]
    for cat in order:
        n = category_counts[cat]
        pct = 100 * n / in_domain_junctions
        print(f"    {cat:<20}  {n:>6}  {pct:>6.1f}%")

    # --- length-normalized density ---
    total_residues = sum(category_residues.values())
    expected_density = in_domain_junctions / total_residues if total_residues > 0 else 0

    print(f"\n  Length-normalized junction density (junctions per residue):")
    print(f"  (Enrichment = observed density / expected uniform density)")
    print(f"    {'Category':<20}  {'Res':>7}  {'J/res':>8}  {'Enrich':>8}")
    print(f"    {'-'*20}  {'-'*7}  {'-'*8}  {'-'*8}")
    for cat in order:
        res = category_residues[cat]
        jct = category_counts[cat]
        density = jct / res if res > 0 else 0
        enrich  = density / expected_density if expected_density > 0 else 0
        print(f"    {cat:<20}  {res:>7}  {density:>8.4f}  {enrich:>8.2f}x")

    # --- inter-motif gap distribution ---
    print(f"\n  Junctions in each inter-motif gap:")
    print(f"    {'Gap':<14}  {'Count':>6}  {'%':>7}  {'Res':>6}  {'J/res':>8}  {'Enrich':>8}")
    print(f"    {'-'*14}  {'-'*6}  {'-'*7}  {'-'*6}  {'-'*8}  {'-'*8}")
    for i in range(1, 8):
        gap = f"inter_{i}_{i+1}"
        n   = inter_gap_counts.get(gap, 0)
        pct = 100 * n / in_domain_junctions if in_domain_junctions else 0
        res = element_residues.get(gap, 0)
        density = n / res if res > 0 else 0
        enrich  = density / expected_density if expected_density > 0 else 0
        lbl = f"motif {i} > {i+1}"
        print(f"    {lbl:<14}  {n:>6}  {pct:>6.1f}%  {res:>6}  {density:>8.4f}  {enrich:>8.2f}x")

    # --- per-motif detail ---
    print(f"\n  Per-motif detail (beta / loop / alpha junctions):")
    print(f"    {'Region':<16}  {'Count':>6}  {'%':>7}  {'Res':>6}  {'Enrich':>8}")
    print(f"    {'-'*16}  {'-'*6}  {'-'*7}  {'-'*6}  {'-'*8}")
    n_motifs = 8 if full_only else max(
        (int(k.split("_")[1]) for k in label_counts
         if k.startswith(("beta_", "loop_", "alpha_"))),
        default=8
    )
    for m in range(1, n_motifs + 1):
        for prefix, tag in [("beta", "beta"), ("loop", "loop"), ("alpha", "alpha")]:
            key = f"{prefix}_{m}"
            n   = label_counts.get(key, 0)
            pct = 100 * n / in_domain_junctions
            res = element_residues.get(key, 0)
            density = n / res if res > 0 else 0
            enrich  = density / expected_density if expected_density > 0 else 0
            print(f"    {tag} motif {m:<8}  {n:>6}  {pct:>6.1f}%  {res:>6}  {enrich:>8.2f}x")

    print(f"{'='*60}")


# ---------------------------------------------------------------------------
# Q2: which canonical exon junctions are exploited by AS?
# ---------------------------------------------------------------------------

def run_as_junction_analysis(conn: sqlite3.Connection, full_only: bool) -> None:
    """
    Question 2: where do AS events hit the domain structure?

    For each VSP domain event, find the canonical exon junctions that fall
    within [can_start, can_end] — the boundaries of the splice event in
    canonical protein coordinates.  Those junctions are the actual splice
    sites exploited by AS.  Classify them and compare to Q1 (all domain
    junctions) to find which elements are over-represented.
    """
    motif_filter_join  = "AND json_array_length(ca.motif_annotations) = 8" if full_only else ""
    motif_filter_plain = "AND json_array_length(motif_annotations) = 8"    if full_only else ""

    # AS-affected isoforms with VSP spans + canonical exon/motif structure
    as_rows = conn.execute(f"""
        SELECT ai.isoform_id, ai.uniprot_id, ai.vsp_domain_events,
               ca.exon_annotations, ca.motif_annotations,
               ca.domain_start, ca.domain_end
        FROM affected_isoforms ai
        JOIN canonical_analysis ca ON ca.uniprot_id = ai.uniprot_id
        WHERE ai.vsp_domain_events IS NOT NULL
          AND ca.exon_annotations IS NOT NULL
          AND ca.motif_annotations IS NOT NULL
          {motif_filter_join}
        ORDER BY ai.uniprot_id, ai.isoform_id
    """).fetchall()

    # Q1 baseline: all domain junctions across all annotated proteins
    q1_rows = conn.execute(f"""
        SELECT domain_start, domain_end, exon_annotations, motif_annotations
        FROM canonical_analysis
        WHERE exon_annotations IS NOT NULL
          AND motif_annotations IS NOT NULL
          {motif_filter_plain}
    """).fetchall()

    subset_label = "8-motif proteins only" if full_only else "all proteins with motif annotations"
    print(f"\n{'='*60}")
    print(f"  Q2: AS-exploited exon junctions - {subset_label}")
    print(f"{'='*60}")
    print(f"  AS-affected isoforms: {len(as_rows)}")

    # --- Build Q1 baseline distribution ---
    q1_cats: Counter = Counter()
    q1_total = 0
    for ds, de, exon_json, motif_json in q1_rows:
        exons  = json.loads(exon_json)
        motifs = json.loads(motif_json)
        for pos in junctions_in_domain(exons, ds, de):
            q1_cats[label_category(classify_junction(pos, motifs))] += 1
            q1_total += 1

    # --- Find canonical junctions that fall within each VSP span ---
    q2_cats:   Counter = Counter()
    q2_labels: Counter = Counter()
    q2_total  = 0
    n_isoforms_with_junctions = 0

    for iso_id, uid, vsp_json, exon_json, motif_json, ds, de in as_rows:
        vsps   = json.loads(vsp_json)
        exons  = json.loads(exon_json)
        motifs = json.loads(motif_json)

        # All canonical exon junctions inside the domain for this protein
        domain_junctions = junctions_in_domain(exons, ds, de)

        iso_junctions_found = False
        for v in vsps:
            vsp_start = v["can_start"]
            vsp_end   = v["can_end"]
            # Junctions inside this VSP's canonical span
            for pos in domain_junctions:
                if vsp_start <= pos <= vsp_end:
                    lbl = classify_junction(pos, motifs)
                    q2_cats[label_category(lbl)] += 1
                    q2_labels[lbl] += 1
                    q2_total += 1
                    iso_junctions_found = True

        if iso_junctions_found:
            n_isoforms_with_junctions += 1

    print(f"  Isoforms where VSP span contains >= 1 exon junction: {n_isoforms_with_junctions}")
    print(f"  Total AS-exploited junction instances: {q2_total}")
    print(f"  (Q1 baseline: {q1_total} domain junctions across all annotated proteins)")

    if q2_total == 0:
        print("  No canonical junctions found inside VSP spans.")
        return

    order = ["between_motifs", "in_motif_beta", "in_motif_loop", "in_motif_alpha", "flanking"]

    print(f"\n  Classification of AS-exploited junctions vs. Q1 baseline:")
    print(f"    {'Category':<20}  {'Q2 n':>6}  {'Q2 %':>7}  {'Q1 %':>7}  {'AS enrich':>10}")
    print(f"    {'-'*20}  {'-'*6}  {'-'*7}  {'-'*7}  {'-'*10}")
    for cat in order:
        q2_n   = q2_cats[cat]
        q2_pct = 100 * q2_n   / q2_total  if q2_total  > 0 else 0
        q1_pct = 100 * q1_cats[cat] / q1_total if q1_total > 0 else 0
        enrich = (q2_pct / q1_pct) if q1_pct > 0 else 0
        print(f"    {cat:<20}  {q2_n:>6}  {q2_pct:>6.1f}%  {q1_pct:>6.1f}%  {enrich:>9.2f}x")

    print(f"\n  Detailed Q2 junction placement (top elements):")
    print(f"    {'Element':<18}  {'Count':>6}  {'%':>7}")
    print(f"    {'-'*18}  {'-'*6}  {'-'*7}")
    for lbl, n in q2_labels.most_common(15):
        pct = 100 * n / q2_total
        print(f"    {lbl:<18}  {n:>6}  {pct:>6.1f}%")

    print(f"{'='*60}")


# ---------------------------------------------------------------------------
# VSP boundary placement
# ---------------------------------------------------------------------------

def run_vsp_analysis(conn: sqlite3.Connection, full_only: bool) -> None:
    """
    For each AS-affected isoform, classify where each VSP boundary (start and
    end of the splice event in canonical coordinates) falls within the motif
    structure.  This directly mirrors the paper's analysis of which structural
    elements are disrupted by alternative splicing.
    """
    motif_filter = "AND json_array_length(ca.motif_annotations) = 8" if full_only else ""
    rows = conn.execute(f"""
        SELECT ai.isoform_id, ai.uniprot_id, ai.vsp_domain_events,
               ca.motif_annotations, ca.domain_start, ca.domain_end
        FROM affected_isoforms ai
        JOIN canonical_analysis ca ON ca.uniprot_id = ai.uniprot_id
        WHERE ai.vsp_domain_events IS NOT NULL
          AND ca.motif_annotations IS NOT NULL
          {motif_filter}
        ORDER BY ai.uniprot_id, ai.isoform_id
    """).fetchall()

    subset_label = "8-motif proteins only" if full_only else "all proteins with motif annotations"
    print(f"\n{'='*60}")
    print(f"  VSP boundary placement - {subset_label}")
    print(f"{'='*60}")
    print(f"  AS-affected isoforms analysed: {len(rows)}")

    start_labels: Counter = Counter()
    end_labels:   Counter = Counter()
    start_cats:   Counter = Counter()
    end_cats:     Counter = Counter()
    n_events = 0

    for iso_id, uid, vsp_json, motif_json, ds, de in rows:
        vsps   = json.loads(vsp_json)
        motifs = json.loads(motif_json)
        for v in vsps:
            n_events += 1
            s_lbl = classify_junction(v["can_start"], motifs)
            e_lbl = classify_junction(v["can_end"],   motifs)
            start_labels[s_lbl] += 1
            end_labels[e_lbl]   += 1
            start_cats[label_category(s_lbl)] += 1
            end_cats[label_category(e_lbl)]   += 1

    if n_events == 0:
        print("  No VSP events found.")
        return

    print(f"  Total VSP domain events: {n_events}")
    print(f"\n  Where VSP events START (canonical position):")
    print(f"    {'Category':<20}  {'Count':>6}  {'%':>7}")
    print(f"    {'-'*20}  {'-'*6}  {'-'*7}")
    order = ["between_motifs", "in_motif_beta", "in_motif_loop", "in_motif_alpha", "flanking"]
    for cat in order:
        n   = start_cats[cat]
        pct = 100 * n / n_events
        print(f"    {cat:<20}  {n:>6}  {pct:>6.1f}%")

    print(f"\n  Where VSP events END (canonical position):")
    print(f"    {'Category':<20}  {'Count':>6}  {'%':>7}")
    print(f"    {'-'*20}  {'-'*6}  {'-'*7}")
    for cat in order:
        n   = end_cats[cat]
        pct = 100 * n / n_events
        print(f"    {cat:<20}  {n:>6}  {pct:>6.1f}%")

    # --- combined: either boundary falls in category ---
    combined: Counter = Counter()
    for cat, n in start_cats.items():
        combined[cat] += n
    for cat, n in end_cats.items():
        combined[cat] += n
    total_boundaries = n_events * 2

    print(f"\n  Both boundaries combined ({total_boundaries} boundary positions):")
    print(f"    {'Category':<20}  {'Count':>6}  {'%':>7}")
    print(f"    {'-'*20}  {'-'*6}  {'-'*7}")
    for cat in order:
        n   = combined[cat]
        pct = 100 * n / total_boundaries
        print(f"    {cat:<20}  {n:>6}  {pct:>6.1f}%")

    # --- per-motif element breakdown (start boundaries) ---
    print(f"\n  Detailed start-boundary placement (top elements):")
    print(f"    {'Element':<18}  {'Count':>6}  {'%':>7}")
    print(f"    {'-'*18}  {'-'*6}  {'-'*7}")
    for lbl, n in start_labels.most_common(15):
        pct = 100 * n / n_events
        print(f"    {lbl:<18}  {n:>6}  {pct:>6.1f}%")

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
    run_as_junction_analysis(conn, full_only=args.full_only)
    run_vsp_analysis(conn, full_only=args.full_only)

    conn.close()


if __name__ == "__main__":
    main()
