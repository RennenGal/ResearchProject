#!/usr/bin/env python3
"""
Build the analysis_proteins table and associated views in the project database.

Table includes every canonical protein in the statistical analysis set
(those with >= 1 motif) plus every isoform from affected_isoforms whose
canonical protein is in that set.

Views created
-------------
  view_canonical               — canonical proteins only
  view_noncanonical            — isoforms only
  view_noncanonical_with_canonical — isoforms joined with their canonical partner
  view_N_variants              — canonical proteins with exactly N affected isoforms
                                  (created for each distinct count that exists)

Summary printed to stdout:
  - count of canonical / isoform rows
  - how many canonicals have >= 1 AS isoform in domain
  - total isoform count and VSP-event count
  - motif-count distribution across canonical proteins

Usage:
    python scripts/create_analysis_table.py
    python scripts/create_analysis_table.py --db path/to/db.sqlite
"""

import argparse
import json
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from protein_data_collector.config import get_config


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _has_pdb(pdb_source):
    """Return (has_pdb, pdb_accession).  pdb_source looks like '6VGO_A_1.82A' or 'no_structure'."""
    if not pdb_source or pdb_source in ("no_structure", "dssp_failed", "download_failed"):
        return 0, None
    return 1, pdb_source.split("_")[0]   # e.g. '6VGO'


def _count_domain_junctions(ea_json, ds, de):
    exons = json.loads(ea_json)
    return sum(1 for e in exons[:-1] if ds <= e["end"] < de)


# ---------------------------------------------------------------------------
# Build rows
# ---------------------------------------------------------------------------

def build_canonical_rows(conn):
    """Return dict uid -> row for all canonical proteins in the analysis set."""
    rows = conn.execute("""
        SELECT ca.uniprot_id, ca.gene_name, ca.domain_start, ca.domain_end,
               ca.exon_annotations, ca.motif_annotations, ca.pdb_source,
               p.protein_name
        FROM   canonical_analysis ca
        LEFT JOIN proteins p ON p.uniprot_id = ca.uniprot_id
        WHERE  ca.exon_annotations  IS NOT NULL
          AND  ca.motif_annotations IS NOT NULL
          AND  ca.domain_start      IS NOT NULL
          AND  ca.domain_end        IS NOT NULL
    """).fetchall()

    canonical = {}
    for uid, gene, ds, de, ea, ma, pdb_src, pname in rows:
        motifs = json.loads(ma)
        if not motifs:
            continue
        n_junctions = _count_domain_junctions(ea, ds, de)
        has_pdb, pdb_acc = _has_pdb(pdb_src)
        canonical[uid] = dict(
            uniprot_id                 = uid,
            isoform_id                 = None,
            gene_name                  = gene or "",
            protein_name               = pname or "",
            is_canonical               = 1,
            canonical_id               = None,
            has_experimental_structure = has_pdb,
            pdb_accession              = pdb_acc,
            num_motifs                 = len(motifs),
            num_domain_junctions       = n_junctions,
            num_as_variants            = 0,   # filled in below
            identity_pct               = None,
            num_vsp_domain_events      = None,
            domain_start               = ds,
            domain_end                 = de,
            exon_annotations           = ea,
            motif_annotations          = ma,
            vsp_domain_events          = None,
        )
    return canonical


def build_isoform_rows(conn, canonical_uids):
    """Return list of dicts for affected_isoforms whose canonical protein is in canonical_uids."""
    rows = conn.execute("""
        SELECT ai.isoform_id, ai.uniprot_id, ai.gene_name,
               ai.identity_percentage, ai.vsp_domain_events
        FROM   affected_isoforms ai
        WHERE  ai.vsp_domain_events IS NOT NULL
    """).fetchall()

    isoforms = []
    for iso_id, uid, gene, identity_pct, vsp_json in rows:
        if uid not in canonical_uids:
            continue
        vsps = json.loads(vsp_json) if vsp_json else []
        if not vsps:
            continue
        isoforms.append(dict(
            uniprot_id                 = uid,
            isoform_id                 = iso_id,
            gene_name                  = gene or "",
            protein_name               = "",
            is_canonical               = 0,
            canonical_id               = uid,
            has_experimental_structure = None,
            pdb_accession              = None,
            num_motifs                 = None,
            num_domain_junctions       = None,
            num_as_variants            = None,
            identity_pct               = identity_pct,
            num_vsp_domain_events      = len(vsps),
            domain_start               = None,
            domain_end                 = None,
            exon_annotations           = None,
            motif_annotations          = None,
            vsp_domain_events          = vsp_json,
        ))
    return isoforms


# ---------------------------------------------------------------------------
# Create table and views
# ---------------------------------------------------------------------------

CREATE_TABLE = """
CREATE TABLE IF NOT EXISTS analysis_proteins (
    row_id                     INTEGER PRIMARY KEY AUTOINCREMENT,
    uniprot_id                 TEXT    NOT NULL,
    isoform_id                 TEXT,
    gene_name                  TEXT,
    protein_name               TEXT,
    is_canonical               INTEGER NOT NULL,
    canonical_id               TEXT,
    has_experimental_structure INTEGER,
    pdb_accession              TEXT,
    num_motifs                 INTEGER,
    num_domain_junctions       INTEGER,
    num_as_variants            INTEGER,
    identity_pct               REAL,
    num_vsp_domain_events      INTEGER,
    domain_start               INTEGER,
    domain_end                 INTEGER,
    exon_annotations           TEXT,
    motif_annotations          TEXT,
    vsp_domain_events          TEXT,
    UNIQUE(uniprot_id, isoform_id)
)
"""

VIEWS = {
    "view_canonical": """
        SELECT * FROM analysis_proteins WHERE is_canonical = 1
    """,
    "view_noncanonical": """
        SELECT * FROM analysis_proteins WHERE is_canonical = 0
    """,
    "view_noncanonical_with_canonical": """
        SELECT
            c.uniprot_id,
            NULL               AS isoform_id,
            c.gene_name,
            c.protein_name,
            1                  AS canonical,
            c.has_experimental_structure AS has_pdb,
            c.pdb_accession,
            c.num_motifs,
            c.num_domain_junctions,
            c.num_as_variants,
            NULL               AS identity_pct,
            NULL               AS num_vsp_domain_events
        FROM analysis_proteins c
        WHERE c.is_canonical = 1
          AND c.uniprot_id IN (
              SELECT canonical_id FROM analysis_proteins WHERE is_canonical = 0
          )

        UNION ALL

        SELECT
            nc.canonical_id    AS uniprot_id,
            nc.isoform_id,
            nc.gene_name,
            c.protein_name,
            0                  AS canonical,
            c.has_experimental_structure AS has_pdb,
            c.pdb_accession,
            c.num_motifs,
            c.num_domain_junctions,
            c.num_as_variants,
            nc.identity_pct,
            nc.num_vsp_domain_events
        FROM analysis_proteins nc
        JOIN analysis_proteins c
          ON nc.canonical_id = c.uniprot_id AND c.is_canonical = 1
        WHERE nc.is_canonical = 0
    """,
}


def create_variant_views(conn, counts):
    """Create one view per distinct num_as_variants value for canonical proteins."""
    for k in sorted(counts):
        vname = f"view_{k}_variant{'s' if k != 1 else ''}"
        conn.execute(f"DROP VIEW IF EXISTS {vname}")
        conn.execute(f"""
            CREATE VIEW {vname} AS
            SELECT * FROM analysis_proteins
            WHERE is_canonical = 1 AND num_as_variants = {k}
        """)
        print(f"  Created view {vname}  ({counts[k]} proteins)")


def setup_db(conn, canonical_rows, isoform_rows):
    conn.execute("DROP TABLE IF EXISTS analysis_proteins")
    # Drop all dependent views first
    view_names = [r[0] for r in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='view'"
    ).fetchall()]
    for vn in view_names:
        conn.execute(f"DROP VIEW IF EXISTS [{vn}]")

    conn.execute(CREATE_TABLE)

    cols = [
        "uniprot_id", "isoform_id", "gene_name", "protein_name",
        "is_canonical", "canonical_id", "has_experimental_structure",
        "pdb_accession", "num_motifs", "num_domain_junctions",
        "num_as_variants", "identity_pct", "num_vsp_domain_events",
        "domain_start", "domain_end", "exon_annotations",
        "motif_annotations", "vsp_domain_events",
    ]
    placeholders = ", ".join("?" for _ in cols)
    col_list     = ", ".join(cols)
    sql          = f"INSERT OR IGNORE INTO analysis_proteins ({col_list}) VALUES ({placeholders})"

    all_rows = list(canonical_rows.values()) + isoform_rows
    conn.executemany(sql, [[r[c] for c in cols] for r in all_rows])

    # Standard views
    for vname, vbody in VIEWS.items():
        conn.execute(f"DROP VIEW IF EXISTS [{vname}]")
        conn.execute(f"CREATE VIEW {vname} AS {vbody}")
        print(f"  Created view {vname}")

    conn.commit()


# ---------------------------------------------------------------------------
# Summary statistics and motif distribution
# ---------------------------------------------------------------------------

def print_summary(canonical_rows, isoform_rows):
    n_canonical = len(canonical_rows)
    n_isoforms  = len(isoform_rows)
    n_with_iso  = sum(1 for r in canonical_rows.values() if r["num_as_variants"] > 0)
    total_vsps  = sum(r["num_vsp_domain_events"] for r in isoform_rows)
    n_with_pdb  = sum(1 for r in canonical_rows.values() if r["has_experimental_structure"])

    print("\n" + "=" * 60)
    print("  Analysis set summary")
    print("=" * 60)
    print(f"  Canonical proteins in analysis      : {n_canonical}")
    print(f"    with experimental PDB structure   : {n_with_pdb}")
    print(f"    with AlphaFold only               : {n_canonical - n_with_pdb}")
    print(f"  Canonicals with >= 1 domain isoform : {n_with_iso}")
    print(f"  Isoforms (affected) in analysis     : {n_isoforms}")
    print(f"  Total VSP domain events             : {total_vsps}")
    print()

    # Motif distribution
    motif_counts = Counter(r["num_motifs"] for r in canonical_rows.values())
    print("  Motif-count distribution (canonical proteins):")
    print(f"  {'Motifs':>6}  {'Proteins':>8}  {'Bar':}")
    for k in sorted(motif_counts):
        bar = "#" * motif_counts[k]
        print(f"  {k:>6}  {motif_counts[k]:>8}  {bar}")
    print()

    # Variant distribution
    variant_counts = Counter(r["num_as_variants"] for r in canonical_rows.values())
    print("  AS-variant distribution (per canonical protein):")
    print(f"  {'Variants':>8}  {'Proteins':>8}  {'Bar':}")
    for k in sorted(variant_counts):
        bar = "#" * min(variant_counts[k], 50)
        print(f"  {k:>8}  {variant_counts[k]:>8}  {bar}")
    print("=" * 60)

    return motif_counts, variant_counts


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", default=None)
    args = parser.parse_args()

    db_path = args.db or get_config().db_path
    conn    = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row

    print("Building canonical rows …")
    canonical_rows = build_canonical_rows(conn)
    print(f"  {len(canonical_rows)} canonical proteins")

    print("Building isoform rows …")
    isoform_rows = build_isoform_rows(conn, set(canonical_rows))
    print(f"  {len(isoform_rows)} isoforms with VSP domain events")

    # Backfill num_as_variants into canonical rows
    variant_counts_by_uid = defaultdict(int)
    for r in isoform_rows:
        variant_counts_by_uid[r["canonical_id"]] += 1
    for uid, row in canonical_rows.items():
        row["num_as_variants"] = variant_counts_by_uid[uid]

    print("\nSetting up analysis_proteins table and views …")
    conn.row_factory = None
    setup_db(conn, canonical_rows, isoform_rows)

    # Variant-count views
    variant_dist = Counter(r["num_as_variants"] for r in canonical_rows.values())
    create_variant_views(conn, variant_dist)
    conn.commit()

    print_summary(canonical_rows, isoform_rows)

    conn.close()
    print(f"Done.  DB: {db_path}")


if __name__ == "__main__":
    main()
