"""CRUD operations for the protein database."""

import json
import sqlite3
from typing import List, Optional

from ..models.entities import TIMBarrelEntry, Protein, Isoform


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _row_to_dict(row: sqlite3.Row) -> dict:
    return dict(row)


# ---------------------------------------------------------------------------
# Domain entries
# ---------------------------------------------------------------------------

def upsert_domain_entry(
    conn: sqlite3.Connection,
    entry: TIMBarrelEntry,
    table: str = "tb_entries",
) -> None:
    conn.execute(
        f"""
        INSERT OR REPLACE INTO {table}
            (accession, entry_type, name, description, domain_annotation)
        VALUES (?, ?, ?, ?, ?)
        """,
        (entry.accession, entry.entry_type, entry.name,
         entry.description, entry.domain_annotation),
    )


def get_all_domain_entries(
    conn: sqlite3.Connection,
    table: str = "tb_entries",
) -> List[dict]:
    rows = conn.execute(f"SELECT * FROM {table}").fetchall()
    return [_row_to_dict(r) for r in rows]


# ---------------------------------------------------------------------------
# Proteins
# ---------------------------------------------------------------------------

def upsert_protein(
    conn: sqlite3.Connection,
    protein: Protein,
    table: str = "tb_proteins",
) -> None:
    conn.execute(
        f"""
        INSERT OR REPLACE INTO {table}
            (uniprot_id, tim_barrel_accession, protein_name, gene_name,
             organism, reviewed, protein_existence, annotation_score, canonical_uniprot_id)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (protein.uniprot_id, protein.tim_barrel_accession, protein.protein_name,
         protein.gene_name, protein.organism,
         int(protein.reviewed) if protein.reviewed is not None else None,
         protein.protein_existence, protein.annotation_score, protein.canonical_uniprot_id),
    )


def get_all_proteins(
    conn: sqlite3.Connection,
    table: str = "tb_proteins",
) -> List[dict]:
    rows = conn.execute(f"SELECT * FROM {table}").fetchall()
    return [_row_to_dict(r) for r in rows]


def get_proteins_without_isoforms(
    conn: sqlite3.Connection,
    protein_table: str = "tb_proteins",
    isoform_table: str = "tb_isoforms",
) -> List[str]:
    """Return uniprot_ids of canonical proteins that have no rows in the isoforms table."""
    rows = conn.execute(
        f"""
        SELECT p.uniprot_id
        FROM {protein_table} p
        LEFT JOIN {isoform_table} i ON p.uniprot_id = i.uniprot_id
        WHERE p.canonical_uniprot_id IS NULL
          AND i.uniprot_id IS NULL
        ORDER BY p.uniprot_id
        """
    ).fetchall()
    return [r[0] for r in rows]


def deduplicate_proteins(
    conn: sqlite3.Connection,
    protein_table: str = "tb_proteins",
    isoform_table: str = "tb_isoforms",
) -> int:
    """
    Group proteins by (protein_name, organism) and mark redundant entries by setting
    canonical_uniprot_id to point to the best representative in each group.

    Ranking within a group: reviewed DESC, alt-spliced isoform count DESC,
    total isoform count DESC, annotation_score DESC, uniprot_id ASC (deterministic).

    Returns the total number of proteins marked as redundant (cumulative).
    """
    conn.execute(f"""
        UPDATE {protein_table}
        SET canonical_uniprot_id = (
            SELECT best.uniprot_id
            FROM (
                SELECT
                    p2.uniprot_id,
                    COALESCE(p2.protein_name, p2.uniprot_id) AS group_name,
                    p2.organism,
                    ROW_NUMBER() OVER (
                        PARTITION BY COALESCE(p2.protein_name, p2.uniprot_id), p2.organism
                        ORDER BY
                            p2.reviewed                                                    DESC,
                            SUM(CASE WHEN i2.is_canonical = 0 THEN 1 ELSE 0 END)          DESC,
                            COUNT(i2.isoform_id)                                           DESC,
                            p2.annotation_score                                            DESC,
                            p2.uniprot_id                                                  ASC
                    ) AS rn
                FROM {protein_table} p2
                LEFT JOIN {isoform_table} i2 ON p2.uniprot_id = i2.uniprot_id
                GROUP BY p2.uniprot_id
            ) best
            WHERE best.group_name = COALESCE({protein_table}.protein_name, {protein_table}.uniprot_id)
              AND best.organism   = {protein_table}.organism
              AND best.rn         = 1
              AND best.uniprot_id != {protein_table}.uniprot_id
        )
        WHERE canonical_uniprot_id IS NULL
    """)
    changed = conn.execute(
        f"SELECT COUNT(*) FROM {protein_table} WHERE canonical_uniprot_id IS NOT NULL"
    ).fetchone()[0]
    conn.commit()
    return changed


# ---------------------------------------------------------------------------
# Isoforms
# ---------------------------------------------------------------------------

def upsert_isoform(
    conn: sqlite3.Connection,
    isoform: Isoform,
    table: str = "tb_isoforms",
) -> None:
    conn.execute(
        f"""
        INSERT OR REPLACE INTO {table}
            (isoform_id, uniprot_id, is_canonical, sequence, sequence_length,
             is_fragment, exon_count, exon_annotations, splice_variants,
             tim_barrel_location, tim_barrel_sequence, ensembl_transcript_id, alphafold_id)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            isoform.isoform_id,
            isoform.uniprot_id,
            int(isoform.is_canonical),
            isoform.sequence,
            isoform.sequence_length,
            int(isoform.is_fragment),
            isoform.exon_count,
            json.dumps(isoform.exon_annotations) if isoform.exon_annotations is not None else None,
            json.dumps(isoform.splice_variants) if isoform.splice_variants is not None else None,
            json.dumps(isoform.tim_barrel_location) if isoform.tim_barrel_location is not None else None,
            isoform.tim_barrel_sequence,
            isoform.ensembl_transcript_id,
            isoform.alphafold_id,
        ),
    )


def get_isoforms_for_protein(
    conn: sqlite3.Connection,
    uniprot_id: str,
    table: str = "tb_isoforms",
) -> List[dict]:
    rows = conn.execute(
        f"SELECT * FROM {table} WHERE uniprot_id = ? ORDER BY is_canonical DESC, isoform_id",
        (uniprot_id,),
    ).fetchall()
    return [_row_to_dict(r) for r in rows]


def get_all_isoforms(
    conn: sqlite3.Connection,
    table: str = "tb_isoforms",
) -> List[dict]:
    rows = conn.execute(
        f"SELECT * FROM {table} ORDER BY uniprot_id, is_canonical DESC"
    ).fetchall()
    return [_row_to_dict(r) for r in rows]


# ---------------------------------------------------------------------------
# Bulk helpers
# ---------------------------------------------------------------------------

def upsert_domain_entries(
    conn: sqlite3.Connection,
    entries: List[TIMBarrelEntry],
    table: str = "tb_entries",
) -> None:
    for entry in entries:
        upsert_domain_entry(conn, entry, table=table)
    conn.commit()


def upsert_proteins(
    conn: sqlite3.Connection,
    proteins: List[Protein],
    table: str = "tb_proteins",
) -> None:
    for protein in proteins:
        upsert_protein(conn, protein, table=table)
    conn.commit()


def upsert_isoforms(
    conn: sqlite3.Connection,
    isoforms: List[Isoform],
    table: str = "tb_isoforms",
) -> None:
    for isoform in isoforms:
        upsert_isoform(conn, isoform, table=table)
    conn.commit()


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def get_counts(
    conn: sqlite3.Connection,
    entries_table: str = "tb_entries",
    protein_table: str = "tb_proteins",
    isoform_table: str = "tb_isoforms",
) -> dict:
    counts = {}
    for table in (entries_table, protein_table, isoform_table):
        row = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()
        counts[table] = row[0]
    row = conn.execute(
        f"SELECT COUNT(*) FROM {isoform_table} WHERE is_canonical = 0"
    ).fetchone()
    counts["alternative_isoforms"] = row[0]
    return counts
