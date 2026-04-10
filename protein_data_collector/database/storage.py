"""CRUD operations for the three-tier protein database."""

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
# TIM barrel entries
# ---------------------------------------------------------------------------

def upsert_tim_barrel_entry(conn: sqlite3.Connection, entry: TIMBarrelEntry) -> None:
    conn.execute(
        """
        INSERT OR REPLACE INTO tim_barrel_entries
            (accession, entry_type, name, description, tim_barrel_annotation)
        VALUES (?, ?, ?, ?, ?)
        """,
        (entry.accession, entry.entry_type, entry.name,
         entry.description, entry.tim_barrel_annotation),
    )


def get_all_tim_barrel_entries(conn: sqlite3.Connection) -> List[dict]:
    rows = conn.execute("SELECT * FROM tim_barrel_entries").fetchall()
    return [_row_to_dict(r) for r in rows]


# ---------------------------------------------------------------------------
# Proteins
# ---------------------------------------------------------------------------

def upsert_protein(conn: sqlite3.Connection, protein: Protein) -> None:
    conn.execute(
        """
        INSERT OR REPLACE INTO proteins
            (uniprot_id, tim_barrel_accession, protein_name, gene_name,
             organism, reviewed, protein_existence, annotation_score)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (protein.uniprot_id, protein.tim_barrel_accession, protein.protein_name,
         protein.gene_name, protein.organism, int(protein.reviewed) if protein.reviewed is not None else None,
         protein.protein_existence, protein.annotation_score),
    )


def get_all_proteins(conn: sqlite3.Connection) -> List[dict]:
    rows = conn.execute("SELECT * FROM proteins").fetchall()
    return [_row_to_dict(r) for r in rows]


def get_proteins_without_isoforms(conn: sqlite3.Connection) -> List[str]:
    """Return uniprot_ids that have no rows in the isoforms table."""
    rows = conn.execute(
        """
        SELECT p.uniprot_id
        FROM proteins p
        LEFT JOIN isoforms i ON p.uniprot_id = i.uniprot_id
        WHERE i.uniprot_id IS NULL
        ORDER BY p.uniprot_id
        """
    ).fetchall()
    return [r[0] for r in rows]


# ---------------------------------------------------------------------------
# Isoforms
# ---------------------------------------------------------------------------

def upsert_isoform(conn: sqlite3.Connection, isoform: Isoform) -> None:
    conn.execute(
        """
        INSERT OR REPLACE INTO isoforms
            (isoform_id, uniprot_id, is_canonical, sequence, sequence_length,
             exon_count, exon_annotations, splice_variants,
             tim_barrel_location, ensembl_gene_id, alphafold_id)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            isoform.isoform_id,
            isoform.uniprot_id,
            int(isoform.is_canonical),
            isoform.sequence,
            isoform.sequence_length,
            isoform.exon_count,
            json.dumps(isoform.exon_annotations) if isoform.exon_annotations is not None else None,
            json.dumps(isoform.splice_variants) if isoform.splice_variants is not None else None,
            json.dumps(isoform.tim_barrel_location) if isoform.tim_barrel_location is not None else None,
            isoform.ensembl_gene_id,
            isoform.alphafold_id,
        ),
    )


def get_isoforms_for_protein(conn: sqlite3.Connection, uniprot_id: str) -> List[dict]:
    rows = conn.execute(
        "SELECT * FROM isoforms WHERE uniprot_id = ? ORDER BY is_canonical DESC, isoform_id",
        (uniprot_id,),
    ).fetchall()
    return [_row_to_dict(r) for r in rows]


def get_all_isoforms(conn: sqlite3.Connection) -> List[dict]:
    rows = conn.execute(
        "SELECT * FROM isoforms ORDER BY uniprot_id, is_canonical DESC"
    ).fetchall()
    return [_row_to_dict(r) for r in rows]


# ---------------------------------------------------------------------------
# Bulk helpers
# ---------------------------------------------------------------------------

def upsert_tim_barrel_entries(conn: sqlite3.Connection, entries: List[TIMBarrelEntry]) -> None:
    for entry in entries:
        upsert_tim_barrel_entry(conn, entry)
    conn.commit()


def upsert_proteins(conn: sqlite3.Connection, proteins: List[Protein]) -> None:
    for protein in proteins:
        upsert_protein(conn, protein)
    conn.commit()


def upsert_isoforms(conn: sqlite3.Connection, isoforms: List[Isoform]) -> None:
    for isoform in isoforms:
        upsert_isoform(conn, isoform)
    conn.commit()


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def get_counts(conn: sqlite3.Connection) -> dict:
    counts = {}
    for table in ("tim_barrel_entries", "proteins", "isoforms"):
        row = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()
        counts[table] = row[0]
    row = conn.execute(
        "SELECT COUNT(*) FROM isoforms WHERE is_canonical = 0"
    ).fetchone()
    counts["alternative_isoforms"] = row[0]
    return counts
