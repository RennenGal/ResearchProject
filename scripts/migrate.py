#!/usr/bin/env python3
"""
Migrate the existing database to the new schema.

Old schema  (produced by the collection scripts)
    proteins      — one row per canonical protein, isoforms stored as JSON
    interpro_proteins — protein ↔ TIM-barrel family mapping
    tim_barrel_entries — unchanged

New schema
    tim_barrel_entries — same structure, trimmed extra columns
    proteins           — one row per UniProt protein (was interpro_proteins)
    isoforms           — one row per isoform (was the JSON blob in old proteins)

Run:
    python scripts/migrate.py [--old db/protein_data.db] [--new db/protein_data_new.db]

The old database is left untouched.  After verifying the new one, swap it in manually.
"""

import argparse
import json
import logging
import sqlite3
import sys
from pathlib import Path

# Allow running from the repo root
sys.path.insert(0, str(Path(__file__).parent.parent))

from protein_data_collector.database.connection import ensure_db, get_connection
from protein_data_collector.database.storage import (
    upsert_domain_entries, upsert_isoforms, upsert_proteins,
)
from protein_data_collector.models.entities import Isoform, Protein, TIMBarrelEntry

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def migrate(old_path: str, new_path: str) -> None:
    logger.info("Migrating %s → %s", old_path, new_path)
    ensure_db(new_path)

    old = sqlite3.connect(old_path)
    old.row_factory = sqlite3.Row

    _migrate_tim_barrel_entries(old, new_path)
    _migrate_proteins(old, new_path)
    _migrate_isoforms(old, new_path)

    old.close()
    _print_counts(new_path)


# ---------------------------------------------------------------------------
# Phase helpers
# ---------------------------------------------------------------------------

def _migrate_tim_barrel_entries(old: sqlite3.Connection, new_path: str) -> None:
    rows = old.execute("SELECT * FROM tim_barrel_entries").fetchall()
    entries = []
    for r in rows:
        try:
            entries.append(TIMBarrelEntry(
                accession=r["accession"],
                entry_type=r["entry_type"],
                name=r["name"],
                description=r["description"],
                domain_annotation=r.get("domain_annotation") or r.get("tim_barrel_annotation") or "TIM barrel",
            ))
        except Exception as e:
            logger.warning("Skipping tim_barrel_entry %s: %s", r["accession"], e)

    with get_connection(new_path) as conn:
        upsert_domain_entries(conn, entries)
    logger.info("Migrated %d tim_barrel_entries", len(entries))


def _migrate_proteins(old: sqlite3.Connection, new_path: str) -> None:
    """Create proteins rows from interpro_proteins + protein metadata."""
    ip_rows = old.execute("SELECT * FROM interpro_proteins").fetchall()
    # Build a lookup of protein metadata from the old proteins table
    meta = {}
    try:
        p_rows = old.execute(
            "SELECT uniprot_id, protein_name, organism, reviewed, protein_existence, annotation_score FROM proteins"
        ).fetchall()
        for r in p_rows:
            meta[r["uniprot_id"]] = dict(r)
    except sqlite3.OperationalError:
        pass  # old proteins table might not exist yet

    proteins = []
    for r in ip_rows:
        uid = r["uniprot_id"]
        m = meta.get(uid, {})
        # Extract gene_name from protein_name heuristic (entry name like "P53_HUMAN")
        old_name = old.execute(
            "SELECT name FROM proteins WHERE uniprot_id = ?", (uid,)
        ).fetchone()
        gene_name = None
        if old_name and old_name[0]:
            gene_name = old_name[0].split("_")[0]  # e.g. "P53_HUMAN" → "P53"

        try:
            proteins.append(Protein(
                uniprot_id=uid,
                tim_barrel_accession=r["tim_barrel_accession"],
                protein_name=m.get("protein_name"),
                gene_name=gene_name,
                organism=r["organism"] or "Homo sapiens",
                reviewed=bool(m["reviewed"]) if m.get("reviewed") is not None else None,
                protein_existence=m.get("protein_existence"),
                annotation_score=m.get("annotation_score"),
            ))
        except Exception as e:
            logger.warning("Skipping protein %s: %s", uid, e)

    with get_connection(new_path) as conn:
        upsert_proteins(conn, proteins)
    logger.info("Migrated %d proteins", len(proteins))


def _migrate_isoforms(old: sqlite3.Connection, new_path: str) -> None:
    """Convert each old proteins row (canonical only) to an isoform row."""
    try:
        rows = old.execute("SELECT * FROM proteins").fetchall()
    except sqlite3.OperationalError:
        logger.warning("No old proteins table found — skipping isoform migration")
        return

    isoforms = []
    for r in rows:
        uid = r["uniprot_id"]
        seq = r["sequence"] if "sequence" in r.keys() else None
        seq_len = r["sequence_length"] if "sequence_length" in r.keys() else None
        if not seq:
            continue

        # TIM barrel location from the nested JSON
        tim_loc = None
        tb_feat = r["tim_barrel_features"] if "tim_barrel_features" in r.keys() else None
        if tb_feat:
            try:
                tb_data = json.loads(tb_feat)
                boundaries = tb_data.get("isoform_boundaries", {})
                # Take first non-empty boundary
                for v in boundaries.values():
                    if v:
                        tim_loc = v
                        break
            except (json.JSONDecodeError, AttributeError):
                pass

        # Ensembl gene ID from the nested JSON
        ensembl_transcript_id = None
        ens_ref = r["ensembl_references"] if "ensembl_references" in r.keys() else None
        if ens_ref:
            try:
                ens_data = json.loads(ens_ref)
                for mapping in ens_data.get("isoform_mappings", {}).values():
                    if mapping and isinstance(mapping, list):
                        ensembl_transcript_id = mapping[0].get("gene_id")
                        break
            except (json.JSONDecodeError, AttributeError):
                pass

        try:
            isoforms.append(Isoform(
                isoform_id=f"{uid}-1",
                uniprot_id=uid,
                is_canonical=True,
                sequence=seq,
                sequence_length=seq_len or len(seq),
                tim_barrel_location=tim_loc,
                ensembl_transcript_id=ensembl_transcript_id,
            ))
        except Exception as e:
            logger.warning("Skipping isoform for %s: %s", uid, e)

    with get_connection(new_path) as conn:
        upsert_isoforms(conn, isoforms)
    logger.info("Migrated %d isoforms (canonical only; re-run collection for alternatives)", len(isoforms))


def _print_counts(db_path: str) -> None:
    conn = sqlite3.connect(db_path)
    for table in ("tim_barrel_entries", "proteins", "isoforms"):
        n = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
        logger.info("  %-25s %d rows", table, n)
    conn.close()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Migrate protein database to new schema")
    parser.add_argument("--old", default="db/protein_data.db", help="Path to old database")
    parser.add_argument("--new", default="db/protein_data_v2.db", help="Path for new database")
    args = parser.parse_args()
    migrate(args.old, args.new)


if __name__ == "__main__":
    main()
