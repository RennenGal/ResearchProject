#!/usr/bin/env python3
"""
Fetch primary gene names from UniProt and populate them across all relevant
TIM barrel tables (human only).

Steps
-----
1. Add gene_name column to tables that don't have it yet (idempotent ALTER TABLE).
2. Fetch gene names from the UniProt search API in batches for all proteins
   in tb_proteins that still have gene_name IS NULL.
3. UPDATE tb_proteins with the fetched names.
4. Propagate tb_proteins.gene_name into all downstream tables via JOIN UPDATE.

Usage
-----
    python scripts/fetch_gene_names.py
    python scripts/fetch_gene_names.py --db db/protein_data.db
    python scripts/fetch_gene_names.py --propagate-only   # skip fetch, just propagate
"""

import argparse
import logging
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from protein_data_collector.api.uniprot_client import UniProtClient
from protein_data_collector.config import get_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)

# Tables that need a gene_name column added (idempotent)
_COLUMNS_TO_ADD = [
    ("tb_isoforms",          "gene_name", "TEXT"),
    ("tb_affected_isoforms", "gene_name", "TEXT"),
    ("tb_ensembl_affected",  "gene_name", "TEXT"),
    # tb_proteins and tb_ensembl_transcripts already have gene_name in schema
]

# Tables to propagate gene_name into via JOIN with tb_proteins
_PROPAGATE_TARGETS = [
    # (target_table, join_column)
    ("tb_isoforms",            "uniprot_id"),
    ("tb_affected_isoforms",   "uniprot_id"),
    ("tb_ensembl_transcripts", "uniprot_id"),
    ("tb_ensembl_affected",    "uniprot_id"),
]


# ---------------------------------------------------------------------------
# Step 1 — ensure gene_name column exists everywhere
# ---------------------------------------------------------------------------

def ensure_gene_name_columns(conn: sqlite3.Connection) -> None:
    for table, col, defn in _COLUMNS_TO_ADD:
        existing = {r[1] for r in conn.execute(f"PRAGMA table_info({table})")}
        if col not in existing:
            conn.execute(f"ALTER TABLE {table} ADD COLUMN {col} {defn}")
            logger.info("Added %s to %s", col, table)
    conn.commit()


# ---------------------------------------------------------------------------
# Step 2 — fetch gene names from UniProt
# ---------------------------------------------------------------------------

def fetch_gene_names(conn: sqlite3.Connection) -> int:
    rows = conn.execute(
        "SELECT uniprot_id FROM tb_proteins WHERE gene_name IS NULL"
    ).fetchall()

    if not rows:
        logger.info("All tb_proteins rows already have gene_name — skipping fetch")
        return 0

    ids = [r[0] for r in rows]
    logger.info("Fetching gene names for %d proteins from UniProt", len(ids))

    client = UniProtClient()
    name_map = client.batch_gene_names(ids)

    found = sum(1 for v in name_map.values() if v)
    logger.info("Gene names found: %d / %d", found, len(ids))

    conn.executemany(
        "UPDATE tb_proteins SET gene_name=? WHERE uniprot_id=?",
        [(name, uid) for uid, name in name_map.items()],
    )
    conn.commit()
    return found


# ---------------------------------------------------------------------------
# Step 3 — propagate from tb_proteins to downstream tables
# ---------------------------------------------------------------------------

def propagate_gene_names(conn: sqlite3.Connection) -> None:
    for table, join_col in _PROPAGATE_TARGETS:
        # Check the table exists
        exists = conn.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", (table,)
        ).fetchone()
        if not exists:
            logger.warning("Table %s not found — skipping propagation", table)
            continue

        existing_cols = {r[1] for r in conn.execute(f"PRAGMA table_info({table})")}
        if "gene_name" not in existing_cols:
            logger.warning("gene_name column missing from %s — skipping", table)
            continue

        conn.execute(f"""
            UPDATE {table}
            SET gene_name = (
                SELECT gene_name FROM tb_proteins
                WHERE uniprot_id = {table}.{join_col}
            )
            WHERE gene_name IS NULL
        """)
        updated = conn.execute(
            f"SELECT changes()"
        ).fetchone()[0]
        conn.commit()
        logger.info("Propagated gene_name to %s: %d rows updated", table, updated)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch and propagate gene names")
    parser.add_argument("--db",              default=None)
    parser.add_argument("--propagate-only",  action="store_true",
                        help="Skip UniProt fetch; only propagate existing tb_proteins.gene_name")
    parser.add_argument("--log-level",       default="INFO")
    args = parser.parse_args()

    logging.getLogger().setLevel(getattr(logging, args.log_level.upper(), logging.INFO))

    db_path = args.db or get_config().db_path
    conn = sqlite3.connect(db_path)

    logger.info("=== Step 1: ensure gene_name columns ===")
    ensure_gene_name_columns(conn)

    if not args.propagate_only:
        logger.info("=== Step 2: fetch gene names from UniProt ===")
        found = fetch_gene_names(conn)

        filled = conn.execute(
            "SELECT COUNT(*) FROM tb_proteins WHERE gene_name IS NOT NULL"
        ).fetchone()[0]
        total = conn.execute("SELECT COUNT(*) FROM tb_proteins").fetchone()[0]
        print(f"\n  Gene names in tb_proteins: {filled} / {total}")

    logger.info("=== Step 3: propagate to downstream tables ===")
    propagate_gene_names(conn)

    conn.close()
    print("\nDone.")


if __name__ == "__main__":
    main()
