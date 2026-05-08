#!/usr/bin/env python3
"""
Backfill protein_name, reviewed, and annotation_score for all proteins that
currently have NULL values for these fields, then re-run deduplication so that
unreviewed / low-quality entries are correctly marked as redundant.

Steps
-----
1. Fetch metadata from UniProt in batches for every protein with any NULL field.
2. UPDATE proteins with the fetched values.
3. Run deduplicate_proteins() so entries sharing the same protein_name are
   consolidated and redundant ones get canonical_uniprot_id set.

Usage
-----
    python scripts/backfill_protein_metadata.py
    python scripts/backfill_protein_metadata.py --db db/protein_data.db
    python scripts/backfill_protein_metadata.py --no-dedup   # skip deduplication
    python scripts/backfill_protein_metadata.py --all        # re-fetch even filled rows
"""

import argparse
import logging
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from protein_data_collector.api.uniprot_client import UniProtClient
from protein_data_collector.config import get_config
from protein_data_collector.database.storage import deduplicate_proteins

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


def fetch_and_update(conn: sqlite3.Connection, refetch_all: bool = False) -> int:
    if refetch_all:
        rows = conn.execute("SELECT uniprot_id FROM proteins").fetchall()
    else:
        rows = conn.execute(
            "SELECT uniprot_id FROM proteins "
            "WHERE protein_name IS NULL OR reviewed IS NULL OR annotation_score IS NULL"
        ).fetchall()

    if not rows:
        logger.info("All proteins already have metadata — skipping fetch")
        return 0

    ids = [r[0] for r in rows]
    logger.info("Fetching metadata for %d proteins from UniProt ...", len(ids))

    client = UniProtClient()
    meta = client.batch_protein_metadata(ids)

    filled = sum(1 for v in meta.values() if v.get("protein_name"))
    logger.info("protein_name resolved for %d / %d", filled, len(ids))

    conn.executemany(
        "UPDATE proteins SET protein_name=?, reviewed=?, annotation_score=? WHERE uniprot_id=?",
        [
            (v["protein_name"], v["reviewed"], v["annotation_score"], uid)
            for uid, v in meta.items()
        ],
    )
    conn.commit()
    return filled


def main() -> None:
    parser = argparse.ArgumentParser(description="Backfill protein metadata and re-deduplicate")
    parser.add_argument("--db",       default=None)
    parser.add_argument("--no-dedup", action="store_true", help="Skip deduplication step")
    parser.add_argument("--all",      action="store_true", help="Re-fetch even already-filled rows")
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()

    logging.getLogger().setLevel(getattr(logging, args.log_level.upper(), logging.INFO))

    db_path = args.db or get_config().db_path
    conn = sqlite3.connect(db_path)

    logger.info("=== Step 1: fetch protein metadata from UniProt ===")
    filled = fetch_and_update(conn, refetch_all=args.all)

    total  = conn.execute("SELECT COUNT(*) FROM proteins").fetchone()[0]
    filled_now = conn.execute(
        "SELECT COUNT(*) FROM proteins WHERE protein_name IS NOT NULL"
    ).fetchone()[0]
    reviewed_now = conn.execute(
        "SELECT COUNT(*) FROM proteins WHERE reviewed=1"
    ).fetchone()[0]
    print(f"\n  protein_name filled : {filled_now} / {total}")
    print(f"  reviewed (Swiss-Prot): {reviewed_now} / {total}")

    if not args.no_dedup:
        logger.info("=== Step 2: deduplicate proteins ===")
        before = conn.execute(
            "SELECT COUNT(*) FROM proteins WHERE canonical_uniprot_id IS NULL"
        ).fetchone()[0]
        deduplicate_proteins(conn)
        after = conn.execute(
            "SELECT COUNT(*) FROM proteins WHERE canonical_uniprot_id IS NULL"
        ).fetchone()[0]
        redundant = conn.execute(
            "SELECT COUNT(*) FROM proteins WHERE canonical_uniprot_id IS NOT NULL"
        ).fetchone()[0]
        print(f"\n  Canonical proteins before dedup : {before}")
        print(f"  Canonical proteins after dedup  : {after}")
        print(f"  Marked redundant                : {redundant}")

    conn.close()
    print("\nDone.")


if __name__ == "__main__":
    main()
