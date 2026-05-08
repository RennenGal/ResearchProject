#!/usr/bin/env python3
"""
Build (or rebuild) affected_isoforms.

Primary detection: VSP overlap — any alternative isoform with at least one
UniProt VSP feature overlapping the canonical domain is included.
Fallback for isoforms with no VSP data: sliding window (12.5% <= id < 95%).

Usage:
    python scripts/build_affected_isoforms.py
    python scripts/build_affected_isoforms.py --db db/protein_data.db
"""

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from protein_data_collector.analysis.tim_barrel_alignment import populate_tim_barrel_isoforms
from protein_data_collector.database.connection import get_connection
from protein_data_collector.database.schema import init_db
from protein_data_collector.config import get_config


def main() -> None:
    parser = argparse.ArgumentParser(description="Build affected_isoforms analysis table")
    parser.add_argument("--db", default=None, help="Database path (default: from config)")
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format="%(asctime)s %(levelname)-8s %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )

    db_path = args.db or get_config().db_path

    with get_connection(db_path) as conn:
        init_db(conn)
        inserted, skipped_no_overlap, skipped_fallback = populate_tim_barrel_isoforms(
            conn,
            isoform_table="isoforms",
            output_table="affected_isoforms",
        )

    print(f"\n{'='*55}")
    print(f"  AS-affected domain isoforms   : {inserted}")
    print(f"  Skipped — no domain overlap   : {skipped_no_overlap}")
    print(f"  Skipped — no VSP + low id     : {skipped_fallback}")
    print(f"  Total alternative isoforms    : {inserted + skipped_no_overlap + skipped_fallback}")
    print(f"{'='*55}")


if __name__ == "__main__":
    main()
