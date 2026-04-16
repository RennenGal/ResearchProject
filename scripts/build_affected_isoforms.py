#!/usr/bin/env python3
"""
Build (or rebuild) the {domain}_affected_isoforms analysis table.

For each alternative, non-fragment isoform, performs an ungapped local
alignment of the canonical domain sequence against the isoform and
inserts rows where 12.5% <= identity < 95%.

Usage:
    python scripts/build_affected_isoforms.py
    python scripts/build_affected_isoforms.py --domain beta_propeller --organism mus_musculus
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
from protein_data_collector.config import get_config, DOMAINS, ORGANISMS


def main() -> None:
    parser = argparse.ArgumentParser(description="Build affected_isoforms analysis table")
    parser.add_argument("--domain", default="tim_barrel",
                        choices=list(DOMAINS),
                        help="Domain to analyse (default: tim_barrel)")
    parser.add_argument("--organism", default="homo_sapiens",
                        choices=list(ORGANISMS),
                        help="Organism to analyse (default: homo_sapiens)")
    parser.add_argument("--db", default=None, help="Database path (default: from config)")
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format="%(asctime)s %(levelname)-8s %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )

    domain_cfg   = DOMAINS[args.domain]
    organism_cfg = ORGANISMS[args.organism]
    db_path      = args.db or get_config().db_path

    isoform_table  = organism_cfg.isoform_table(domain_cfg)
    output_table   = organism_cfg.affected_isoforms_table(domain_cfg)

    with get_connection(db_path) as conn:
        init_db(conn)
        inserted, skipped_identical, skipped_absent, insertions = populate_tim_barrel_isoforms(
            conn,
            isoform_table=isoform_table,
            output_table=output_table,
        )

    print(f"\n{'='*55}")
    print(f"  Domain   : {domain_cfg.display_name}")
    print(f"  Organism : {organism_cfg.display_name}")
    print(f"{'='*55}")
    print(f"  AS-affected domain isoforms   : {inserted}")
    print(f"    of which insertion-detected : {insertions}")
    print(f"  Skipped — domain >= 95%       : {skipped_identical}")
    print(f"  Skipped — domain absent       : {skipped_absent}")
    print(f"  Total alternative isoforms    : {inserted + skipped_identical + skipped_absent}")
    print(f"{'='*55}")


if __name__ == "__main__":
    main()
