#!/usr/bin/env python3
"""
Run the protein data collection pipeline.

Usage
-----
Full collection from scratch (TIM barrel, human):
    python scripts/collect.py

Collect beta propeller data for mouse:
    python scripts/collect.py --domain beta_propeller --organism mus_musculus

Resume isoform collection for proteins not yet processed:
    python scripts/collect.py --resume

Custom database path:
    python scripts/collect.py --db db/my_db.db

Log to file:
    python scripts/collect.py --log-file logs/collect.log
"""

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from protein_data_collector.collector.data_collector import DataCollector
from protein_data_collector.config import get_config, DOMAINS, ORGANISMS


def setup_logging(log_file: str = None, level: str = "INFO") -> None:
    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file))
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s %(levelname)-8s %(name)s — %(message)s",
        handlers=handlers,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Collect domain protein data")
    parser.add_argument("--domain", default="tim_barrel",
                        choices=list(DOMAINS),
                        help="Domain to collect data for (default: tim_barrel)")
    parser.add_argument("--organism", default="homo_sapiens",
                        choices=list(ORGANISMS),
                        help="Organism to collect data for (default: homo_sapiens)")
    parser.add_argument("--resume", action="store_true",
                        help="Only collect isoforms for proteins not yet processed")
    parser.add_argument("--collect-proteins", action="store_true",
                        help="Phase 1+2 only: update entries and proteins; skip isoform collection. "
                             "Run --resume afterwards to collect isoforms for new proteins.")
    parser.add_argument("--recollect-isoforms", action="store_true",
                        help="Delete all isoforms and re-fetch from UniProt")
    parser.add_argument("--backfill-domains", action="store_true",
                        help="Fetch domain location for canonical isoforms where it is NULL")
    parser.add_argument("--db", default=None, help="Override database path from config")
    parser.add_argument("--log-file", default=None, help="Also write logs to this file")
    parser.add_argument("--log-level", default="INFO", help="Logging level (default: INFO)")
    args = parser.parse_args()

    setup_logging(args.log_file, args.log_level)
    logger = logging.getLogger(__name__)

    db_path = args.db or get_config().db_path
    logger.info("Database : %s", db_path)
    logger.info("Domain   : %s", args.domain)
    logger.info("Organism : %s", args.organism)

    collector = DataCollector(db_path=db_path, domain=args.domain, organism=args.organism)

    if args.collect_proteins:
        logger.info("Collecting entries and proteins (Phase 1+2 only)...")
        report = collector.collect_entries_and_proteins()
        print("\n" + report.summary())
    elif args.backfill_domains:
        logger.info("Backfilling domain locations for canonical isoforms...")
        updated = collector.backfill_domain_locations()
        print(f"\nUpdated {updated} isoforms with domain location.")
    elif args.recollect_isoforms:
        logger.info("Re-collecting all isoforms from UniProt...")
        report = collector.recollect_all_isoforms()
        print("\n" + report.summary())
    elif args.resume:
        logger.info("Resuming isoform collection...")
        report = collector.resume_isoform_collection()
        print("\n" + report.summary())
    else:
        logger.info("Starting full collection pipeline...")
        report = collector.run_full_collection()
        print("\n" + report.summary())


if __name__ == "__main__":
    main()
