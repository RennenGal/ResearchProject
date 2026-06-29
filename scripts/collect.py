#!/usr/bin/env python3
"""
Run the protein data collection pipeline end-to-end.

Phases
------
Phase 1  — Fetch InterPro / UniProt entry list (DataCollector)
Phase 2  — Collect protein records (DataCollector)
Phase 3  — Collect isoforms from UniProt (DataCollector)
Phase 4  — Backfill protein metadata (protein_name, reviewed, annotation_score)
Phase 5  — Fetch and propagate gene names
Phase 6  — Backfill TIM barrel domain locations
Phase 7  — Deduplicate proteins by gene_name
Phase 8  — Build affected_isoforms table + backfill fragment isoforms
Phase 9  — Build canonical_analysis table
Phase 10 — Annotate TIM barrel motifs (AlphaFold + pydssp)
Phase 11 — Collect Ensembl transcripts + alignment analysis + backfill exon data
Phase 12 — Backfill isoform exon junction data (UniProt isoforms)
Phase 13 — Build analysis_proteins table and views
Phase 14 — (Reserved for future use)
Phase 15 — (Reserved for future use)

Usage
-----
Full collection from scratch:
    python scripts/collect.py

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
from protein_data_collector.config import get_config

from backfill_protein_metadata import run as run_backfill_protein_metadata
from fetch_gene_names import run as run_fetch_gene_names
from backfill_domain_locations import run as run_backfill_domain_locations
from dedup_by_gene import run as run_dedup_by_gene
from build_affected_isoforms import run as run_build_affected_isoforms
from build_canonical_analysis import run as run_build_canonical_analysis
from annotate_motifs import run as run_annotate_motifs
from collect_ensembl import run as run_collect_ensembl
from backfill_isoform_exons import run as run_backfill_isoform_exons
from create_analysis_table import run as run_create_analysis_table


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
    parser = argparse.ArgumentParser(description="Collect TIM barrel protein data (Homo sapiens)")
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

    collector = DataCollector(db_path=db_path, domain="tim_barrel", organism="homo_sapiens")

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

        # Phases 1–3: UniProt data collection
        report = collector.run_full_collection()
        print("\n" + report.summary())

        # Phase 4: Backfill protein metadata
        logger.info("=== Phase 4: Backfill protein metadata ===")
        run_backfill_protein_metadata(db_path)

        # Phase 5: Fetch gene names
        logger.info("=== Phase 5: Fetch gene names ===")
        run_fetch_gene_names(db_path)

        # Phase 6: Backfill domain locations
        logger.info("=== Phase 6: Backfill domain locations ===")
        run_backfill_domain_locations(db_path)

        # Phase 7: Deduplicate by gene name
        logger.info("=== Phase 7: Deduplicate by gene name ===")
        run_dedup_by_gene(db_path)

        # Phase 8: Build affected_isoforms + backfill fragment isoforms
        logger.info("=== Phase 8: Build affected_isoforms ===")
        run_build_affected_isoforms(db_path)

        # Phase 9: Build canonical_analysis
        logger.info("=== Phase 9: Build canonical_analysis ===")
        run_build_canonical_analysis(db_path)

        # Phase 10: Annotate motifs
        logger.info("=== Phase 10: Annotate motifs ===")
        run_annotate_motifs(db_path)

        # Phase 11: Collect Ensembl + backfill exon data
        logger.info("=== Phase 11: Collect Ensembl transcripts ===")
        run_collect_ensembl(db_path)

        # Phase 12: Backfill isoform exon junctions
        logger.info("=== Phase 12: Backfill isoform exon junctions ===")
        run_backfill_isoform_exons(db_path)

        # Phase 13: Build analysis_proteins table
        logger.info("=== Phase 13: Build analysis_proteins table ===")
        run_create_analysis_table(db_path)

        logger.info("Pipeline complete.")


if __name__ == "__main__":
    main()
