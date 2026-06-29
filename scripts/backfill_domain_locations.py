#!/usr/bin/env python3
"""
Re-query InterPro for ALL canonical proteins and store the complete list of
domain instances in isoforms.{location_col} (replacing the existing
single-dict format with a list).

This is required before rebuilding canonical_analysis with the new
(uniprot_id, domain_index) primary key so that proteins with two domain
instances are represented by two rows rather than one.

Usage
-----
    python scripts/backfill_domain_locations.py
    python scripts/backfill_domain_locations.py --db db/protein_data.db
    python scripts/backfill_domain_locations.py --dry-run   # print only
"""

import json
import logging
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from protein_data_collector.api.interpro_client import InterProClient
from protein_data_collector.config import DOMAINS

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


def _run_backfill(conn: sqlite3.Connection, accession_col: str, location_col: str,
                  dry_run: bool = False) -> dict:
    rows = conn.execute(f"""
        SELECT iso.isoform_id, iso.uniprot_id, p.{accession_col},
               iso.{location_col}
        FROM   isoforms iso
        JOIN   proteins p ON p.uniprot_id = iso.uniprot_id
        WHERE  iso.is_canonical          = 1
          AND  iso.is_fragment           = 0
          AND  p.canonical_uniprot_id   IS NULL
          AND  p.{accession_col}        IS NOT NULL
        ORDER BY iso.uniprot_id
    """).fetchall()

    total = len(rows)
    logger.info("Canonical non-fragment proteins to check: %d", total)

    client = InterProClient()
    stats = {"updated": 0, "no_locations": 0, "multi_domain": 0, "errors": 0}

    for i, row in enumerate(rows, 1):
        isoform_id  = row[0]
        uniprot_id  = row[1]
        accession   = row[2]
        existing_raw = row[3]
        try:
            locations = client.get_domain_boundaries(uniprot_id, accession)
        except Exception as e:
            logger.error("API error for %s: %s", uniprot_id, e)
            stats["errors"] += 1
            continue

        if not locations:
            logger.warning("No domain locations returned for %s (%s)", uniprot_id, accession)
            stats["no_locations"] += 1
            continue

        if len(locations) > 1:
            stats["multi_domain"] += 1
            logger.info(
                "Multi-domain: %s has %d domain instances: %s",
                uniprot_id,
                len(locations),
                [(d["start"], d["end"]) for d in locations],
            )

        new_json = json.dumps(locations)

        if not dry_run:
            conn.execute(
                f"UPDATE isoforms SET {location_col} = ? WHERE isoform_id = ?",
                (new_json, isoform_id),
            )

        stats["updated"] += 1

        if i % 50 == 0:
            if not dry_run:
                conn.commit()
            logger.info("  %d / %d processed", i, total)

    if not dry_run:
        conn.commit()

    return stats


def run(db_path: str, domain: str = "tim_barrel") -> None:
    domain_cfg = DOMAINS[domain]
    accession_col = domain_cfg.accession_col
    location_col  = domain_cfg.location_col

    conn = sqlite3.connect(db_path)

    stats = _run_backfill(conn, accession_col=accession_col, location_col=location_col,
                          dry_run=False)

    print(f"\n{'='*55}")
    print("  Domain location backfill summary")
    print(f"{'='*55}")
    print(f"  Proteins processed          : {stats['updated']}")
    print(f"  Multi-domain proteins found : {stats['multi_domain']}")
    print(f"  No locations returned       : {stats['no_locations']}")
    print(f"  API errors                  : {stats['errors']}")
    print(f"{'='*55}")

    conn.close()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Backfill domain locations")
    parser.add_argument("--db", default="db/protein_data.db")
    parser.add_argument("--domain", default="tim_barrel")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    if args.dry_run:
        domain_cfg = DOMAINS[args.domain]
        conn = sqlite3.connect(args.db)
        stats = _run_backfill(conn, accession_col=domain_cfg.accession_col,
                              location_col=domain_cfg.location_col, dry_run=True)
        conn.close()
        print(stats)
    else:
        run(args.db, domain=args.domain)
