#!/usr/bin/env python3
"""
Backfill exon boundary data for Ensembl transcripts and flag AS-affected
transcripts whose exon junctions fall inside the TIM barrel domain.

Phase 1 — Fetch exon boundaries
    For each row in tb_ensembl_transcripts without exon_annotations, calls the
    Ensembl REST API (/lookup/id/{ENST}?expand=1) to compute protein-space exon
    boundary positions and stores them as a JSON array in exon_annotations.

    Each element is the 1-based amino-acid position of the last residue of that
    exon (ceiling of cumulative CDS bases / 3). The final exon is omitted.

Phase 2 — Flag domain-disrupting exon junctions
    For each row in tb_ensembl_affected, checks whether any exon boundary from
    the corresponding tb_ensembl_transcripts row falls strictly inside the
    domain region [domain_start, domain_end - 1].  Writes the result to
    exon_boundary_in_domain (1 or 0).

Usage:
    python scripts/backfill_exons.py
    python scripts/backfill_exons.py --phase1-only    # fetch exon data, skip flagging
    python scripts/backfill_exons.py --phase2-only    # flag only (requires Phase 1 done)
    python scripts/backfill_exons.py --limit 50       # first N transcripts (testing)
    python scripts/backfill_exons.py --db db/protein_data.db
"""

import argparse
import json
import logging
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from protein_data_collector.api.ensembl_client import transcript_exon_boundaries
from protein_data_collector.config import get_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)

_TRANSCRIPT_TABLE = "tb_ensembl_transcripts"
_AFFECTED_TABLE   = "tb_ensembl_affected"


# ---------------------------------------------------------------------------
# Schema migration — add columns to existing DB if absent
# ---------------------------------------------------------------------------

def _ensure_columns(conn: sqlite3.Connection) -> None:
    existing_enst = {r[1] for r in conn.execute(f"PRAGMA table_info({_TRANSCRIPT_TABLE})")}
    if "exon_annotations" not in existing_enst:
        conn.execute(f"ALTER TABLE {_TRANSCRIPT_TABLE} ADD COLUMN exon_annotations TEXT")
        logger.info("Added exon_annotations column to %s", _TRANSCRIPT_TABLE)

    existing_aff = {r[1] for r in conn.execute(f"PRAGMA table_info({_AFFECTED_TABLE})")}
    if "exon_boundary_in_domain" not in existing_aff:
        conn.execute(
            f"ALTER TABLE {_AFFECTED_TABLE} "
            f"ADD COLUMN exon_boundary_in_domain INTEGER NOT NULL DEFAULT 0"
        )
        logger.info("Added exon_boundary_in_domain column to %s", _AFFECTED_TABLE)
    if "exon_boundaries_in_domain_count" not in existing_aff:
        conn.execute(
            f"ALTER TABLE {_AFFECTED_TABLE} "
            f"ADD COLUMN exon_boundaries_in_domain_count INTEGER NOT NULL DEFAULT 0"
        )
        logger.info("Added exon_boundaries_in_domain_count column to %s", _AFFECTED_TABLE)

    conn.commit()


# ---------------------------------------------------------------------------
# Phase 1: fetch and store exon boundary annotations
# ---------------------------------------------------------------------------

def backfill_exon_annotations(conn: sqlite3.Connection, limit: int | None = None) -> int:
    """
    Fetch exon boundaries from Ensembl for transcripts missing exon_annotations.
    Returns the number of rows updated.
    """
    query = (
        f"SELECT enst_id FROM {_TRANSCRIPT_TABLE} "
        f"WHERE exon_annotations IS NULL"
    )
    if limit:
        query += f" LIMIT {limit}"

    rows = conn.execute(query).fetchall()
    total = len(rows)
    logger.info("Fetching exon boundaries for %d transcripts", total)

    updated = 0
    for i, (enst_id,) in enumerate(rows, 1):
        boundaries = transcript_exon_boundaries(enst_id)
        annotation = json.dumps(boundaries)
        conn.execute(
            f"UPDATE {_TRANSCRIPT_TABLE} SET exon_annotations=? WHERE enst_id=?",
            (annotation, enst_id),
        )
        updated += 1

        if i % 100 == 0 or i == total:
            conn.commit()
            logger.info("  %d / %d transcripts processed", i, total)

    conn.commit()
    return updated


# ---------------------------------------------------------------------------
# Phase 2: flag exon boundaries that fall inside the domain
# ---------------------------------------------------------------------------

def flag_exon_boundary_in_domain(conn: sqlite3.Connection) -> tuple[int, int]:
    """
    For each row in tb_ensembl_affected, check whether any exon boundary from
    the corresponding transcript's exon_annotations falls strictly inside the
    domain region [domain_start, domain_end - 1].

    Returns (flagged, total).
    """
    rows = conn.execute(f"""
        SELECT
            aff.id,
            aff.domain_location,
            et.exon_annotations
        FROM {_AFFECTED_TABLE} aff
        JOIN {_TRANSCRIPT_TABLE} et ON et.enst_id = aff.enst_id
        WHERE et.exon_annotations IS NOT NULL
    """).fetchall()

    logger.info("Evaluating exon boundaries for %d AS-affected transcripts", len(rows))

    flagged = 0
    for row_id, domain_loc_raw, exon_ann_raw in rows:
        flag = 0
        count = 0
        try:
            domain_loc  = json.loads(domain_loc_raw)
            dom_start   = domain_loc["start"]
            dom_end     = domain_loc["end"]
            boundaries  = json.loads(exon_ann_raw)
            count = sum(1 for b in boundaries if dom_start <= b < dom_end)
            if count > 0:
                flag = 1
                flagged += 1
        except (TypeError, KeyError, json.JSONDecodeError):
            pass

        conn.execute(
            f"UPDATE {_AFFECTED_TABLE} "
            f"SET exon_boundary_in_domain=?, exon_boundaries_in_domain_count=? WHERE id=?",
            (flag, count, row_id),
        )

    conn.commit()
    logger.info("Flagged %d / %d AS-affected transcripts as exon_boundary_in_domain=1",
                flagged, len(rows))
    return flagged, len(rows)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Backfill exon boundary data and flag domain-disrupting junctions"
    )
    parser.add_argument("--db",           default=None, help="Override DB path")
    parser.add_argument("--limit",        type=int, default=None,
                        help="Process only first N transcripts in Phase 1 (testing)")
    parser.add_argument("--phase1-only",  action="store_true",
                        help="Fetch exon data only; skip domain-boundary flagging")
    parser.add_argument("--phase2-only",  action="store_true",
                        help="Flag domain boundaries only (requires Phase 1 complete)")
    parser.add_argument("--log-level",    default="INFO")
    args = parser.parse_args()

    logging.getLogger().setLevel(getattr(logging, args.log_level.upper(), logging.INFO))

    db_path = args.db or get_config().db_path
    conn = sqlite3.connect(db_path)

    _ensure_columns(conn)

    if not args.phase2_only:
        logger.info("=== Phase 1: fetching exon annotations ===")
        updated = backfill_exon_annotations(conn, limit=args.limit)
        print(f"\n  Exon annotations written: {updated}")

    if not args.phase1_only:
        logger.info("=== Phase 2: flagging exon boundaries in domain ===")
        flagged, total = flag_exon_boundary_in_domain(conn)

        # Summary stats
        total_aff      = conn.execute(f"SELECT COUNT(*) FROM {_AFFECTED_TABLE}").fetchone()[0]
        with_ann       = conn.execute(
            f"SELECT COUNT(*) FROM {_AFFECTED_TABLE} WHERE exon_boundary_in_domain IS NOT NULL"
        ).fetchone()[0]

        print(f"\n{'='*60}")
        print(f"  Exon boundary analysis — TIM barrel (Homo sapiens)")
        print(f"{'='*60}")
        print(f"  AS-affected transcripts (total)          : {total_aff}")
        print(f"  Evaluated (exon data available)          : {total}")
        print(f"  Exon boundary falls inside domain        : {flagged}")
        if total > 0:
            print(f"  Fraction with intra-domain junction      : {flagged/total:.1%}")
        print(f"{'='*60}")

    conn.close()


if __name__ == "__main__":
    main()
