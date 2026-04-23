#!/usr/bin/env python3
"""
Backfill exon junction data for tb_isoforms and flag domain-disrupting
junctions in tb_affected_isoforms.

Phase 1 — Canonical isoforms
    For each canonical isoform with an ensembl_transcript_id, calls the
    Ensembl REST API to obtain protein-space exon boundary positions and
    stores them as a JSON array in tb_isoforms.exon_annotations.

Phase 2 — Alternative isoforms
    Derives exon boundary positions for each alternative isoform by
    transforming the canonical boundaries through the splice-variant (VSP)
    coordinate map.

    Strategy: the canonical and alternative sequences are identical outside
    VSP regions.  For each unchanged segment (between/before/after VSPs), the
    script finds the segment as an exact substring of the alternative sequence
    (searching sequentially, so it works for both deletions and substitutions
    without knowing the replacement sequence).  Exon boundaries inside VSP
    regions are dropped; boundaries outside are shifted to their alternative-
    space positions.

Phase 3 — Flag domain-disrupting junctions in tb_affected_isoforms
    Copies exon_annotations from tb_isoforms into tb_affected_isoforms and
    checks whether any exon junction falls strictly inside the domain region
    [domain_start, domain_end).

Usage:
    python scripts/backfill_isoform_exons.py
    python scripts/backfill_isoform_exons.py --phase1-only
    python scripts/backfill_isoform_exons.py --phase2-only
    python scripts/backfill_isoform_exons.py --phase3-only
    python scripts/backfill_isoform_exons.py --limit 20   # test: first N canonical isoforms
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

_ISOFORM_TABLE  = "tb_isoforms"
_AFFECTED_TABLE = "tb_affected_isoforms"


# ---------------------------------------------------------------------------
# Schema migration
# ---------------------------------------------------------------------------

def _ensure_columns(conn: sqlite3.Connection) -> None:
    for tbl, col, defn in [
        (_AFFECTED_TABLE, "exon_boundary_in_domain",         "INTEGER NOT NULL DEFAULT 0"),
        (_AFFECTED_TABLE, "exon_boundaries_in_domain_count", "INTEGER NOT NULL DEFAULT 0"),
    ]:
        cols = {r[1] for r in conn.execute(f"PRAGMA table_info({tbl})")}
        if col not in cols:
            conn.execute(f"ALTER TABLE {tbl} ADD COLUMN {col} {defn}")
            logger.info("Added %s to %s", col, tbl)
    conn.commit()


# ---------------------------------------------------------------------------
# Phase 1 — canonical isoform exon boundaries from Ensembl
# ---------------------------------------------------------------------------

def backfill_canonical(conn: sqlite3.Connection, limit: int | None = None) -> int:
    query = (
        f"SELECT isoform_id, ensembl_transcript_id FROM {_ISOFORM_TABLE} "
        f"WHERE is_canonical=1 AND ensembl_transcript_id IS NOT NULL "
        f"AND exon_annotations IS NULL"
    )
    if limit:
        query += f" LIMIT {limit}"

    rows = conn.execute(query).fetchall()
    logger.info("Phase 1: fetching exon boundaries for %d canonical isoforms", len(rows))

    updated = 0
    for i, (iso_id, enst_raw) in enumerate(rows, 1):
        enst = enst_raw.split(".")[0]
        boundaries = transcript_exon_boundaries(enst)
        conn.execute(
            f"UPDATE {_ISOFORM_TABLE} SET exon_annotations=? WHERE isoform_id=?",
            (json.dumps(boundaries), iso_id),
        )
        updated += 1
        if i % 100 == 0 or i == len(rows):
            conn.commit()
            logger.info("  %d / %d canonical isoforms processed", i, len(rows))

    conn.commit()
    return updated


# ---------------------------------------------------------------------------
# Phase 2 — derive alternative isoform boundaries via sequence matching
# ---------------------------------------------------------------------------

def _transform_boundaries(
    can_seq: str,
    alt_seq: str,
    vsps: list[dict],
    can_boundaries: list[int],
) -> list[int]:
    """
    Map canonical exon boundary positions (1-based) to alternative-isoform
    positions using exact-substring matching on the unchanged flanking segments.

    Boundaries inside VSP regions are omitted (that sequence is replaced).
    """
    if not vsps or not can_boundaries:
        return can_boundaries

    vsps_sorted = sorted(vsps, key=lambda v: v["location"]["start"]["value"])

    # Build (can_start, can_end, alt_start) for each unchanged segment
    segments: list[tuple[int, int, int]] = []
    prev_can_end = 0
    search_alt   = 0

    def _locate_segment(seg_can_start: int, seg_can_end: int) -> int | None:
        """Find the position of canonical[seg_can_start-1:seg_can_end] in alt,
        searching forward from search_alt. Returns alt_start (1-based) or None."""
        seg = can_seq[seg_can_start - 1:seg_can_end]
        if not seg:
            return None
        pos = alt_seq.find(seg, search_alt)
        return (pos + 1) if pos >= 0 else None

    for vsp in vsps_sorted:
        vs = vsp["location"]["start"]["value"]  # 1-based in canonical
        ve = vsp["location"]["end"]["value"]

        if vs > prev_can_end + 1:
            seg_s = prev_can_end + 1
            seg_e = vs - 1
            alt_start = _locate_segment(seg_s, seg_e)
            if alt_start is not None:
                segments.append((seg_s, seg_e, alt_start))
                search_alt = (alt_start - 1) + (seg_e - seg_s + 1)

        prev_can_end = ve

    # Suffix after last VSP
    if prev_can_end < len(can_seq):
        seg_s = prev_can_end + 1
        seg_e = len(can_seq)
        alt_start = _locate_segment(seg_s, seg_e)
        if alt_start is not None:
            segments.append((seg_s, seg_e, alt_start))

    # Transform each canonical boundary
    alt_boundaries: list[int] = []
    for b in can_boundaries:
        for seg_s, seg_e, seg_alt_s in segments:
            if seg_s <= b <= seg_e:
                alt_boundaries.append(b + (seg_alt_s - seg_s))
                break
        # If no segment matched: b is inside a VSP region — omit

    return alt_boundaries


def backfill_alternative(conn: sqlite3.Connection) -> int:
    """
    Derive and store exon boundaries for all alternative isoforms whose
    canonical protein has exon data.  Returns the number of rows updated.
    """
    rows = conn.execute(f"""
        SELECT
            alt.isoform_id,
            alt.sequence,
            alt.splice_variants,
            can.sequence        AS can_seq,
            can.exon_annotations AS can_exon_ann
        FROM {_ISOFORM_TABLE} alt
        JOIN {_ISOFORM_TABLE} can
          ON  can.uniprot_id   = alt.uniprot_id
          AND can.is_canonical = 1
        WHERE alt.is_canonical  = 0
          AND alt.exon_annotations IS NULL
          AND can.exon_annotations IS NOT NULL
    """).fetchall()

    logger.info("Phase 2: deriving exon boundaries for %d alternative isoforms", len(rows))

    updated = 0
    for iso_id, alt_seq, sv_raw, can_seq, can_ann_raw in rows:
        can_boundaries = json.loads(can_ann_raw) if can_ann_raw else []
        vsps = json.loads(sv_raw) if sv_raw else []

        if not can_boundaries:
            conn.execute(
                f"UPDATE {_ISOFORM_TABLE} SET exon_annotations='[]' WHERE isoform_id=?",
                (iso_id,),
            )
            updated += 1
            continue

        alt_boundaries = _transform_boundaries(can_seq, alt_seq, vsps, can_boundaries)
        conn.execute(
            f"UPDATE {_ISOFORM_TABLE} SET exon_annotations=? WHERE isoform_id=?",
            (json.dumps(alt_boundaries), iso_id),
        )
        updated += 1

    conn.commit()
    return updated


# ---------------------------------------------------------------------------
# Phase 3 — copy exon data into tb_affected_isoforms and flag domain junctions
# ---------------------------------------------------------------------------

def flag_domain_boundaries(conn: sqlite3.Connection) -> tuple[int, int]:
    """
    1. Copy exon_annotations from tb_isoforms into tb_affected_isoforms.
    2. For each row, check whether any boundary falls in [domain_start, domain_end).
    Returns (flagged, total_evaluated).
    """
    # Copy exon annotations
    conn.execute(f"""
        UPDATE {_AFFECTED_TABLE}
        SET exon_annotations = (
            SELECT exon_annotations FROM {_ISOFORM_TABLE}
            WHERE isoform_id = {_AFFECTED_TABLE}.isoform_id
        )
        WHERE EXISTS (
            SELECT 1 FROM {_ISOFORM_TABLE}
            WHERE isoform_id = {_AFFECTED_TABLE}.isoform_id
              AND exon_annotations IS NOT NULL
        )
    """)
    conn.commit()

    # Compute flags
    rows = conn.execute(f"""
        SELECT isoform_id, domain_location, exon_annotations
        FROM {_AFFECTED_TABLE}
        WHERE exon_annotations IS NOT NULL
    """).fetchall()

    logger.info("Phase 3: flagging domain boundaries for %d affected isoforms", len(rows))
    flagged = 0
    for iso_id, dom_raw, ann_raw in rows:
        flag = count = 0
        try:
            dom   = json.loads(dom_raw)
            dom_s = dom["start"]
            dom_e = dom["end"]
            bounds = json.loads(ann_raw)
            count = sum(1 for b in bounds if dom_s <= b < dom_e)
            flag  = 1 if count > 0 else 0
            if flag:
                flagged += 1
        except (TypeError, KeyError, json.JSONDecodeError):
            pass

        conn.execute(
            f"UPDATE {_AFFECTED_TABLE} "
            f"SET exon_boundary_in_domain=?, exon_boundaries_in_domain_count=? "
            f"WHERE isoform_id=?",
            (flag, count, iso_id),
        )

    conn.commit()
    return flagged, len(rows)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Backfill exon junction data for UniProt isoforms"
    )
    parser.add_argument("--db",           default=None)
    parser.add_argument("--limit",        type=int, default=None,
                        help="Process only first N canonical isoforms in Phase 1 (testing)")
    parser.add_argument("--phase1-only",  action="store_true")
    parser.add_argument("--phase2-only",  action="store_true")
    parser.add_argument("--phase3-only",  action="store_true")
    parser.add_argument("--log-level",    default="INFO")
    args = parser.parse_args()

    logging.getLogger().setLevel(getattr(logging, args.log_level.upper(), logging.INFO))
    db_path = args.db or get_config().db_path
    conn = sqlite3.connect(db_path)

    _ensure_columns(conn)

    run_all = not (args.phase1_only or args.phase2_only or args.phase3_only)

    if run_all or args.phase1_only:
        n = backfill_canonical(conn, limit=args.limit)
        print(f"\n  Phase 1 complete: {n} canonical isoforms updated")

    if run_all or args.phase2_only:
        n = backfill_alternative(conn)
        print(f"  Phase 2 complete: {n} alternative isoforms updated")

    if run_all or args.phase3_only:
        flagged, total = flag_domain_boundaries(conn)

        total_aff  = conn.execute(f"SELECT COUNT(*) FROM {_AFFECTED_TABLE}").fetchone()[0]
        evaluated  = conn.execute(
            f"SELECT COUNT(*) FROM {_AFFECTED_TABLE} WHERE exon_annotations IS NOT NULL"
        ).fetchone()[0]

        print(f"\n{'='*60}")
        print(f"  Exon junction analysis — UniProt isoforms")
        print(f"{'='*60}")
        print(f"  AS-affected isoforms (total)        : {total_aff}")
        print(f"  With exon data (coverage)           : {evaluated} / {total_aff}")
        print(f"  Exon junction inside domain         : {flagged}")
        if evaluated:
            print(f"  Fraction with intra-domain junction : {flagged/evaluated:.1%}")
        print(f"{'='*60}")

    conn.close()


if __name__ == "__main__":
    main()
