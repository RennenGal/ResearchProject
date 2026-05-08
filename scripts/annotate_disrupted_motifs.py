#!/usr/bin/env python3
"""
Annotate disrupted_motifs for AS-affected isoforms.

For affected_isoforms (UniProt VSP data)
---------------------------------------------
Each row has vsp_domain_events with the canonical overlap coordinates of
every splice event.  A motif N is disrupted by a VSP if the VSP's
[overlap_start, overlap_end] intersects with [motif.beta_start, motif.alpha_end].

Output per row (JSON list, one entry per VSP):
  [{"feature_id": ..., "overlap_start": ..., "overlap_end": ...,
    "disrupted_motifs": [1, 4, 5]}, ...]

For ensembl_affected (sliding-window data, no VSP)
------------------------------------------------------
canonical_domain_sequence and domain_sequence are the same length (ungapped
sliding-window alignment).  Consecutive positions where they differ are grouped
into "changed regions" in canonical coordinates, which are then mapped to motifs.

Output per row (JSON list, one entry per changed region):
  [{"can_start": ..., "can_end": ..., "disrupted_motifs": [3]}, ...]

Usage
-----
    python scripts/annotate_disrupted_motifs.py
    python scripts/annotate_disrupted_motifs.py --db db/protein_data.db
"""

import argparse
import json
import logging
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from protein_data_collector.config import get_config

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Core: which motifs does an interval [start, end] overlap?
# ---------------------------------------------------------------------------

def overlapping_motifs(start: int, end: int, motifs: list[dict]) -> list[int]:
    """
    Return sorted list of motif numbers whose full span [beta_start, alpha_end]
    overlaps with [start, end].
    """
    hit = []
    for m in motifs:
        if m["beta_start"] <= end and m["alpha_end"] >= start:
            hit.append(m["motif"])
    return hit


# ---------------------------------------------------------------------------
# UniProt: VSP-based annotation
# ---------------------------------------------------------------------------

def annotate_uniprot(conn: sqlite3.Connection) -> tuple[int, int]:
    """
    Populate disrupted_motifs in affected_isoforms.
    Returns (updated, skipped_no_motifs).
    """
    _ensure_column(conn, "affected_isoforms", "disrupted_motifs", "TEXT")

    rows = conn.execute("""
        SELECT ai.isoform_id, ai.vsp_domain_events,
               ca.motif_annotations
        FROM affected_isoforms ai
        JOIN canonical_analysis ca ON ca.uniprot_id = ai.uniprot_id
        WHERE ai.vsp_domain_events IS NOT NULL
          AND ca.motif_annotations IS NOT NULL
    """).fetchall()

    updated = 0
    skipped = 0

    for iso_id, vsp_json, motif_json in rows:
        vsps   = json.loads(vsp_json)
        motifs = json.loads(motif_json)
        if not motifs:
            skipped += 1
            continue

        events = []
        for v in vsps:
            hit = overlapping_motifs(v["overlap_start"], v["overlap_end"], motifs)
            events.append({
                "feature_id":      v["feature_id"],
                "overlap_start":   v["overlap_start"],
                "overlap_end":     v["overlap_end"],
                "disrupted_motifs": hit,
            })

        conn.execute(
            "UPDATE affected_isoforms SET disrupted_motifs = ? WHERE isoform_id = ?",
            (json.dumps(events), iso_id),
        )
        updated += 1

    conn.commit()
    return updated, skipped


# ---------------------------------------------------------------------------
# Ensembl: sequence-comparison annotation
# ---------------------------------------------------------------------------

def _changed_regions(can_seq: str, alt_seq: str, domain_start: int) -> list[dict]:
    """
    Walk two ungapped sequences of equal length and group consecutive differing
    positions into changed regions in canonical coordinates.
    """
    if len(can_seq) != len(alt_seq):
        return []

    regions = []
    in_change = False
    reg_start = 0

    for i, (c, a) in enumerate(zip(can_seq, alt_seq)):
        canon_pos = domain_start + i
        if c != a:
            if not in_change:
                in_change  = True
                reg_start  = canon_pos
            reg_end = canon_pos
        else:
            if in_change:
                regions.append({"can_start": reg_start, "can_end": reg_end})
                in_change = False

    if in_change:
        regions.append({"can_start": reg_start, "can_end": reg_end})

    return regions


def annotate_ensembl(conn: sqlite3.Connection) -> tuple[int, int]:
    """
    Populate disrupted_motifs in ensembl_affected.
    Returns (updated, skipped_no_motifs).
    """
    _ensure_column(conn, "ensembl_affected", "disrupted_motifs", "TEXT")

    rows = conn.execute("""
        SELECT ea.id, ea.canonical_domain_sequence, ea.domain_sequence,
               ea.canonical_domain_location, ca.motif_annotations
        FROM ensembl_affected ea
        JOIN canonical_analysis ca ON ca.uniprot_id = ea.uniprot_id
        WHERE ea.domain_sequence IS NOT NULL
          AND ea.canonical_domain_sequence IS NOT NULL
          AND ca.motif_annotations IS NOT NULL
    """).fetchall()

    updated = 0
    skipped = 0

    for row_id, can_seq, alt_seq, can_loc_json, motif_json in rows:
        motifs = json.loads(motif_json)
        if not motifs:
            skipped += 1
            continue

        can_loc     = json.loads(can_loc_json)
        domain_start = can_loc["start"]

        regions = _changed_regions(can_seq, alt_seq, domain_start)
        events  = []
        for reg in regions:
            hit = overlapping_motifs(reg["can_start"], reg["can_end"], motifs)
            events.append({**reg, "disrupted_motifs": hit})

        conn.execute(
            "UPDATE ensembl_affected SET disrupted_motifs = ? WHERE id = ?",
            (json.dumps(events), row_id),
        )
        updated += 1

    conn.commit()
    return updated, skipped


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _ensure_column(conn: sqlite3.Connection, table: str, col: str, typedef: str) -> None:
    existing = {r[1] for r in conn.execute(f"PRAGMA table_info({table})")}
    if col not in existing:
        conn.execute(f"ALTER TABLE {table} ADD COLUMN {col} {typedef}")
        logger.info("Added column %s to %s", col, table)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Annotate disrupted_motifs for AS-affected isoforms"
    )
    parser.add_argument("--db", default=None)
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format="%(asctime)s %(levelname)-8s %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )

    db_path = args.db or get_config().db_path
    conn    = sqlite3.connect(db_path)

    logger.info("Annotating UniProt affected isoforms (VSP-based)...")
    u_updated, u_skipped = annotate_uniprot(conn)
    logger.info("  Updated: %d  |  Skipped (no motifs): %d", u_updated, u_skipped)

    logger.info("Annotating Ensembl affected transcripts (sequence comparison)...")
    e_updated, e_skipped = annotate_ensembl(conn)
    logger.info("  Updated: %d  |  Skipped (no motifs): %d", e_updated, e_skipped)

    conn.close()

    print(f"\n{'='*50}")
    print(f"  disrupted_motifs annotation complete")
    print(f"{'='*50}")
    print(f"  UniProt isoforms annotated  : {u_updated}")
    print(f"  Ensembl transcripts annotated: {e_updated}")
    print(f"  Skipped (no motif data)     : {u_skipped + e_skipped}")
    print(f"{'='*50}")


if __name__ == "__main__":
    main()
