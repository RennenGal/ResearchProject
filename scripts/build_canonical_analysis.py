#!/usr/bin/env python3
"""
Build the canonical_analysis table from isoforms + proteins.

One row per canonical, non-fragment human TIM barrel protein.

Fields populated here
---------------------
  uniprot_id       — from isoforms
  gene_name        — from proteins (populated by fetch_gene_names.py first)
  sequence         — full protein sequence (from isoforms)
  domain_start     — from tim_barrel_location JSON
  domain_end       — from tim_barrel_location JSON
  domain_sequence  — from tim_barrel_sequence (or sliced if missing)
  exon_annotations — reformatted from flat boundary list to
                     [{exon, start, end}] (1-based, inclusive)
  motif_annotations — NULL (filled later by annotate_motifs.py)

Run after fetch_gene_names.py to get gene_name populated.

Usage
-----
    python scripts/build_canonical_analysis.py
    python scripts/build_canonical_analysis.py --db db/protein_data.db
    python scripts/build_canonical_analysis.py --rebuild   # DROP and recreate
"""

import argparse
import json
import logging
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from protein_data_collector.config import get_config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)

_TABLE = "canonical_analysis"

_CREATE_TABLE = f"""
CREATE TABLE IF NOT EXISTS {_TABLE} (
    uniprot_id          TEXT PRIMARY KEY,
    gene_name           TEXT,
    sequence            TEXT NOT NULL,
    domain_start        INTEGER,
    domain_end          INTEGER,
    domain_sequence     TEXT,
    exon_annotations    TEXT,
    motif_annotations   TEXT,
    dssp_source         TEXT,
    hmmer_source        TEXT,
    created_at          DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (uniprot_id) REFERENCES proteins(uniprot_id) ON DELETE CASCADE
);
CREATE INDEX IF NOT EXISTS idx_can_gene ON {_TABLE}(gene_name);
"""


# ---------------------------------------------------------------------------
# Exon annotation reformatting
# ---------------------------------------------------------------------------

def _reformat_exon_annotations(boundaries: list[int], seq_len: int) -> list[dict]:
    """
    Convert flat boundary list (1-based end of each exon, final omitted) to
    [{exon, start, end}] with 1-based inclusive coordinates.

    Example: boundaries=[85, 116, 209], seq_len=300
      → [{exon:1, start:1, end:85},
         {exon:2, start:86, end:116},
         {exon:3, start:117, end:209},
         {exon:4, start:210, end:300}]
    """
    exons: list[dict] = []
    prev_end = 0
    for i, end_pos in enumerate(sorted(boundaries), 1):
        exons.append({"exon": i, "start": prev_end + 1, "end": end_pos})
        prev_end = end_pos
    exons.append({"exon": len(boundaries) + 1, "start": prev_end + 1, "end": seq_len})
    return exons


# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------

def ensure_table(conn: sqlite3.Connection) -> None:
    conn.executescript(_CREATE_TABLE)
    conn.commit()


def build(conn: sqlite3.Connection, rebuild: bool = False) -> int:
    if rebuild:
        conn.execute(f"DELETE FROM {_TABLE}")
        conn.commit()
        logger.info("Cleared existing rows from %s", _TABLE)

    rows = conn.execute("""
        SELECT
            iso.uniprot_id,
            p.gene_name,
            iso.sequence,
            iso.sequence_length,
            iso.tim_barrel_location,
            iso.tim_barrel_sequence,
            iso.exon_annotations
        FROM isoforms iso
        JOIN proteins p ON p.uniprot_id = iso.uniprot_id
        WHERE iso.is_canonical = 1
          AND iso.is_fragment  = 0
          AND p.canonical_uniprot_id IS NULL
        ORDER BY iso.uniprot_id
    """).fetchall()

    logger.info("Processing %d canonical non-fragment proteins", len(rows))

    inserted = 0
    skipped  = 0

    for (uid, gene_name, seq, seq_len,
         tb_loc_raw, tb_seq, exon_ann_raw) in rows:

        # Parse domain location
        domain_start = domain_end = None
        domain_sequence = tb_seq
        if tb_loc_raw:
            try:
                loc = json.loads(tb_loc_raw)
                domain_start = loc.get("start")
                domain_end   = loc.get("end")
                if domain_sequence is None and domain_start and domain_end and seq:
                    domain_sequence = seq[domain_start - 1:domain_end]
            except (json.JSONDecodeError, TypeError):
                logger.warning("Bad tim_barrel_location for %s: %s", uid, tb_loc_raw)

        # Reformat exon annotations
        exon_ann_out = None
        if exon_ann_raw:
            try:
                boundaries = json.loads(exon_ann_raw)
                if isinstance(boundaries, list) and boundaries:
                    reformatted = _reformat_exon_annotations(boundaries, seq_len)
                    exon_ann_out = json.dumps(reformatted)
            except (json.JSONDecodeError, TypeError):
                logger.warning("Bad exon_annotations for %s", uid)

        try:
            conn.execute(f"""
                INSERT OR IGNORE INTO {_TABLE}
                    (uniprot_id, gene_name, sequence, domain_start, domain_end,
                     domain_sequence, exon_annotations, motif_annotations)
                VALUES (?, ?, ?, ?, ?, ?, ?, NULL)
            """, (uid, gene_name, seq, domain_start, domain_end,
                  domain_sequence, exon_ann_out))
            inserted += 1
        except sqlite3.IntegrityError as e:
            logger.warning("Skipping %s: %s", uid, e)
            skipped += 1

        if inserted % 200 == 0 and inserted > 0:
            conn.commit()
            logger.info("  %d / %d inserted", inserted, len(rows))

    conn.commit()
    logger.info("Done: %d inserted, %d skipped", inserted, skipped)
    return inserted


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def print_summary(conn: sqlite3.Connection) -> None:
    total  = conn.execute(f"SELECT COUNT(*) FROM {_TABLE}").fetchone()[0]
    w_gene = conn.execute(f"SELECT COUNT(*) FROM {_TABLE} WHERE gene_name IS NOT NULL").fetchone()[0]
    w_exon = conn.execute(f"SELECT COUNT(*) FROM {_TABLE} WHERE exon_annotations IS NOT NULL").fetchone()[0]
    w_dom  = conn.execute(f"SELECT COUNT(*) FROM {_TABLE} WHERE domain_start IS NOT NULL").fetchone()[0]

    sample = conn.execute(f"""
        SELECT exon_annotations FROM {_TABLE}
        WHERE exon_annotations IS NOT NULL LIMIT 1
    """).fetchone()
    n_exons = None
    if sample:
        try:
            n_exons = len(json.loads(sample[0]))
        except Exception:
            pass

    print(f"\n{'='*60}")
    print(f"  canonical_analysis summary")
    print(f"{'='*60}")
    print(f"  Total rows             : {total}")
    print(f"  With gene_name         : {w_gene} / {total}")
    print(f"  With domain location   : {w_dom} / {total}")
    print(f"  With exon_annotations  : {w_exon} / {total}")
    if n_exons is not None:
        print(f"  Exons in sample row    : {n_exons}")
    print(f"{'='*60}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build canonical_analysis from isoforms + proteins"
    )
    parser.add_argument("--db",       default=None)
    parser.add_argument("--rebuild",  action="store_true",
                        help="Clear existing rows and repopulate")
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()

    logging.getLogger().setLevel(getattr(logging, args.log_level.upper(), logging.INFO))

    db_path = args.db or get_config().db_path
    conn = sqlite3.connect(db_path)

    ensure_table(conn)
    build(conn, rebuild=args.rebuild)
    print_summary(conn)

    conn.close()


if __name__ == "__main__":
    main()
