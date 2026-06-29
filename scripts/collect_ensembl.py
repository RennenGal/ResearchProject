#!/usr/bin/env python3
"""
Expand TIM barrel isoform coverage using Ensembl protein-coding transcripts.

For every canonical protein in proteins that has an Ensembl mapping
(ENST stored in isoforms.ensembl_transcript_id), this script:

  1. Resolves ENST → ENSG (gene ID) via the Ensembl REST API.
  2. Fetches all protein-coding transcripts for the gene.
  3. Downloads the translated amino-acid sequence for each transcript.
  4. Stores new transcripts in ensembl_transcripts, flagging duplicates
     (sequences already present in isoforms).
  5. Runs the sliding-window alignment analysis against the canonical
     TIM barrel sequence and stores AS-affected hits in ensembl_affected.

Results are written to the two new tables only — isoforms is not touched.

Usage:
    python scripts/collect_ensembl.py
    python scripts/collect_ensembl.py --limit 50        # first N proteins (for testing)
    python scripts/collect_ensembl.py --skip-analysis   # collect only, no alignment
    python scripts/collect_ensembl.py --rebuild         # drop and re-collect everything
"""

import json
import logging
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from protein_data_collector.analysis.tim_barrel_alignment import (
    find_tim_barrel_span,
    sliding_window_align,
    _IDENTITY_MIN,
    _IDENTITY_MAX,
)
from protein_data_collector.api.ensembl_client import (
    ensg_for_enst,
    transcripts_for_gene,
    protein_sequence,
    transcript_exon_boundaries,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)

_ISOFORM_TABLE     = "isoforms"
_TRANSCRIPT_TABLE  = "ensembl_transcripts"
_AFFECTED_TABLE    = "ensembl_affected"
_MIN_SEQ_LEN       = 200   # fragment threshold (consistent with isoforms table)


# ---------------------------------------------------------------------------
# Step 1: build ENST → ENSG mapping for all proteins
# ---------------------------------------------------------------------------

def _build_ensg_map(conn: sqlite3.Connection) -> dict[str, tuple[str, str]]:
    """
    Return {uniprot_id: (enst_id, ensg_id)} for canonical proteins that
    have an ensembl_transcript_id set.  ENSG is resolved via the Ensembl API.
    """
    rows = conn.execute(
        f"SELECT iso.uniprot_id, iso.ensembl_transcript_id FROM {_ISOFORM_TABLE} iso "
        f"JOIN proteins p ON p.uniprot_id = iso.uniprot_id "
        f"WHERE iso.is_canonical=1 AND iso.ensembl_transcript_id IS NOT NULL "
        f"  AND p.canonical_uniprot_id IS NULL"
    ).fetchall()
    logger.info("Found %d proteins with an Ensembl transcript ID", len(rows))

    mapping: dict[str, tuple[str, str]] = {}
    for i, (uid, enst_raw) in enumerate(rows, 1):
        enst = enst_raw.split(".")[0]
        ensg = ensg_for_enst(enst)
        if ensg:
            mapping[uid] = (enst, ensg)
        else:
            logger.debug("Could not resolve ENSG for %s (ENST %s)", uid, enst)
        if i % 100 == 0:
            logger.info("  Resolved %d / %d gene IDs", i, len(rows))

    logger.info("ENSG resolved for %d / %d proteins", len(mapping), len(rows))
    return mapping


# ---------------------------------------------------------------------------
# Step 2+3: fetch transcripts and sequences
# ---------------------------------------------------------------------------

def _existing_uniprot_sequences(conn: sqlite3.Connection) -> dict[str, str]:
    """Return {sequence: isoform_id} for all sequences in isoforms."""
    rows = conn.execute(
        f"SELECT sequence, isoform_id FROM {_ISOFORM_TABLE} WHERE sequence IS NOT NULL"
    ).fetchall()
    return {seq: iso_id for seq, iso_id in rows}


def _existing_enst_sequences(conn: sqlite3.Connection) -> dict[str, str]:
    """Return {sequence: enst_id} for all sequences already in ensembl_transcripts."""
    rows = conn.execute(
        f"SELECT sequence, enst_id FROM {_TRANSCRIPT_TABLE} WHERE sequence IS NOT NULL"
    ).fetchall()
    return {seq: enst_id for seq, enst_id in rows}


def _existing_enst_ids(conn: sqlite3.Connection) -> set[str]:
    rows = conn.execute(f"SELECT enst_id FROM {_TRANSCRIPT_TABLE}").fetchall()
    return {r[0] for r in rows}


def collect_transcripts(
    conn: sqlite3.Connection,
    ensg_map: dict[str, tuple[str, str]],
    uniprot_seqs: dict[str, str],
    enst_seqs: dict[str, str],
    existing_ensts: set[str],
) -> int:
    """Fetch all transcripts and upsert into ensembl_transcripts. Returns insert count."""
    inserted = 0
    total_proteins = len(ensg_map)

    for i, (uid, (canonical_enst, ensg)) in enumerate(ensg_map.items(), 1):
        transcripts = transcripts_for_gene(ensg)
        if not transcripts:
            logger.debug("No transcripts returned for %s (ENSG %s)", uid, ensg)
            continue

        row = conn.execute("SELECT gene_name FROM proteins WHERE uniprot_id=?", (uid,)).fetchone()
        gene_name = row[0] if row else None

        for tx in transcripts:
            enst = tx["enst_id"]
            if enst in existing_ensts:
                continue

            seq = protein_sequence(enst)
            if not seq:
                continue

            seq_len     = len(seq)
            is_frag     = 1 if seq_len < _MIN_SEQ_LEN else 0
            dup_iso_id  = uniprot_seqs.get(seq)          # UniProt isoform duplicate
            dup_enst_id = enst_seqs.get(seq) if not dup_iso_id else None  # Ensembl internal duplicate

            conn.execute(f"""
                INSERT OR IGNORE INTO {_TRANSCRIPT_TABLE}
                    (enst_id, ensg_id, ensp_id, uniprot_id, gene_name,
                     sequence, sequence_length, is_fragment, is_mane_select,
                     biotype, duplicate_isoform_id, duplicate_enst_id)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                enst, ensg, tx.get("ensp_id"), uid, gene_name,
                seq, seq_len, is_frag, tx["is_mane_select"],
                tx["biotype"], dup_iso_id, dup_enst_id,
            ))
            existing_ensts.add(enst)
            if seq not in uniprot_seqs and seq not in enst_seqs:
                enst_seqs[seq] = enst  # register as representative for subsequent duplicates
            inserted += 1

        conn.commit()

        if i % 50 == 0 or i == total_proteins:
            logger.info("  Processed %d / %d proteins (%d transcripts inserted so far)",
                        i, total_proteins, inserted)

    return inserted


# ---------------------------------------------------------------------------
# Step 4: alignment analysis
# ---------------------------------------------------------------------------

def run_alignment_analysis(conn: sqlite3.Connection) -> tuple[int, int, int]:
    """
    Align all non-duplicate, non-fragment Ensembl transcripts against their
    protein's canonical TIM barrel sequence.

    Returns (affected, skipped_identical, skipped_absent).
    """
    conn.execute(f"DELETE FROM {_AFFECTED_TABLE}")
    conn.commit()

    rows = conn.execute(f"""
        SELECT
            et.enst_id,
            et.uniprot_id,
            et.sequence,
            et.sequence_length,
            can.tim_barrel_sequence   AS canonical_tb_seq,
            can.tim_barrel_location   AS canonical_tb_loc,
            can.sequence              AS canonical_sequence
        FROM {_TRANSCRIPT_TABLE} et
        JOIN {_ISOFORM_TABLE} can
          ON  can.uniprot_id   = et.uniprot_id
          AND can.is_canonical = 1
        WHERE et.duplicate_isoform_id IS NULL
          AND et.duplicate_enst_id IS NULL
          AND et.is_fragment = 0
          AND can.tim_barrel_sequence IS NOT NULL
    """).fetchall()

    logger.info("Running alignment analysis on %d non-duplicate Ensembl transcripts", len(rows))

    affected = skipped_identical = skipped_absent = 0

    for row in rows:
        enst_id    = row[0]
        uid        = row[1]
        iso_seq    = row[2]
        tb_seq     = row[4]
        tb_loc_raw = row[5]
        can_seq    = row[6]

        if not iso_seq or not tb_seq:
            continue

        tb_len = len(tb_seq)
        score, win_start, win_end = sliding_window_align(tb_seq, iso_seq)
        identity  = score / tb_len if tb_len > 0 else 0.0
        loc_start = win_start
        loc_end   = win_end
        span_len  = tb_len
        insertion = 0

        if identity >= _IDENTITY_MAX:
            tb_loc_dict  = json.loads(tb_loc_raw)
            can_tb_start = tb_loc_dict["start"]
            can_tb_end   = tb_loc_dict["end"]

            span = find_tim_barrel_span(can_seq, can_tb_start, can_tb_end, iso_seq)
            if span is not None:
                span_start, span_end = span
                detected_span_len = span_end - span_start + 1
                if detected_span_len > tb_len:
                    identity  = score / detected_span_len
                    loc_start = span_start
                    loc_end   = span_end
                    span_len  = detected_span_len
                    insertion = 1
                else:
                    skipped_identical += 1
                    continue
            else:
                skipped_identical += 1
                continue

        if identity < _IDENTITY_MIN:
            skipped_absent += 1
            continue

        domain_loc = json.dumps({
            "start":  loc_start,
            "end":    loc_end,
            "length": span_len,
            "source": "local_alignment" if span_len == tb_len else "local_alignment_span",
        })
        domain_seq = iso_seq[loc_start - 1:loc_end]

        conn.execute(f"""
            INSERT INTO {_AFFECTED_TABLE}
                (enst_id, uniprot_id, domain_location, domain_sequence,
                 canonical_domain_location, canonical_domain_sequence,
                 alignment_identity, alignment_score, insertion_detected)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            enst_id, uid, domain_loc, domain_seq,
            tb_loc_raw, tb_seq,
            round(identity * 100, 2), score, insertion,
        ))
        affected += 1

    conn.commit()
    logger.info(
        "Analysis done: %d AS-affected | %d skipped (>=95%%) | %d skipped (<12.5%%)",
        affected, skipped_identical, skipped_absent,
    )
    return affected, skipped_identical, skipped_absent


# ---------------------------------------------------------------------------
# Merge A: backfill_exons logic (merged from backfill_exons.py)
# ---------------------------------------------------------------------------

_BACKFILL_TRANSCRIPT_TABLE = "ensembl_transcripts"
_BACKFILL_AFFECTED_TABLE   = "ensembl_affected"


def _ensure_exon_columns(conn: sqlite3.Connection) -> None:
    existing_enst = {r[1] for r in conn.execute(f"PRAGMA table_info({_BACKFILL_TRANSCRIPT_TABLE})")}
    if "exon_annotations" not in existing_enst:
        conn.execute(f"ALTER TABLE {_BACKFILL_TRANSCRIPT_TABLE} ADD COLUMN exon_annotations TEXT")
        logger.info("Added exon_annotations column to %s", _BACKFILL_TRANSCRIPT_TABLE)
    if "duplicate_enst_id" not in existing_enst:
        conn.execute(f"ALTER TABLE {_BACKFILL_TRANSCRIPT_TABLE} ADD COLUMN duplicate_enst_id TEXT")
        logger.info("Added duplicate_enst_id column to %s", _BACKFILL_TRANSCRIPT_TABLE)

    existing_aff = {r[1] for r in conn.execute(f"PRAGMA table_info({_BACKFILL_AFFECTED_TABLE})")}
    if "exon_boundary_in_domain" not in existing_aff:
        conn.execute(
            f"ALTER TABLE {_BACKFILL_AFFECTED_TABLE} "
            f"ADD COLUMN exon_boundary_in_domain INTEGER NOT NULL DEFAULT 0"
        )
        logger.info("Added exon_boundary_in_domain column to %s", _BACKFILL_AFFECTED_TABLE)
    if "exon_boundaries_in_domain_count" not in existing_aff:
        conn.execute(
            f"ALTER TABLE {_BACKFILL_AFFECTED_TABLE} "
            f"ADD COLUMN exon_boundaries_in_domain_count INTEGER NOT NULL DEFAULT 0"
        )
        logger.info("Added exon_boundaries_in_domain_count column to %s", _BACKFILL_AFFECTED_TABLE)

    conn.commit()


def _backfill_exon_annotations(conn: sqlite3.Connection) -> int:
    """Fetch exon boundaries from Ensembl for transcripts missing exon_annotations."""
    rows = conn.execute(
        f"SELECT enst_id FROM {_BACKFILL_TRANSCRIPT_TABLE} WHERE exon_annotations IS NULL"
    ).fetchall()
    total = len(rows)
    logger.info("Fetching exon boundaries for %d transcripts", total)

    updated = 0
    for i, (enst_id,) in enumerate(rows, 1):
        boundaries = transcript_exon_boundaries(enst_id)
        annotation = json.dumps(boundaries)
        conn.execute(
            f"UPDATE {_BACKFILL_TRANSCRIPT_TABLE} SET exon_annotations=? WHERE enst_id=?",
            (annotation, enst_id),
        )
        updated += 1
        if i % 100 == 0 or i == total:
            conn.commit()
            logger.info("  %d / %d transcripts processed", i, total)

    conn.commit()
    return updated


def _flag_exon_boundary_in_domain(conn: sqlite3.Connection) -> tuple[int, int]:
    """Flag exon boundaries that fall inside the domain in ensembl_affected."""
    rows = conn.execute(f"""
        SELECT
            aff.id,
            aff.domain_location,
            et.exon_annotations
        FROM {_BACKFILL_AFFECTED_TABLE} aff
        JOIN {_BACKFILL_TRANSCRIPT_TABLE} et ON et.enst_id = aff.enst_id
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
            f"UPDATE {_BACKFILL_AFFECTED_TABLE} "
            f"SET exon_boundary_in_domain=?, exon_boundaries_in_domain_count=? WHERE id=?",
            (flag, count, row_id),
        )

    conn.commit()
    logger.info("Flagged %d / %d AS-affected transcripts as exon_boundary_in_domain=1",
                flagged, len(rows))
    return flagged, len(rows)


# ---------------------------------------------------------------------------
# Pipeline entry point
# ---------------------------------------------------------------------------

def run(db_path: str) -> None:
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row

    # Step 1: resolve ENSG IDs
    ensg_map = _build_ensg_map(conn)

    # Step 2+3: collect transcripts
    uniprot_seqs   = _existing_uniprot_sequences(conn)
    enst_seqs      = _existing_enst_sequences(conn)
    existing_ensts = _existing_enst_ids(conn)
    logger.info("Collecting Ensembl transcripts (%d UniProt sequences, %d Ensembl sequences cached)",
                len(uniprot_seqs), len(enst_seqs))
    inserted = collect_transcripts(conn, ensg_map, uniprot_seqs, enst_seqs, existing_ensts)

    # Summary so far
    total_enst  = conn.execute(f"SELECT COUNT(*) FROM {_TRANSCRIPT_TABLE}").fetchone()[0]
    dup_uniprot = conn.execute(f"SELECT COUNT(*) FROM {_TRANSCRIPT_TABLE} WHERE duplicate_isoform_id IS NOT NULL").fetchone()[0]
    dup_enst    = conn.execute(f"SELECT COUNT(*) FROM {_TRANSCRIPT_TABLE} WHERE duplicate_enst_id IS NOT NULL").fetchone()[0]
    dup_count   = dup_uniprot + dup_enst
    frag_count  = conn.execute(f"SELECT COUNT(*) FROM {_TRANSCRIPT_TABLE} WHERE is_fragment=1 AND duplicate_isoform_id IS NULL AND duplicate_enst_id IS NULL").fetchone()[0]
    novel_count = total_enst - dup_count

    print(f"\n{'='*60}")
    print(f"  Ensembl transcript collection")
    print(f"{'='*60}")
    print(f"  Proteins processed            : {len(ensg_map)}")
    print(f"  Transcripts inserted          : {inserted}")
    print(f"  Total in {_TRANSCRIPT_TABLE:<20}: {total_enst}")
    print(f"  Duplicates — same seq as UniProt: {dup_uniprot}")
    print(f"  Duplicates — internal Ensembl   : {dup_enst}")
    print(f"  Fragments (< {_MIN_SEQ_LEN} aa, novel)    : {frag_count}")
    print(f"  Novel unique transcripts        : {novel_count}")

    # Step 4: alignment analysis
    print(f"\n  Running alignment analysis...")
    affected, skipped_id, skipped_ab = run_alignment_analysis(conn)

    print(f"\n{'='*60}")
    print(f"  Alignment analysis")
    print(f"{'='*60}")
    print(f"  Novel non-fragment transcripts: {novel_count - frag_count}")
    print(f"  AS-affected (12.5%–95%)       : {affected}")
    print(f"  Skipped (>= 95%% identical)   : {skipped_id}")
    print(f"  Skipped (< 12.5%% identity)   : {skipped_ab}")
    print(f"{'='*60}")

    # Merge A: backfill exon boundary data for Ensembl transcripts
    logger.info("=== Backfilling exon boundary data for Ensembl transcripts ===")
    conn.row_factory = None
    _ensure_exon_columns(conn)

    logger.info("=== Phase 1: fetching exon annotations ===")
    updated = _backfill_exon_annotations(conn)
    print(f"\n  Exon annotations written: {updated}")

    logger.info("=== Phase 2: flagging exon boundaries in domain ===")
    flagged, total_eval = _flag_exon_boundary_in_domain(conn)

    total_aff = conn.execute(f"SELECT COUNT(*) FROM {_BACKFILL_AFFECTED_TABLE}").fetchone()[0]
    print(f"\n{'='*60}")
    print(f"  Exon boundary analysis — TIM barrel (Homo sapiens)")
    print(f"{'='*60}")
    print(f"  AS-affected transcripts (total)          : {total_aff}")
    print(f"  Evaluated (exon data available)          : {total_eval}")
    print(f"  Exon boundary falls inside domain        : {flagged}")
    if total_eval > 0:
        print(f"  Fraction with intra-domain junction      : {flagged/total_eval:.1%}")
    print(f"{'='*60}")

    conn.close()
