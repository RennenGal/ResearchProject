"""
TIM barrel local alignment analysis.

For each alternative (non-canonical), non-fragment isoform, slides the canonical
TIM barrel sequence along the isoform sequence using ungapped alignment
(match = 1, mismatch = 0) to find the best-matching window.

Insertion detection:
    For isoforms that score >= 95% on the window alignment (TIM barrel apparently
    identical), the conserved flanking sequences (20 aa on each side of the canonical
    TIM barrel) are searched in the alternative isoform.  If both flanks are found
    exactly, the span between them measures the actual TIM barrel region length in the
    alternative.  A span longer than tb_len means extra sequence was inserted into the
    TIM barrel region; identity is recomputed as window_score / span_len.

Inclusion thresholds:
    identity >= 12.5%   — at least one beta-alpha motif is present
    identity <  95%     — meaningful AS effect on the TIM barrel
"""

import json
import logging
import sqlite3
from dataclasses import dataclass
from typing import Optional

logger = logging.getLogger(__name__)

# Thresholds (fraction, not percent)
_IDENTITY_MIN = 0.125   # 12.5% — one beta-alpha motif out of 8
_IDENTITY_MAX = 0.95    # 95%   — meaningful AS effect cutoff

_FLANK_LEN = 20         # aa of conserved context used for span detection


# ---------------------------------------------------------------------------
# Alignment helpers
# ---------------------------------------------------------------------------

def sliding_window_align(tim_barrel_seq: str, isoform_seq: str) -> tuple[int, int, int]:
    """
    Ungapped local alignment of ``tim_barrel_seq`` against ``isoform_seq``.

    Slides a window of len(tim_barrel_seq) along isoform_seq and scores each
    position with match=1, mismatch=0.

    Returns
    -------
    (best_score, start_1based, end_1based)
        Coordinates are 1-based, inclusive, in isoform space.
        Returns (0, 0, 0) if isoform_seq is shorter than tim_barrel_seq.
    """
    L = len(tim_barrel_seq)
    n = len(isoform_seq)

    if n < L:
        return 0, 0, 0

    best_score = -1
    best_start = 0

    for i in range(n - L + 1):
        score = sum(
            1 for j in range(L) if tim_barrel_seq[j] == isoform_seq[i + j]
        )
        if score > best_score:
            best_score = score
            best_start = i

    return best_score, best_start + 1, best_start + L  # 1-based, inclusive


def _find_exact(query: str, target: str) -> int:
    """
    Return the 0-based start position of the first exact occurrence of
    ``query`` in ``target``, or -1 if not found.
    """
    if not query or len(query) > len(target):
        return -1
    for i in range(len(target) - len(query) + 1):
        if target[i:i + len(query)] == query:
            return i
    return -1


def find_tim_barrel_span(
    canonical_seq: str,
    can_tb_start: int,   # 1-based
    can_tb_end: int,     # 1-based
    alt_seq: str,
    flank_len: int = _FLANK_LEN,
) -> Optional[tuple[int, int]]:
    """
    Use conserved flanking sequences to locate the full TIM barrel region
    (including any inserted sequence) in the alternative isoform.

    Extracts ``flank_len`` aa immediately before and after the canonical TIM
    barrel, searches for each as an exact substring in the alternative, and
    returns the span between them.

    Returns
    -------
    (span_start_1based, span_end_1based)  or  None if either flank is missing
    or the located span is malformed (end < start).
    """
    n_flank = canonical_seq[max(0, can_tb_start - 1 - flank_len): can_tb_start - 1]
    c_flank = canonical_seq[can_tb_end: can_tb_end + flank_len]

    # Need at least 5 aa of context on each side
    if len(n_flank) < 5 or len(c_flank) < 5:
        return None

    n_pos = _find_exact(n_flank, alt_seq)   # 0-based start of N-flank in alt
    c_pos = _find_exact(c_flank, alt_seq)   # 0-based start of C-flank in alt

    if n_pos == -1 or c_pos == -1:
        return None

    # Span: immediately after N-flank up to immediately before C-flank (1-based)
    span_start = n_pos + len(n_flank) + 1
    span_end   = c_pos                       # c_flank starts at c_pos+1 (1-based)

    if span_end < span_start:
        return None

    return span_start, span_end


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class AlignmentResult:
    isoform_id: str
    uniprot_id: str
    sequence: str
    sequence_length: int
    is_fragment: int
    exon_count: Optional[int]
    exon_annotations: Optional[str]
    splice_variants: Optional[str]
    ensembl_transcript_id: Optional[str]
    alphafold_id: Optional[str]
    # Alignment-derived
    domain_location: str            # JSON — span in this isoform
    domain_sequence: str            # isoform subsequence at span
    canonical_domain_location: str  # JSON — original canonical location
    canonical_domain_sequence: str  # canonical domain sequence used as query
    identity_percentage: float
    alignment_score: int


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def build_tim_barrel_isoforms(
    conn: sqlite3.Connection,
    isoform_table: str = "tb_isoforms",
) -> tuple[list, int, int, int]:
    """
    Run the alignment analysis and return
    (results, skipped_identical, skipped_absent, insertions_detected).

    Queries ``isoform_table`` for all alternative, non-fragment isoforms
    whose protein has a canonical isoform with a known tim_barrel_sequence.
    """
    conn.row_factory = sqlite3.Row

    rows = conn.execute(f"""
        SELECT
            alt.isoform_id,
            alt.uniprot_id,
            alt.sequence,
            alt.sequence_length,
            alt.is_fragment,
            alt.exon_count,
            alt.exon_annotations,
            alt.splice_variants,
            alt.ensembl_transcript_id,
            alt.alphafold_id,
            can.sequence              AS canonical_sequence,
            can.tim_barrel_sequence   AS canonical_tb_seq,
            can.tim_barrel_location   AS canonical_tb_loc
        FROM {isoform_table} alt
        JOIN {isoform_table} can
          ON  can.uniprot_id   = alt.uniprot_id
          AND can.is_canonical = 1
        WHERE alt.is_canonical = 0
          AND alt.is_fragment  = 0
          AND can.tim_barrel_sequence IS NOT NULL
        ORDER BY alt.uniprot_id, alt.isoform_id
    """).fetchall()

    logger.info("Aligning %d alternative isoforms against their canonical TIM barrel", len(rows))

    results: list[AlignmentResult] = []
    skipped_identical  = 0
    skipped_absent     = 0
    insertions_detected = 0

    for row in rows:
        tb_seq  = row["canonical_tb_seq"]
        iso_seq = row["sequence"]
        tb_len  = len(tb_seq)

        score, win_start, win_end = sliding_window_align(tb_seq, iso_seq)

        identity  = score / tb_len if tb_len > 0 else 0.0
        loc_start = win_start
        loc_end   = win_end
        span_len  = tb_len

        if identity >= _IDENTITY_MAX:
            # Window says TIM barrel is essentially unchanged.
            # Check whether an exon was inserted into the TIM barrel region
            # using conserved flanking sequences to measure the actual span.
            tb_loc_dict  = json.loads(row["canonical_tb_loc"])
            can_tb_start = tb_loc_dict["start"]
            can_tb_end   = tb_loc_dict["end"]

            span = find_tim_barrel_span(
                row["canonical_sequence"], can_tb_start, can_tb_end, iso_seq
            )

            if span is not None:
                span_start, span_end = span
                detected_span_len = span_end - span_start + 1
                if detected_span_len > tb_len:
                    # Insertion confirmed: recompute identity over the full span
                    identity  = score / detected_span_len
                    loc_start = span_start
                    loc_end   = span_end
                    span_len  = detected_span_len
                    insertions_detected += 1
                    # Fall through to threshold checks below
                else:
                    skipped_identical += 1
                    continue
            else:
                # Flanks not conserved — cannot determine span; treat as identical
                skipped_identical += 1
                continue

        if identity < _IDENTITY_MIN:
            skipped_absent += 1
            continue

        tb_location = json.dumps({
            "start":  loc_start,
            "end":    loc_end,
            "length": span_len,
            "source": "local_alignment" if span_len == tb_len else "local_alignment_span",
        })
        tb_subsequence = iso_seq[loc_start - 1:loc_end]

        results.append(AlignmentResult(
            isoform_id=row["isoform_id"],
            uniprot_id=row["uniprot_id"],
            sequence=iso_seq,
            sequence_length=row["sequence_length"],
            is_fragment=row["is_fragment"],
            exon_count=row["exon_count"],
            exon_annotations=row["exon_annotations"],
            splice_variants=row["splice_variants"],
            ensembl_transcript_id=row["ensembl_transcript_id"],
            alphafold_id=row["alphafold_id"],
            domain_location=tb_location,
            domain_sequence=tb_subsequence,
            canonical_domain_location=row["canonical_tb_loc"],
            canonical_domain_sequence=tb_seq,
            identity_percentage=round(identity * 100, 2),
            alignment_score=score,
        ))

    logger.info(
        "Results: %d inserted (%d by insertion detection) | "
        "%d skipped (>= 95%%) | %d skipped (< 12.5%%)",
        len(results), insertions_detected, skipped_identical, skipped_absent,
    )
    return results, skipped_identical, skipped_absent, insertions_detected


def populate_tim_barrel_isoforms(
    conn: sqlite3.Connection,
    isoform_table: str = "tb_isoforms",
    output_table: str = "tb_affected_isoforms",
) -> tuple[int, int, int, int]:
    """
    Rebuild ``output_table`` from scratch and return
    (inserted, skipped_identical, skipped_absent, insertions_detected).
    """
    conn.execute(f"DELETE FROM {output_table}")

    results, skipped_identical, skipped_absent, insertions_detected = build_tim_barrel_isoforms(
        conn, isoform_table=isoform_table
    )

    conn.executemany(f"""
        INSERT OR REPLACE INTO {output_table} (
            isoform_id, uniprot_id, is_canonical, sequence, sequence_length,
            is_fragment, exon_count, exon_annotations, splice_variants,
            domain_location, domain_sequence,
            canonical_domain_location, canonical_domain_sequence,
            identity_percentage, alignment_score,
            ensembl_transcript_id, alphafold_id
        ) VALUES (
            :isoform_id, :uniprot_id, 0, :sequence, :sequence_length,
            :is_fragment, :exon_count, :exon_annotations, :splice_variants,
            :domain_location, :domain_sequence,
            :canonical_domain_location, :canonical_domain_sequence,
            :identity_percentage, :alignment_score,
            :ensembl_transcript_id, :alphafold_id
        )
    """, [r.__dict__ for r in results])

    conn.commit()
    return len(results), skipped_identical, skipped_absent, insertions_detected
