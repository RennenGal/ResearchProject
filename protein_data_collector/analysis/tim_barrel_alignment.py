"""
TIM barrel AS-affected isoform detection.

Primary method: VSP-based overlap
----------------------------------
Each alternative isoform in UniProt carries one or more VSP (Variable
Sequence Position) feature annotations that describe exactly which residues
of the canonical sequence are changed.  An isoform is "domain-affecting" when
at least one VSP overlaps (by any number of residues) with the canonical TIM
barrel domain [domain_start, domain_end].

This is more principled than a sliding-window alignment because:
  - It uses the explicitly annotated splice boundaries rather than inferring
    them from sequence similarity.
  - It detects small loop changes (~3 aa) that a similarity threshold would
    silently exclude.
  - It does not require arbitrary identity thresholds.

Overlap stored per affected VSP:
  feature_id, can_start, can_end (canonical coords),
  overlap_start, overlap_end (intersection with domain),
  overlap_residues, overlap_fraction (overlap_residues / domain_len)

Fallback: sliding-window alignment
------------------------------------
For the 6 isoforms whose splice_variants list is empty (no VSP data), the
original ungapped sliding-window is used with thresholds 12.5% < id < 95%.
"""

import json
import logging
import sqlite3
from dataclasses import dataclass, field
from typing import Optional

logger = logging.getLogger(__name__)

# Sliding-window fallback thresholds (fraction)
_IDENTITY_MIN = 0.125
_IDENTITY_MAX = 0.95
_FLANK_LEN    = 20


# ---------------------------------------------------------------------------
# VSP overlap detection
# ---------------------------------------------------------------------------

def _vsp_domain_events(
    vsps: list[dict],
    domain_start: int,
    domain_end: int,
) -> list[dict]:
    """
    Return the subset of VSPs whose canonical coordinates overlap with
    [domain_start, domain_end].  Each returned dict includes precomputed
    overlap metrics.
    """
    domain_len = domain_end - domain_start + 1
    overlapping = []
    for v in vsps:
        vs = v["location"]["start"]["value"]
        ve = v["location"]["end"]["value"]
        if vs <= domain_end and ve >= domain_start:
            ov_start = max(vs, domain_start)
            ov_end   = min(ve, domain_end)
            ov_res   = ov_end - ov_start + 1
            overlapping.append({
                "feature_id":       v["featureId"],
                "can_start":        vs,
                "can_end":          ve,
                "overlap_start":    ov_start,
                "overlap_end":      ov_end,
                "overlap_residues": ov_res,
                "overlap_fraction": round(ov_res / domain_len, 3),
            })
    return overlapping


# ---------------------------------------------------------------------------
# Sliding-window helpers (fallback + domain location in alt isoform)
# ---------------------------------------------------------------------------

def sliding_window_align(tim_barrel_seq: str, isoform_seq: str) -> tuple[int, int, int]:
    """
    Ungapped local alignment.  Returns (best_score, start_1based, end_1based).
    Returns (0, 0, 0) if isoform_seq is shorter than tim_barrel_seq.
    """
    L = len(tim_barrel_seq)
    n = len(isoform_seq)
    if n < L:
        return 0, 0, 0
    best_score = -1
    best_start = 0
    for i in range(n - L + 1):
        score = sum(1 for j in range(L) if tim_barrel_seq[j] == isoform_seq[i + j])
        if score > best_score:
            best_score = score
            best_start = i
    return best_score, best_start + 1, best_start + L


def _find_exact(query: str, target: str) -> int:
    if not query or len(query) > len(target):
        return -1
    for i in range(len(target) - len(query) + 1):
        if target[i:i + len(query)] == query:
            return i
    return -1


def find_tim_barrel_span(
    canonical_seq: str,
    can_tb_start: int,
    can_tb_end: int,
    alt_seq: str,
    flank_len: int = _FLANK_LEN,
) -> Optional[tuple[int, int]]:
    n_flank = canonical_seq[max(0, can_tb_start - 1 - flank_len): can_tb_start - 1]
    c_flank = canonical_seq[can_tb_end: can_tb_end + flank_len]
    if len(n_flank) < 5 or len(c_flank) < 5:
        return None
    n_pos = _find_exact(n_flank, alt_seq)
    c_pos = _find_exact(c_flank, alt_seq)
    if n_pos == -1 or c_pos == -1:
        return None
    span_start = n_pos + len(n_flank) + 1
    span_end   = c_pos
    if span_end < span_start:
        return None
    return span_start, span_end


def _sliding_window_result(
    tb_seq: str, iso_seq: str, can_seq: str, can_tb_start: int, can_tb_end: int
) -> tuple[float, Optional[str], Optional[str]]:
    """
    Compute sliding-window identity, domain_location JSON, and domain_sequence
    for one alternative isoform.  Returns (identity, loc_json_or_None, subseq_or_None).
    """
    tb_len = len(tb_seq)
    score, win_start, win_end = sliding_window_align(tb_seq, iso_seq)
    identity = score / tb_len if tb_len > 0 else 0.0

    if identity >= _IDENTITY_MAX:
        span = find_tim_barrel_span(can_seq, can_tb_start, can_tb_end, iso_seq)
        if span is not None:
            span_start, span_end = span
            span_len = span_end - span_start + 1
            if span_len > tb_len:
                identity = score / span_len
                win_start, win_end = span_start, span_end

    if identity < 0.05:
        return identity, None, None

    loc = json.dumps({
        "start":  win_start,
        "end":    win_end,
        "length": win_end - win_start + 1,
        "source": "sliding_window",
    })
    seq = iso_seq[win_start - 1: win_end]
    return identity, loc, seq


# ---------------------------------------------------------------------------
# Data class
# ---------------------------------------------------------------------------

@dataclass
class AlignmentResult:
    isoform_id:                str
    uniprot_id:                str
    sequence:                  str
    sequence_length:           int
    is_fragment:               int
    exon_count:                Optional[int]
    exon_annotations:          Optional[str]
    splice_variants:           Optional[str]
    ensembl_transcript_id:     Optional[str]
    alphafold_id:              Optional[str]
    domain_location:           Optional[str]   # JSON, alt isoform coords
    domain_sequence:           Optional[str]
    canonical_domain_location: str             # JSON, canonical coords
    canonical_domain_sequence: str
    identity_percentage:       float
    alignment_score:           int
    vsp_domain_events:         str             # JSON list of overlapping VSPs
    detection_method:          str             # 'vsp' or 'sliding_window'


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def build_tim_barrel_isoforms(
    conn: sqlite3.Connection,
    isoform_table: str = "isoforms",
) -> tuple[list, int, int]:
    """
    Detect domain-affecting isoforms using VSP overlap as the primary method.

    Returns
    -------
    (results, skipped_no_overlap, skipped_fallback_absent)
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

    logger.info("Checking %d alternative isoforms for domain-affecting events", len(rows))

    results: list[AlignmentResult] = []
    skipped_no_overlap   = 0
    skipped_fallback     = 0

    for row in rows:
        tb_seq    = row["canonical_tb_seq"]
        iso_seq   = row["sequence"]
        can_seq   = row["canonical_sequence"]
        tb_loc    = json.loads(row["canonical_tb_loc"])
        can_ds    = tb_loc["start"]
        can_de    = tb_loc["end"]
        tb_len    = len(tb_seq)

        sv_raw = row["splice_variants"]
        vsps   = json.loads(sv_raw) if sv_raw else []

        # --- Primary: VSP-based detection ---
        if vsps:
            events = _vsp_domain_events(vsps, can_ds, can_de)
            if not events:
                skipped_no_overlap += 1
                continue

            identity, loc_json, dom_subseq = _sliding_window_result(
                tb_seq, iso_seq, can_seq, can_ds, can_de
            )
            # identity_percentage stored for reference; detection is VSP-driven
            score = round(identity * tb_len)
            method = "vsp"

        # --- Fallback: sliding-window for isoforms with no VSP annotations ---
        else:
            score, win_start, win_end = sliding_window_align(tb_seq, iso_seq)
            identity = score / tb_len if tb_len > 0 else 0.0
            span_len = tb_len

            if identity >= _IDENTITY_MAX:
                span = find_tim_barrel_span(can_seq, can_ds, can_de, iso_seq)
                if span:
                    span_start, span_end = span
                    detected_span_len = span_end - span_start + 1
                    if detected_span_len > tb_len:
                        identity = score / detected_span_len
                        win_start, win_end = span_start, span_end
                        span_len = detected_span_len

                if identity >= _IDENTITY_MAX:
                    skipped_no_overlap += 1
                    continue

            if identity < _IDENTITY_MIN:
                skipped_fallback += 1
                continue

            loc_json = json.dumps({
                "start": win_start, "end": win_end,
                "length": span_len, "source": "sliding_window",
            })
            dom_subseq = iso_seq[win_start - 1: win_end]
            events = []
            method = "sliding_window"

        results.append(AlignmentResult(
            isoform_id=row["isoform_id"],
            uniprot_id=row["uniprot_id"],
            sequence=iso_seq,
            sequence_length=row["sequence_length"],
            is_fragment=row["is_fragment"],
            exon_count=row["exon_count"],
            exon_annotations=row["exon_annotations"],
            splice_variants=sv_raw,
            ensembl_transcript_id=row["ensembl_transcript_id"],
            alphafold_id=row["alphafold_id"],
            domain_location=loc_json,
            domain_sequence=dom_subseq,
            canonical_domain_location=row["canonical_tb_loc"],
            canonical_domain_sequence=tb_seq,
            identity_percentage=round(identity * 100, 2),
            alignment_score=score,
            vsp_domain_events=json.dumps(events),
            detection_method=method,
        ))

    logger.info(
        "Results: %d domain-affecting isoforms | %d no overlap | %d fallback absent",
        len(results), skipped_no_overlap, skipped_fallback,
    )
    return results, skipped_no_overlap, skipped_fallback


def populate_tim_barrel_isoforms(
    conn: sqlite3.Connection,
    isoform_table: str  = "isoforms",
    output_table: str   = "affected_isoforms",
) -> tuple[int, int, int]:
    """Rebuild output_table and return (inserted, skipped_no_overlap, skipped_fallback)."""
    # Add new columns if missing
    existing = {r[1] for r in conn.execute(f"PRAGMA table_info({output_table})")}
    for col, typedef in [("vsp_domain_events", "TEXT"), ("detection_method", "TEXT")]:
        if col not in existing:
            conn.execute(f"ALTER TABLE {output_table} ADD COLUMN {col} {typedef}")
            logger.info("Added %s to %s", col, output_table)

    conn.execute(f"DELETE FROM {output_table}")
    results, skipped_no_overlap, skipped_fallback = build_tim_barrel_isoforms(
        conn, isoform_table=isoform_table
    )

    conn.executemany(f"""
        INSERT OR REPLACE INTO {output_table} (
            isoform_id, uniprot_id, is_canonical, sequence, sequence_length,
            is_fragment, exon_count, exon_annotations, splice_variants,
            domain_location, domain_sequence,
            canonical_domain_location, canonical_domain_sequence,
            identity_percentage, alignment_score,
            ensembl_transcript_id, alphafold_id,
            vsp_domain_events, detection_method
        ) VALUES (
            :isoform_id, :uniprot_id, 0, :sequence, :sequence_length,
            :is_fragment, :exon_count, :exon_annotations, :splice_variants,
            :domain_location, :domain_sequence,
            :canonical_domain_location, :canonical_domain_sequence,
            :identity_percentage, :alignment_score,
            :ensembl_transcript_id, :alphafold_id,
            :vsp_domain_events, :detection_method
        )
    """, [r.__dict__ for r in results])

    conn.commit()
    return len(results), skipped_no_overlap, skipped_fallback
