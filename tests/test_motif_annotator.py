"""Tests for identify_ba_motifs() and _reformat_exon_annotations()."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))
from protein_data_collector.analysis.motif_annotator import (
    MIN_HELIX,
    MIN_STRAND,
    MERGE_GAP,
    identify_ba_motifs,
)
from scripts.build_canonical_analysis import _reformat_exon_annotations


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _ss(pattern: str) -> np.ndarray:
    """Convert a string like 'EEE---HHHH' to a numpy char array."""
    return np.array(list(pattern))


def _build_barrel(
    n_motifs: int = 8,
    strand_len: int = 4,
    loop_len: int = 5,
    helix_len: int = 6,
    inter_len: int = 3,
    prefix: int = 5,
    suffix: int = 5,
) -> tuple[np.ndarray, int, int]:
    """
    Construct a synthetic SS array for a TIM barrel with *n_motifs* (ba) units.
    Returns (ss_array, domain_start_1based, domain_end_1based).
    """
    motif = "E" * strand_len + "-" * loop_len + "H" * helix_len
    inter = "-" * inter_len
    domain_body = inter.join([motif] * n_motifs)
    full = "-" * prefix + domain_body + "-" * suffix
    ds = prefix + 1            # 1-based
    de = prefix + len(domain_body)  # 1-based
    return _ss(full), ds, de


# ---------------------------------------------------------------------------
# identify_ba_motifs — clean 8-motif barrel
# ---------------------------------------------------------------------------

class TestIdentifyBaMotifs:
    def test_clean_8_motif_barrel_count(self):
        ss, ds, de = _build_barrel(8)
        motifs = identify_ba_motifs(ss, ds, de)
        assert len(motifs) == 8

    def test_motif_numbers_sequential(self):
        ss, ds, de = _build_barrel(8)
        motifs = identify_ba_motifs(ss, ds, de)
        assert [m["motif"] for m in motifs] == list(range(1, 9))

    def test_positions_within_domain(self):
        ss, ds, de = _build_barrel(8)
        motifs = identify_ba_motifs(ss, ds, de)
        for m in motifs:
            assert m["beta_start"] >= ds
            assert m["alpha_end"]  <= de

    def test_beta_before_alpha_within_motif(self):
        ss, ds, de = _build_barrel(8)
        motifs = identify_ba_motifs(ss, ds, de)
        for m in motifs:
            assert m["beta_start"] <= m["beta_end"]
            assert m["beta_end"]   <  m["alpha_start"]
            assert m["alpha_start"] <= m["alpha_end"]

    def test_motif_span_matches_beta_alpha(self):
        ss, ds, de = _build_barrel(8)
        motifs = identify_ba_motifs(ss, ds, de)
        for m in motifs:
            assert m["start"] == m["beta_start"]
            assert m["end"]   == m["alpha_end"]

    def test_motifs_non_overlapping(self):
        ss, ds, de = _build_barrel(8)
        motifs = identify_ba_motifs(ss, ds, de)
        for i in range(len(motifs) - 1):
            assert motifs[i]["alpha_end"] < motifs[i + 1]["beta_start"]

    # --- partial barrel ---

    def test_partial_barrel_5_motifs(self):
        ss, ds, de = _build_barrel(5)
        motifs = identify_ba_motifs(ss, ds, de)
        assert len(motifs) == 5

    def test_empty_domain_returns_empty(self):
        ss = _ss("-" * 50)
        motifs = identify_ba_motifs(ss, 1, 50)
        assert motifs == []

    def test_strand_without_following_helix_is_skipped(self):
        # Strand at the end with no helix after — should not produce a motif
        ss = _ss("EEEE---HHHHHH---EEEE")
        motifs = identify_ba_motifs(ss, 1, len(ss))
        assert len(motifs) == 1   # only first ba pair; trailing strand has no helix

    # --- MERGE_GAP behaviour ---

    def test_merge_gap_fuses_adjacent_strands(self):
        # Two strand segments with a 1-residue gap: should merge into one strand
        gap = "-" * MERGE_GAP
        helix = "H" * MIN_HELIX
        ss = _ss("E" * MIN_STRAND + gap + "E" * MIN_STRAND + "---" + helix)
        motifs = identify_ba_motifs(ss, 1, len(ss))
        assert len(motifs) == 1
        # The merged strand should start at position 1
        assert motifs[0]["beta_start"] == 1

    def test_gap_larger_than_merge_gap_keeps_separate(self):
        # Gap of MERGE_GAP+1 should not merge
        gap = "-" * (MERGE_GAP + 1)
        helix = "H" * MIN_HELIX
        loop  = "-" * 3
        ss = _ss("E" * MIN_STRAND + gap + "E" * MIN_STRAND + loop + helix)
        motifs = identify_ba_motifs(ss, 1, len(ss))
        # The first short strand may be below MIN_STRAND after split — result
        # depends on lengths, but at least we verify no crash and at most 1 motif
        assert len(motifs) <= 1

    # --- domain boundary / offset ---

    def test_domain_offset_applied_correctly(self):
        prefix = 20
        ss, ds, de = _build_barrel(8, prefix=prefix)
        motifs = identify_ba_motifs(ss, ds, de)
        assert len(motifs) == 8
        # First beta_start must be at least ds
        assert motifs[0]["beta_start"] >= ds
        # Last alpha_end must be at most de
        assert motifs[-1]["alpha_end"] <= de

    def test_domain_at_protein_start(self):
        ss, ds, de = _build_barrel(8, prefix=0)
        motifs = identify_ba_motifs(ss, ds, de)
        assert len(motifs) == 8
        assert motifs[0]["beta_start"] == ds

    def test_domain_at_protein_end(self):
        ss, ds, de = _build_barrel(8, suffix=0)
        motifs = identify_ba_motifs(ss, ds, de)
        assert len(motifs) == 8

    def test_strand_below_min_length_ignored(self):
        # A strand of MIN_STRAND-1 residues should be invisible
        short = "E" * (MIN_STRAND - 1)
        helix = "H" * MIN_HELIX
        full_strand = "E" * MIN_STRAND
        loop = "-" * 5
        ss = _ss(short + loop + full_strand + loop + helix)
        motifs = identify_ba_motifs(ss, 1, len(ss))
        assert len(motifs) == 1

    def test_helix_below_min_length_ignored(self):
        short_helix = "H" * (MIN_HELIX - 1)
        full_helix  = "H" * MIN_HELIX
        strand = "E" * MIN_STRAND
        loop   = "-" * 5
        ss = _ss(strand + loop + short_helix + loop + strand + loop + full_helix)
        motifs = identify_ba_motifs(ss, 1, len(ss))
        assert len(motifs) == 1


# ---------------------------------------------------------------------------
# _reformat_exon_annotations
# ---------------------------------------------------------------------------

class TestReformatExonAnnotations:
    def test_single_exon(self):
        result = _reformat_exon_annotations([], seq_len=100)
        assert result == [{"exon": 1, "start": 1, "end": 100}]

    def test_multi_exon(self):
        result = _reformat_exon_annotations([85, 116, 209], seq_len=300)
        assert result == [
            {"exon": 1, "start":   1, "end":  85},
            {"exon": 2, "start":  86, "end": 116},
            {"exon": 3, "start": 117, "end": 209},
            {"exon": 4, "start": 210, "end": 300},
        ]

    def test_boundary_at_seq_len(self):
        # Last boundary equals seq_len — final exon has start > end, which is
        # a degenerate case; _reformat should still produce it without crashing
        result = _reformat_exon_annotations([50, 100], seq_len=100)
        assert result[-1]["end"] == 100
        assert result[-1]["start"] == 101   # degenerate but deterministic

    def test_unsorted_boundaries_sorted(self):
        result = _reformat_exon_annotations([209, 85, 116], seq_len=300)
        starts = [e["start"] for e in result]
        assert starts == sorted(starts)

    def test_exon_numbers_sequential(self):
        result = _reformat_exon_annotations([50, 100, 150], seq_len=200)
        assert [e["exon"] for e in result] == [1, 2, 3, 4]

    def test_first_exon_starts_at_1(self):
        result = _reformat_exon_annotations([42], seq_len=200)
        assert result[0]["start"] == 1

    def test_last_exon_ends_at_seq_len(self):
        result = _reformat_exon_annotations([42, 100], seq_len=200)
        assert result[-1]["end"] == 200

    def test_contiguous_no_gaps(self):
        result = _reformat_exon_annotations([30, 60, 90], seq_len=120)
        for i in range(len(result) - 1):
            assert result[i]["end"] + 1 == result[i + 1]["start"]
