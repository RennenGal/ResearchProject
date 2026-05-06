"""
TIM barrel (beta-alpha)_8 motif identification via DSSP secondary structure.

Given a numpy array of per-residue SS assignments (H / E / -) from pydssp,
and the 1-based [domain_start, domain_end] boundaries, identifies the 8
beta-alpha repeat units that define the TIM barrel.

Strategy
--------
1. Slice the SS array to the domain region.
2. Find all contiguous E-runs (beta-strands, min 3 residues) and H-runs
   (alpha-helices, min 4 residues) within the domain.
3. Merge strands or helices separated by <= 2 loop residues (handles DSSP
   fragmentation of long secondary structure elements).
4. Pair each beta-strand with the immediately following alpha-helix.
5. Return the first 8 pairs as motifs (or fewer if the domain is incomplete).

All positions are 1-based, inclusive, in full-protein coordinates.
"""

from __future__ import annotations

import numpy as np

MIN_STRAND  = 3   # minimum residues for a beta-strand to be counted
MIN_HELIX   = 4   # minimum residues for an alpha-helix to be counted
MERGE_GAP   = 2   # loop gap <= this many residues merges adjacent SS elements


def _runs(ss_slice: np.ndarray, char: str, domain_offset: int) -> list[tuple[int, int]]:
    """
    Return list of (start, end) 1-based full-protein positions for all
    contiguous runs of *char* in *ss_slice*.  Runs separated by at most
    MERGE_GAP other characters are merged before filtering by minimum length.
    """
    n = len(ss_slice)
    raw: list[tuple[int, int]] = []
    i = 0
    while i < n:
        if ss_slice[i] == char:
            j = i + 1
            while j < n and ss_slice[j] == char:
                j += 1
            raw.append((i, j - 1))  # 0-based in domain space
            i = j
        else:
            i += 1

    # Merge segments separated by <= MERGE_GAP positions
    merged: list[tuple[int, int]] = []
    for seg in raw:
        if merged and seg[0] - merged[-1][1] - 1 <= MERGE_GAP:
            merged[-1] = (merged[-1][0], seg[1])
        else:
            merged.append(list(seg))

    min_len = MIN_STRAND if char == 'E' else MIN_HELIX
    result = []
    for s, e in merged:
        if e - s + 1 >= min_len:
            # Convert to 1-based full-protein coordinates
            result.append((domain_offset + s + 1, domain_offset + e + 1))
    return result


def identify_ba_motifs(
    ss: np.ndarray,
    domain_start: int,
    domain_end: int,
) -> list[dict]:
    """
    Identify up to 8 beta-alpha motifs in *ss* within [domain_start, domain_end].

    Parameters
    ----------
    ss          : 1-D array of characters ('H', 'E', '-'), length = protein length.
    domain_start: 1-based start of TIM barrel domain in full protein.
    domain_end  : 1-based end   of TIM barrel domain in full protein.

    Returns
    -------
    List of up to 8 dicts:
        {motif, start, end, beta_start, beta_end, alpha_start, alpha_end}
    All positions 1-based, inclusive, in full-protein coordinates.
    """
    # Slice to domain (0-based Python indexing)
    dom_slice = ss[domain_start - 1 : domain_end]
    offset    = domain_start - 1  # add to 0-based domain index to get 0-based protein index

    strands = _runs(dom_slice, 'E', offset)
    helices = _runs(dom_slice, 'H', offset)

    motifs: list[dict] = []
    helix_idx = 0

    for beta_start, beta_end in strands:
        # Find the first helix that starts after this strand ends
        while helix_idx < len(helices) and helices[helix_idx][0] <= beta_end:
            helix_idx += 1

        if helix_idx >= len(helices):
            break

        alpha_start, alpha_end = helices[helix_idx]
        motifs.append({
            "motif":       len(motifs) + 1,
            "start":       beta_start,
            "end":         alpha_end,
            "beta_start":  beta_start,
            "beta_end":    beta_end,
            "alpha_start": alpha_start,
            "alpha_end":   alpha_end,
        })
        helix_idx += 1

        if len(motifs) == 8:
            break

    return motifs
