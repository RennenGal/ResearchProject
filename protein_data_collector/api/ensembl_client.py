"""
Ensembl REST API client.

Covers the three operations needed for transcript expansion:
  1. ensg_for_enst(enst_id)          ENST → ENSG (gene lookup)
  2. ensg_for_uniprot(uniprot_id)     UniProt accession → ENSG (xref lookup)
  3. transcripts_for_gene(ensg_id)    ENSG → list of protein-coding transcripts
  4. protein_sequence(enst_id)        ENST → amino-acid sequence

All IDs are returned without version suffixes (ENST… not ENST….4).
Rate limit: 15 req/s on the public HTTPS endpoint; the client sleeps 0.07 s
between calls automatically.
"""

import logging
import time
from typing import Optional

import requests

logger = logging.getLogger(__name__)

_BASE = "https://rest.ensembl.org"
_HEADERS = {"Content-Type": "application/json"}
_DELAY = 0.07   # ~14 req/s, safely under the 15 req/s limit


_MAX_RETRIES = 5
_RETRY_BACKOFF = [5, 15, 30, 60, 120]   # seconds between retries


def _get(path: str, params: Optional[dict] = None, timeout: int = 30) -> Optional[dict | list]:
    url = f"{_BASE}{path}"
    for attempt in range(_MAX_RETRIES):
        try:
            r = requests.get(url, headers=_HEADERS, params=params, timeout=timeout)
            time.sleep(_DELAY)
            if r.status_code == 200:
                return r.json()
            if r.status_code == 429:
                retry = int(r.headers.get("Retry-After", 10))
                logger.warning("Rate limited — sleeping %d s", retry)
                time.sleep(retry)
                continue
            logger.debug("Ensembl %s → HTTP %d", path, r.status_code)
            return None
        except requests.RequestException as e:
            wait = _RETRY_BACKOFF[min(attempt, len(_RETRY_BACKOFF) - 1)]
            logger.warning("Ensembl request failed for %s (attempt %d/%d): %s — retrying in %ds",
                           path, attempt + 1, _MAX_RETRIES, e, wait)
            time.sleep(wait)
    logger.error("Ensembl request permanently failed for %s after %d attempts", path, _MAX_RETRIES)
    return None


def _strip_version(eid: Optional[str]) -> Optional[str]:
    if not eid:
        return None
    return eid.split(".")[0]


# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------

def ensg_for_enst(enst_id: str) -> Optional[str]:
    """Return the ENSG gene ID for a given ENST transcript ID."""
    enst = _strip_version(enst_id)
    data = _get(f"/lookup/id/{enst}", params={"content-type": "application/json"})
    if isinstance(data, dict):
        return _strip_version(data.get("Parent"))
    return None


def ensg_for_uniprot(uniprot_id: str, species: str = "homo_sapiens") -> Optional[str]:
    """Return the first ENSG gene ID found for a UniProt accession via xrefs."""
    data = _get(
        f"/xrefs/symbol/{species}/{uniprot_id}",
        params={"object_type": "gene", "content-type": "application/json"},
    )
    if isinstance(data, list) and data:
        return _strip_version(data[0].get("id"))
    # fallback: direct xref lookup
    data2 = _get(
        f"/xrefs/id/{uniprot_id}",
        params={"content-type": "application/json", "external_db": "UniProtKB/Swiss-Prot"},
    )
    if isinstance(data2, list):
        for ref in data2:
            if ref.get("type") == "gene":
                return _strip_version(ref.get("id"))
    return None


def transcripts_for_gene(ensg_id: str) -> list[dict]:
    """
    Return all protein-coding transcripts for a gene.

    Each dict has:
        enst_id, ensp_id, biotype, is_mane_select, length
    """
    ensg = _strip_version(ensg_id)
    data = _get(
        f"/lookup/id/{ensg}",
        params={"expand": 1, "content-type": "application/json"},
    )
    if not isinstance(data, dict):
        return []

    results = []
    for tx in data.get("Transcript", []):
        if tx.get("biotype") != "protein_coding":
            continue
        enst = _strip_version(tx.get("id"))
        if not enst:
            continue
        translation = tx.get("Translation")
        ensp = _strip_version(translation.get("id")) if isinstance(translation, dict) else None
        is_canonical = 1 if tx.get("is_canonical") else 0
        results.append({
            "enst_id":       enst,
            "ensp_id":       ensp,
            "biotype":       tx.get("biotype"),
            "is_mane_select": is_canonical,
            "length":        tx.get("length", 0),
        })
    return results


def protein_sequence(enst_id: str) -> Optional[str]:
    """Return the translated amino-acid sequence for a transcript, or None."""
    enst = _strip_version(enst_id)
    data = _get(
        f"/sequence/id/{enst}",
        params={"type": "protein", "content-type": "application/json"},
    )
    if isinstance(data, dict):
        seq = data.get("seq", "")
        # Ensembl sometimes returns sequences ending with * (stop codon)
        return seq.rstrip("*") if seq else None
    return None


def transcript_exon_boundaries(enst_id: str) -> list[int]:
    """
    Return protein-space exon boundary positions for a transcript.

    Each integer is the 1-based amino-acid position of the last residue
    contributed by that exon (ceiling division of cumulative CDS bases).
    The final exon is excluded — its downstream junction is the end of the protein.

    Returns [] if data is unavailable or the transcript has no Translation.

    Algorithm
    ---------
    1. Fetch /lookup/id/{ENST}?expand=1 to get Exon list + Translation CDS bounds.
    2. Sort exons by rank (transcription order, works for both strands because
       Ensembl rank reflects transcription direction, not genomic direction).
    3. Clip each exon to [Translation.start, Translation.end] to skip UTR bases.
    4. Accumulate CDS bases; after each exon (except the last) compute
       boundary = ceil(cumulative / 3).
    """
    enst = _strip_version(enst_id)
    data = _get(f"/lookup/id/{enst}", params={"expand": 1, "content-type": "application/json"})
    if not isinstance(data, dict):
        return []

    translation = data.get("Translation")
    if not isinstance(translation, dict):
        return []

    cds_start = translation.get("start")
    cds_end   = translation.get("end")
    if cds_start is None or cds_end is None:
        return []

    exons = data.get("Exon", [])
    if not exons:
        return []

    exons_sorted = sorted(exons, key=lambda e: e.get("rank", 0))

    # Filter to only coding exons (overlapping the CDS)
    coding_exons = []
    for ex in exons_sorted:
        ex_start = ex.get("start")
        ex_end   = ex.get("end")
        if ex_start is None or ex_end is None:
            continue
        clipped_start = max(ex_start, cds_start)
        clipped_end   = min(ex_end,   cds_end)
        if clipped_end < clipped_start:
            continue
        coding_exons.append(clipped_end - clipped_start + 1)

    boundaries = []
    cumulative = 0
    for i, cds_bases in enumerate(coding_exons):
        cumulative += cds_bases
        if i < len(coding_exons) - 1:
            # Ceiling division: last AA that has any bases in this exon
            boundaries.append((cumulative + 2) // 3)

    return boundaries
