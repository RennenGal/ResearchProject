"""Synchronous UniProt REST API client."""

import logging
import time
from typing import Any, Dict, List, Optional

import requests

from ..config import get_config
from ..errors import APIError, NetworkError
from ..retry import with_retry

logger = logging.getLogger(__name__)

_UNIPROT_FIELDS = ",".join([
    "accession", "id", "protein_name", "gene_names", "organism_name",
    "sequence", "length", "reviewed", "protein_existence", "annotation_score",
    "cc_alternative_products", "ft_var_seq",
    "xref_interpro", "xref_ensembl", "xref_alphafolddb",
])

_GENE_NAME_BATCH = 100  # UniProt limits OR conditions to 100 per search query
_METADATA_BATCH  = 100


class UniProtClient:
    """Fetch protein and isoform data from the UniProt REST API."""

    def __init__(self):
        self.cfg = get_config()
        self.session = requests.Session()
        self.session.headers.update({"Accept": "application/json"})

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def get_protein(self, uniprot_id: str) -> Optional[Dict[str, Any]]:
        """Return the full UniProt JSON entry for *uniprot_id*, or None if not found."""
        return self._fetch(uniprot_id)

    def get_isoform_sequence(self, isoform_id: str) -> Optional[str]:
        """Return the amino acid sequence for *isoform_id* (e.g. 'P04637-2'), or None."""
        data = self._fetch(isoform_id)
        if data:
            return data.get("sequence", {}).get("value")
        return None

    def batch_gene_names(self, uniprot_ids: List[str]) -> Dict[str, Optional[str]]:
        """
        Fetch the primary gene name for each UniProt accession in *uniprot_ids*.

        Returns {uniprot_id: gene_name_or_None}.  Proteins not found in UniProt
        are mapped to None.  Requests are chunked to stay within URL length limits.
        """
        result: Dict[str, Optional[str]] = {uid: None for uid in uniprot_ids}
        search_url = f"{self.cfg.uniprot_base_url}/search"

        for i in range(0, len(uniprot_ids), _GENE_NAME_BATCH):
            chunk = uniprot_ids[i:i + _GENE_NAME_BATCH]
            query = " OR ".join(f"accession:{uid}" for uid in chunk)
            params = {
                "query": query,
                "fields": "accession,gene_names",
                "format": "json",
                "size": len(chunk),
            }
            try:
                resp = self.session.get(search_url, params=params, timeout=self.cfg.request_timeout)
                time.sleep(self.cfg.request_delay)
            except requests.exceptions.RequestException as e:
                logger.warning("Network error fetching gene names (chunk %d): %s", i, e)
                continue

            if not resp.ok:
                logger.warning("UniProt search returned %d for gene name chunk %d", resp.status_code, i)
                continue

            for entry in resp.json().get("results", []):
                uid = entry.get("primaryAccession")
                if uid and uid in result:
                    genes = entry.get("genes", [])
                    if genes and genes[0].get("geneName"):
                        result[uid] = genes[0]["geneName"]["value"]

            logger.debug("Gene name batch %d/%d done", min(i + _GENE_NAME_BATCH, len(uniprot_ids)), len(uniprot_ids))

        return result

    def batch_protein_metadata(
        self, uniprot_ids: List[str]
    ) -> Dict[str, Dict[str, Any]]:
        """
        Fetch protein_name, reviewed, and annotation_score for each accession.

        Returns {uniprot_id: {"protein_name": str|None, "reviewed": int, "annotation_score": int|None}}.
        Missing accessions are mapped to {"protein_name": None, "reviewed": 0, "annotation_score": None}.
        """
        default: Dict[str, Any] = {"protein_name": None, "reviewed": 0, "annotation_score": None}
        result: Dict[str, Dict[str, Any]] = {uid: dict(default) for uid in uniprot_ids}
        search_url = f"{self.cfg.uniprot_base_url}/search"

        for i in range(0, len(uniprot_ids), _METADATA_BATCH):
            chunk = uniprot_ids[i:i + _METADATA_BATCH]
            query = " OR ".join(f"accession:{uid}" for uid in chunk)
            params = {
                "query":  query,
                "fields": "accession,protein_name,reviewed,annotation_score",
                "format": "json",
                "size":   len(chunk),
            }
            try:
                resp = self.session.get(search_url, params=params, timeout=self.cfg.request_timeout)
                time.sleep(self.cfg.request_delay)
            except requests.exceptions.RequestException as e:
                logger.warning("Network error fetching metadata (chunk %d): %s", i, e)
                continue

            if not resp.ok:
                logger.warning("UniProt search returned %d for metadata chunk %d", resp.status_code, i)
                continue

            for entry in resp.json().get("results", []):
                uid = entry.get("primaryAccession")
                if not uid or uid not in result:
                    continue

                # protein_name: prefer recommendedName, fall back to submissionNames (TrEMBL)
                desc = entry.get("proteinDescription", {})
                rec  = desc.get("recommendedName", {})
                name: Optional[str] = None
                if rec:
                    name = rec.get("fullName", {}).get("value")
                if not name:
                    submitted = desc.get("submissionNames", [])
                    if submitted:
                        name = submitted[0].get("fullName", {}).get("value")

                # reviewed: only "UniProtKB reviewed (Swiss-Prot)" → 1; "unreviewed" must NOT match
                entry_type = entry.get("entryType", "")
                reviewed = 1 if "swiss-prot" in entry_type.lower() else 0

                # annotation_score: integer 1-5
                score_raw = entry.get("annotationScore")
                ann_score: Optional[int] = int(score_raw) if score_raw is not None else None

                result[uid] = {
                    "protein_name":    name,
                    "reviewed":        reviewed,
                    "annotation_score": ann_score,
                }

            logger.debug(
                "Metadata batch %d/%d done",
                min(i + _METADATA_BATCH, len(uniprot_ids)),
                len(uniprot_ids),
            )

        return result

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    @with_retry()
    def _fetch(self, entry_id: str) -> Optional[Dict[str, Any]]:
        url = f"{self.cfg.uniprot_base_url}/{entry_id}"
        try:
            resp = self.session.get(url, params={"format": "json"}, timeout=self.cfg.request_timeout)
            time.sleep(self.cfg.request_delay)
        except requests.exceptions.ConnectionError as e:
            raise NetworkError(f"Connection error fetching {entry_id}: {e}") from e
        except requests.exceptions.Timeout as e:
            raise NetworkError(f"Timeout fetching {entry_id}: {e}") from e

        if resp.status_code == 404:
            logger.warning("UniProt entry not found: %s", entry_id)
            return None
        if resp.status_code == 429:
            raise APIError("UniProt rate limit exceeded", status_code=429)
        if not resp.ok:
            raise APIError(
                f"UniProt returned {resp.status_code} for {entry_id}",
                status_code=resp.status_code,
            )
        return resp.json()
