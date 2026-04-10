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
