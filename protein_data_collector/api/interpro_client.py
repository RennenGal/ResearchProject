"""Synchronous InterPro REST API client."""

import logging
import time
from typing import Any, Dict, List, Optional

import requests

from ..config import get_config
from ..errors import APIError, NetworkError
from ..retry import with_retry

logger = logging.getLogger(__name__)


class InterProClient:
    """Fetch TIM barrel family and protein data from the InterPro REST API."""

    def __init__(self):
        self.cfg = get_config()
        self.session = requests.Session()
        self.session.headers.update({"Accept": "application/json"})

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def get_domain_pfam_entries(self, annotation: str) -> List[Dict[str, Any]]:
        """Return all PFAM entries matching *annotation* (e.g. 'TIM barrel')."""
        return self._paginate("entry/pfam/", params={"annotation": annotation})

    def get_domain_interpro_entries(self, annotation: str) -> List[Dict[str, Any]]:
        """Return all InterPro entries matching *annotation*."""
        return self._paginate("entry/interpro/", params={"annotation": annotation})

    def search_pfam_entries(self, search: str) -> List[Dict[str, Any]]:
        """Return all PFAM entries whose name/description contains *search*."""
        return self._paginate("entry/pfam/", params={"search": search})

    def search_interpro_entries(self, search: str) -> List[Dict[str, Any]]:
        """Return all InterPro entries whose name/description contains *search*."""
        return self._paginate("entry/interpro/", params={"search": search})

    def get_entry(self, accession: str) -> Optional[Dict[str, Any]]:
        """Fetch a single entry by accession (pfam or interpro)."""
        db = "pfam" if accession.startswith("PF") else "interpro"
        return self._get(f"entry/{db}/{accession}")

    def get_proteins_for_entry(self, accession: str, taxon_id: int) -> List[str]:
        """Return UniProt IDs of proteins for *taxon_id* belonging to *accession*."""
        db = "pfam" if accession.startswith("PF") else "interpro"
        endpoint = f"protein/uniprot/taxonomy/uniprot/{taxon_id}/entry/{db}/{accession}/"
        results = self._paginate(endpoint)
        return [r.get("metadata", {}).get("accession") for r in results
                if r.get("metadata", {}).get("accession")]

    def get_human_proteins_for_entry(self, accession: str) -> List[str]:
        """Return UniProt IDs of human (taxon 9606) proteins for *accession*."""
        return self.get_proteins_for_entry(accession, taxon_id=9606)

    def get_domain_boundaries(
        self, uniprot_id: str, tim_barrel_accession: str
    ) -> Optional[Dict[str, Any]]:
        """
        Return the TIM barrel domain boundary for *uniprot_id* from InterPro.

        Uses the entry-centric endpoint:
            entry/{db}/{accession}/protein/uniprot/{uid}
        which directly returns the protein's entry_protein_locations for that entry.

        Returns a dict {domain_id, start, end, length, source} or None.
        """
        db = "pfam" if tim_barrel_accession.startswith("PF") else "interpro"
        endpoint = f"entry/{db}/{tim_barrel_accession}/protein/uniprot/{uniprot_id}"
        data = self._get(endpoint)
        if not data:
            return None

        for protein in data.get("proteins", []):
            for loc in protein.get("entry_protein_locations", []):
                frags = loc.get("fragments", [])
                if frags:
                    start = frags[0].get("start")
                    end = frags[-1].get("end")
                    if start and end:
                        return {
                            "domain_id": tim_barrel_accession,
                            "start": start,
                            "end": end,
                            "length": end - start + 1,
                            "source": "interpro_api",
                        }
        return None

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    @with_retry()
    def _get(self, endpoint: str, params: Dict[str, Any] = None) -> Optional[Dict[str, Any]]:
        url = f"{self.cfg.interpro_base_url}/{endpoint.lstrip('/')}"
        try:
            resp = self.session.get(url, params=params or {}, timeout=self.cfg.request_timeout)
            time.sleep(self.cfg.request_delay)
        except requests.exceptions.ConnectionError as e:
            raise NetworkError(f"Connection error: {e}") from e
        except requests.exceptions.Timeout as e:
            raise NetworkError(f"Timeout: {e}") from e

        if resp.status_code == 404:
            return None
        if resp.status_code == 429:
            raise APIError("InterPro rate limit exceeded", status_code=429)
        if not resp.ok:
            raise APIError(f"InterPro returned {resp.status_code}", status_code=resp.status_code)
        if not resp.content:
            return None
        return resp.json()

    def _paginate(self, endpoint: str, params: Dict[str, Any] = None) -> List[Dict[str, Any]]:
        """Collect all pages from a paginated InterPro endpoint."""
        all_results: List[Dict[str, Any]] = []
        next_url = f"{self.cfg.interpro_base_url}/{endpoint.lstrip('/')}"
        p = params or {}

        while next_url:
            try:
                resp = self.session.get(next_url, params=p, timeout=self.cfg.request_timeout)
                time.sleep(self.cfg.request_delay)
            except requests.exceptions.RequestException as e:
                raise NetworkError(f"Pagination error at {next_url}: {e}") from e

            if not resp.ok:
                raise APIError(f"InterPro returned {resp.status_code}", status_code=resp.status_code)

            data = resp.json()
            all_results.extend(data.get("results", []))
            next_url = data.get("next")
            p = {}  # params are encoded in next_url after first page

        return all_results
