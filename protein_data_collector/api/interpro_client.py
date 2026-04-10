"""Synchronous InterPro REST API client."""

import logging
import time
from typing import Any, Dict, List, Optional, Set

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

    def get_tim_barrel_pfam_entries(self) -> List[Dict[str, Any]]:
        """Return all PFAM entries annotated as TIM barrel."""
        return self._paginate("entry/pfam/", params={"annotation": "TIM barrel"})

    def get_tim_barrel_interpro_entries(self) -> List[Dict[str, Any]]:
        """Return all InterPro entries annotated as TIM barrel."""
        return self._paginate("entry/interpro/", params={"annotation": "TIM barrel"})

    def get_human_proteins_for_entry(self, accession: str) -> List[str]:
        """Return UniProt IDs of human proteins belonging to *accession*."""
        endpoint = f"protein/uniprot/entry/interpro/{accession}/"
        results = self._paginate(endpoint, params={"taxonomy_id": "9606"})
        return [r.get("metadata", {}).get("accession") for r in results
                if r.get("metadata", {}).get("accession")]

    def get_domain_boundaries(
        self, uniprot_id: str, tim_barrel_accessions: Set[str]
    ) -> Optional[Dict[str, Any]]:
        """
        Return the first TIM barrel domain boundary for *uniprot_id* from InterPro.

        Matches entries whose accession is in *tim_barrel_accessions* (the set
        collected during Phase 1).  Returns a dict {domain_id, start, end,
        length, source} or None.
        """
        data = self._get(f"protein/uniprot/{uniprot_id}")
        if not data or "results" not in data:
            return None

        for result in data["results"]:
            for entry in result.get("entries", []):
                meta = entry.get("metadata", {})
                if meta.get("accession") in tim_barrel_accessions:
                    for loc in entry.get("entry_protein_locations", []):
                        frags = loc.get("fragments", [])
                        if frags:
                            start = frags[0].get("start")
                            end = frags[-1].get("end")
                            if start and end:
                                return {
                                    "domain_id": meta.get("accession"),
                                    "domain_name": meta.get("name"),
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
