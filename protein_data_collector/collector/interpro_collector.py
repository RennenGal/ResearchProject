"""Collect domain family entries and their proteins from InterPro."""

import logging
from typing import List, Optional

from ..api.interpro_client import InterProClient
from ..models.entities import Protein, TIMBarrelEntry

logger = logging.getLogger(__name__)


class InterProCollector:
    def __init__(self, client: Optional[InterProClient] = None):
        self.client = client or InterProClient()

    # ------------------------------------------------------------------
    # Domain entries (Phase 1)
    # ------------------------------------------------------------------

    def collect_domain_entries(
        self,
        annotation: str,
        search: str = "",
        extra_accessions: tuple = (),
    ) -> List[TIMBarrelEntry]:
        """
        Return all PFAM and InterPro entries for a domain using up to three strategies:
        1. annotation= query (exact InterPro annotation term, e.g. 'TIM barrel')
        2. search= query (text search, e.g. 'propeller')
        3. extra_accessions — mandatory accessions to always include

        Results are deduplicated by accession.
        """
        seen: set = set()
        entries: List[TIMBarrelEntry] = []

        def _add(raw: dict, entry_type: str, fallback: str) -> None:
            entry = _parse_domain_entry(raw, entry_type=entry_type, fallback_annotation=fallback)
            if entry and entry.accession not in seen:
                seen.add(entry.accession)
                entries.append(entry)

        if annotation:
            logger.info("Fetching PFAM entries via annotation='%s'...", annotation)
            for raw in self.client.get_domain_pfam_entries(annotation):
                _add(raw, "pfam", annotation)
            logger.info("Fetching InterPro entries via annotation='%s'...", annotation)
            for raw in self.client.get_domain_interpro_entries(annotation):
                _add(raw, "interpro", annotation)

        if search:
            logger.info("Fetching PFAM entries via search='%s'...", search)
            for raw in self.client.search_pfam_entries(search):
                _add(raw, "pfam", search)
            logger.info("Fetching InterPro entries via search='%s'...", search)
            for raw in self.client.search_interpro_entries(search):
                _add(raw, "interpro", search)

        for accession in extra_accessions:
            if accession in seen:
                continue
            logger.info("Fetching extra entry %s...", accession)
            db_type = "pfam" if accession.startswith("PF") else "interpro"
            raw = self.client.get_entry(accession)
            if raw:
                _add({"metadata": raw.get("metadata", raw)}, db_type, accession)

        logger.info("Collected %d domain entries total", len(entries))
        return entries

    # ------------------------------------------------------------------
    # Proteins (Phase 2)
    # ------------------------------------------------------------------

    def collect_proteins(
        self,
        entries: List[TIMBarrelEntry],
        organism: str,
        taxon_id: int,
    ) -> List[Protein]:
        """Return proteins for *organism* (by *taxon_id*) across all *entries*, deduplicated."""
        seen: set = set()
        proteins: List[Protein] = []

        for entry in entries:
            logger.info("Fetching %s proteins for %s (%s)", organism, entry.accession, entry.name)
            try:
                uniprot_ids = self.client.get_proteins_for_entry(entry.accession, taxon_id)
            except Exception as e:
                logger.error("Failed to fetch proteins for %s: %s", entry.accession, e)
                continue

            for uid in uniprot_ids:
                if not uid or uid in seen:
                    continue
                seen.add(uid)
                proteins.append(
                    Protein(
                        uniprot_id=uid,
                        tim_barrel_accession=entry.accession,
                        organism=organism,
                    )
                )

        logger.info("Collected %d unique %s proteins", len(proteins), organism)
        return proteins


# ---------------------------------------------------------------------------
# Pure parsing helpers
# ---------------------------------------------------------------------------

def _parse_domain_entry(
    raw: dict,
    entry_type: str,
    fallback_annotation: str = "domain",
) -> Optional[TIMBarrelEntry]:
    meta = raw.get("metadata", {})
    accession = meta.get("accession", "")
    if not accession:
        return None

    annotation = (
        meta.get("integrated", "")
        or meta.get("description", "")
        or fallback_annotation
    )

    try:
        return TIMBarrelEntry(
            accession=accession,
            entry_type=entry_type,
            name=meta.get("name", accession),
            description=meta.get("description"),
            tim_barrel_annotation=annotation,
        )
    except Exception as e:
        logger.warning("Skipping malformed entry %s: %s", accession, e)
        return None
