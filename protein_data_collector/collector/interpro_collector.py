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
        cathgene3d_search: str = "",
        extra_accessions: tuple = (),
    ) -> List[TIMBarrelEntry]:
        """
        Return all domain entries using up to four strategies:
        1. annotation= query  (exact InterPro annotation term, e.g. 'TIM barrel')
        2. search= query      (text search on pfam+interpro, e.g. 'propeller')
        3. cathgene3d_search  (search CATH Gene3D for structurally-classified entries
                               that have no Pfam/InterPro parent, e.g. '3.20.20')
        4. extra_accessions   (mandatory accessions; fetched individually, never skipped)

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

        if cathgene3d_search:
            logger.info("Fetching CATH Gene3D entries via search='%s'...", cathgene3d_search)
            for raw in self.client.search_cathgene3d_entries(cathgene3d_search):
                acc = raw.get("metadata", {}).get("accession", "")
                # Only keep entries whose accession matches the search prefix
                # (e.g. '3.20.20' filters to TIM barrel superfamilies only)
                if cathgene3d_search in acc:
                    _add(raw, "cathgene3d", cathgene3d_search)

        for accession in extra_accessions:
            if accession in seen:
                continue
            logger.info("Fetching extra entry %s...", accession)
            from ..api.interpro_client import _db_for_accession
            db_type = _db_for_accession(accession)
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

    # `name` can be a plain string (paginated results) or a dict with 'name'/'short'
    # keys (individual GET responses).
    raw_name = meta.get("name", accession)
    if isinstance(raw_name, dict):
        name = raw_name.get("name") or raw_name.get("short") or accession
    else:
        name = raw_name or accession

    # `description` can be a plain string or a list of rich-text dicts from the API.
    raw_desc = meta.get("description")
    if isinstance(raw_desc, list):
        description = " ".join(
            d.get("text", "") for d in raw_desc if isinstance(d, dict)
        ).strip() or None
    else:
        description = raw_desc or None

    # `integrated` is the parent IPR accession for PFAM/Gene3D entries.
    integrated = meta.get("integrated")
    if isinstance(integrated, dict):
        integrated = integrated.get("accession", "")

    annotation = integrated or description or fallback_annotation

    try:
        return TIMBarrelEntry(
            accession=accession,
            entry_type=entry_type,
            name=name,
            description=description,
            domain_annotation=annotation,
        )
    except Exception as e:
        logger.warning("Skipping malformed entry %s: %s", accession, e)
        return None
