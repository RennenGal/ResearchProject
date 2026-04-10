"""Collect TIM barrel families and their human proteins from InterPro."""

import logging
from typing import List, Optional

from ..api.interpro_client import InterProClient
from ..models.entities import Protein, TIMBarrelEntry

logger = logging.getLogger(__name__)


class InterProCollector:
    def __init__(self, client: Optional[InterProClient] = None):
        self.client = client or InterProClient()

    # ------------------------------------------------------------------
    # TIM barrel entries (Phase 1)
    # ------------------------------------------------------------------

    def collect_tim_barrel_entries(self) -> List[TIMBarrelEntry]:
        """Return all PFAM and InterPro entries annotated as TIM barrel."""
        entries: List[TIMBarrelEntry] = []

        logger.info("Fetching PFAM TIM barrel entries from InterPro...")
        for raw in self.client.get_tim_barrel_pfam_entries():
            entry = _parse_tim_barrel_entry(raw, entry_type="pfam")
            if entry:
                entries.append(entry)

        logger.info("Fetching InterPro TIM barrel entries...")
        for raw in self.client.get_tim_barrel_interpro_entries():
            entry = _parse_tim_barrel_entry(raw, entry_type="interpro")
            if entry:
                entries.append(entry)

        logger.info("Collected %d TIM barrel entries", len(entries))
        return entries

    # ------------------------------------------------------------------
    # Human proteins (Phase 2)
    # ------------------------------------------------------------------

    def collect_human_proteins(
        self, entries: List[TIMBarrelEntry]
    ) -> List[Protein]:
        """Return human proteins for every entry in *entries*, deduplicated."""
        seen: set = set()
        proteins: List[Protein] = []

        for entry in entries:
            logger.info("Fetching proteins for %s (%s)", entry.accession, entry.name)
            try:
                uniprot_ids = self.client.get_human_proteins_for_entry(entry.accession)
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
                        organism="Homo sapiens",
                    )
                )

        logger.info("Collected %d unique human proteins", len(proteins))
        return proteins


# ---------------------------------------------------------------------------
# Pure parsing helpers
# ---------------------------------------------------------------------------

def _parse_tim_barrel_entry(raw: dict, entry_type: str) -> Optional[TIMBarrelEntry]:
    meta = raw.get("metadata", {})
    accession = meta.get("accession", "")
    if not accession:
        return None

    # The annotation field differs between PFAM and InterPro responses
    annotation = (
        meta.get("integrated", "")
        or meta.get("description", "")
        or "TIM barrel"
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
