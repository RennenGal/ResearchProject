"""
Orchestrate the three-phase collection pipeline.

Phase 1 — TIM barrel families  : InterPro → tim_barrel_entries
Phase 2 — Human proteins       : InterPro → proteins
Phase 3 — Isoforms             : UniProt  → isoforms
"""

import json
import logging
from dataclasses import dataclass, field
from typing import List, Optional

from ..database.connection import ensure_db, get_connection
from ..database.storage import (
    deduplicate_proteins,
    get_all_proteins,
    get_proteins_without_isoforms,
    get_counts,
    upsert_isoforms,
    upsert_proteins,
    upsert_tim_barrel_entries,
)
from ..models.entities import Isoform, Protein, TIMBarrelEntry
from .interpro_collector import InterProCollector
from .uniprot_collector import UniProtCollector

logger = logging.getLogger(__name__)


@dataclass
class CollectionReport:
    tim_barrel_entries: int = 0
    proteins_collected: int = 0
    isoforms_collected: int = 0
    failed_proteins: List[str] = field(default_factory=list)

    @property
    def alternative_isoforms(self) -> int:
        return self.isoforms_collected - self.proteins_collected

    def summary(self) -> str:
        return (
            f"TIM barrel entries : {self.tim_barrel_entries}\n"
            f"Proteins           : {self.proteins_collected}\n"
            f"Isoforms           : {self.isoforms_collected}\n"
            f"  (canonical       : {self.proteins_collected})\n"
            f"  (alternative     : {self.alternative_isoforms})\n"
            f"Failed proteins    : {len(self.failed_proteins)}"
        )


class DataCollector:
    def __init__(
        self,
        interpro_collector: Optional[InterProCollector] = None,
        uniprot_collector: Optional[UniProtCollector] = None,
        db_path: Optional[str] = None,
    ):
        self.interpro = interpro_collector or InterProCollector()
        self.uniprot = uniprot_collector or UniProtCollector()
        self.db_path = db_path

    # ------------------------------------------------------------------
    # Full pipeline
    # ------------------------------------------------------------------

    def run_full_collection(self) -> CollectionReport:
        """Run all three phases from scratch."""
        ensure_db(self.db_path)
        report = CollectionReport()

        entries = self._phase1_tim_barrel_entries()
        report.tim_barrel_entries = len(entries)

        proteins = self._phase2_human_proteins(entries)
        proteins = self._phase2b_deduplicate(proteins)
        report.proteins_collected = len(proteins)

        isoforms = self._phase3_isoforms(proteins, report)
        report.isoforms_collected = len(isoforms)

        logger.info("Collection complete.\n%s", report.summary())
        return report

    def recollect_all_isoforms(self) -> CollectionReport:
        """Delete all isoform rows and re-fetch from UniProt (gets alternative isoforms)."""
        ensure_db(self.db_path)
        with get_connection(self.db_path) as conn:
            conn.execute("DELETE FROM isoforms")
            conn.commit()
            logger.info("Cleared isoforms table; re-collecting from UniProt")
            proteins_raw = get_all_proteins(conn)

        proteins = [
            Protein(**p) for p in proteins_raw
            if p.get("canonical_uniprot_id") is None
        ]
        report = CollectionReport()
        isoforms = self._phase3_isoforms(proteins, report)
        report.isoforms_collected = len(isoforms)
        logger.info("Re-collection complete.\n%s", report.summary())
        return report

    def resume_isoform_collection(self) -> CollectionReport:
        """Phase 3 only — collect isoforms for proteins not yet in the database."""
        ensure_db(self.db_path)
        report = CollectionReport()

        with get_connection(self.db_path) as conn:
            remaining_ids = get_proteins_without_isoforms(conn)
            all_proteins = {p["uniprot_id"]: p for p in get_all_proteins(conn)}

        proteins = [
            Protein(**all_proteins[uid])
            for uid in remaining_ids
            if uid in all_proteins
        ]
        logger.info("%d proteins still need isoform collection", len(proteins))

        isoforms = self._phase3_isoforms(proteins, report)
        report.isoforms_collected = len(isoforms)
        logger.info("Resume complete.\n%s", report.summary())
        return report

    def backfill_domain_locations(self) -> int:
        """
        Fetch and store tim_barrel_location for canonical isoforms that currently have NULL.

        Returns the number of isoforms updated.
        """
        ensure_db(self.db_path)

        with get_connection(self.db_path) as conn:
            rows = conn.execute(
                """
                SELECT i.isoform_id, i.uniprot_id, p.tim_barrel_accession
                FROM isoforms i
                JOIN proteins p ON i.uniprot_id = p.uniprot_id
                WHERE i.is_canonical = 1 AND i.tim_barrel_location IS NULL
                """
            ).fetchall()

        updated = 0
        total = len(rows)
        logger.info("Backfilling domain locations for %d canonical isoforms", total)

        for idx, row in enumerate(rows, 1):
            loc = self.uniprot._get_tim_barrel_location(
                row["uniprot_id"], row["tim_barrel_accession"]
            )
            if loc:
                with get_connection(self.db_path) as conn:
                    conn.execute(
                        "UPDATE isoforms SET tim_barrel_location = ? WHERE isoform_id = ?",
                        (json.dumps(loc), row["isoform_id"]),
                    )
                    conn.commit()
                updated += 1
            if idx % 25 == 0 or idx == total:
                logger.info("  %d/%d — %d updated", idx, total, updated)

        logger.info("Backfill complete: %d/%d isoforms have domain location", updated, total)
        return updated

    # ------------------------------------------------------------------
    # Phases
    # ------------------------------------------------------------------

    def _phase1_tim_barrel_entries(self) -> List[TIMBarrelEntry]:
        logger.info("Phase 1: collecting TIM barrel entries from InterPro")
        entries = self.interpro.collect_tim_barrel_entries()
        with get_connection(self.db_path) as conn:
            upsert_tim_barrel_entries(conn, entries)
        return entries

    def _phase2_human_proteins(self, entries: List[TIMBarrelEntry]) -> List[Protein]:
        logger.info("Phase 2: collecting human proteins from InterPro")
        proteins = self.interpro.collect_human_proteins(entries)
        with get_connection(self.db_path) as conn:
            upsert_proteins(conn, proteins)
        return proteins

    def _phase2b_deduplicate(self, proteins: List[Protein]) -> List[Protein]:
        """Mark redundant proteins in the DB and return only the canonical representatives."""
        with get_connection(self.db_path) as conn:
            n_redundant = deduplicate_proteins(conn)
            redundant_ids = {
                row[0] for row in conn.execute(
                    "SELECT uniprot_id FROM proteins WHERE canonical_uniprot_id IS NOT NULL"
                ).fetchall()
            }
        canonical = [p for p in proteins if p.uniprot_id not in redundant_ids]
        logger.info(
            "Phase 2b: deduplication — %d redundant proteins marked; %d canonical proteins proceed to phase 3",
            n_redundant, len(canonical),
        )
        return canonical

    def _phase3_isoforms(
        self, proteins: List[Protein], report: CollectionReport
    ) -> List[Isoform]:
        total = len(proteins)
        logger.info("Phase 3: collecting isoforms for %d proteins from UniProt", total)
        all_isoforms: List[Isoform] = []

        for idx, protein in enumerate(proteins, 1):
            if idx % 25 == 0 or idx == 1 or idx == total:
                logger.info("  %d/%d proteins — %d isoforms so far", idx, total, len(all_isoforms))
            isoforms = self.uniprot.collect_isoforms(protein)
            if not isoforms:
                report.failed_proteins.append(protein.uniprot_id)
                continue
            with get_connection(self.db_path) as conn:
                upsert_isoforms(conn, isoforms)
            all_isoforms.extend(isoforms)

        return all_isoforms
