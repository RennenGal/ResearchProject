"""
Orchestrate the three-phase collection pipeline.

Phase 1  — Domain entries    : InterPro → {domain}_entries
Phase 2  — Proteins          : InterPro → {domain}_proteins[_{organism}]
Phase 2b — Deduplication     : mark redundant proteins
Phase 3  — Isoforms          : UniProt  → {domain}_isoforms[_{organism}]
"""

import json
import logging
from dataclasses import dataclass, field
from typing import List, Optional

from ..config import DOMAINS, ORGANISMS, DomainConfig, OrganismConfig
from ..database.connection import ensure_db, get_connection
from ..database.storage import (
    deduplicate_proteins,
    get_all_domain_entries,
    get_all_proteins,
    get_proteins_without_isoforms,
    upsert_domain_entries,
    upsert_isoforms,
    upsert_proteins,
)
from ..models.entities import Isoform, Protein, TIMBarrelEntry
from .interpro_collector import InterProCollector
from .uniprot_collector import UniProtCollector

logger = logging.getLogger(__name__)


@dataclass
class CollectionReport:
    domain_entries: int = 0
    proteins_collected: int = 0
    isoforms_collected: int = 0
    failed_proteins: List[str] = field(default_factory=list)

    @property
    def alternative_isoforms(self) -> int:
        return self.isoforms_collected - self.proteins_collected

    def summary(self) -> str:
        return (
            f"Domain entries     : {self.domain_entries}\n"
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
        domain: str = "tim_barrel",
        organism: str = "homo_sapiens",
    ):
        if domain not in DOMAINS:
            raise ValueError(f"Unknown domain {domain!r}. Supported: {list(DOMAINS)}")
        if organism not in ORGANISMS:
            raise ValueError(f"Unknown organism {organism!r}. Supported: {list(ORGANISMS)}")

        self.interpro = interpro_collector or InterProCollector()
        self.uniprot  = uniprot_collector  or UniProtCollector()
        self.db_path  = db_path
        self.domain: DomainConfig   = DOMAINS[domain]
        self.org:    OrganismConfig = ORGANISMS[organism]

        # Derive table names from domain + organism
        self.protein_table          = self.org.protein_table(self.domain)
        self.isoform_table          = self.org.isoform_table(self.domain)
        self.affected_isoforms_table = self.org.affected_isoforms_table(self.domain)

    # ------------------------------------------------------------------
    # Full pipeline
    # ------------------------------------------------------------------

    def run_full_collection(self) -> CollectionReport:
        """Run all phases from scratch."""
        ensure_db(self.db_path)
        report = CollectionReport()

        entries = self._phase1_domain_entries()
        report.domain_entries = len(entries)

        proteins = self._phase2_proteins(entries)
        proteins = self._phase2b_deduplicate(proteins)
        report.proteins_collected = len(proteins)

        isoforms = self._phase3_isoforms(proteins, report)
        report.isoforms_collected = len(isoforms)

        logger.info("Collection complete.\n%s", report.summary())
        return report

    def recollect_all_isoforms(self) -> CollectionReport:
        """Delete all isoform rows and re-fetch from UniProt."""
        ensure_db(self.db_path)
        with get_connection(self.db_path) as conn:
            conn.execute(f"DELETE FROM {self.isoform_table}")
            conn.commit()
            logger.info("Cleared %s; re-collecting from UniProt", self.isoform_table)
            proteins_raw = get_all_proteins(conn, table=self.protein_table)

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
        """Phase 3 only — collect isoforms for canonical proteins not yet in the DB."""
        ensure_db(self.db_path)
        report = CollectionReport()

        with get_connection(self.db_path) as conn:
            remaining_ids = get_proteins_without_isoforms(
                conn,
                protein_table=self.protein_table,
                isoform_table=self.isoform_table,
            )
            all_proteins = {
                p["uniprot_id"]: p
                for p in get_all_proteins(conn, table=self.protein_table)
            }

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
        """Populate tim_barrel_location for canonical isoforms where it is NULL."""
        ensure_db(self.db_path)

        with get_connection(self.db_path) as conn:
            rows = conn.execute(
                f"""
                SELECT i.isoform_id, i.uniprot_id, p.tim_barrel_accession
                FROM {self.isoform_table} i
                JOIN {self.protein_table} p ON i.uniprot_id = p.uniprot_id
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
                        f"UPDATE {self.isoform_table} "
                        f"SET tim_barrel_location = ? WHERE isoform_id = ?",
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

    def _phase1_domain_entries(self) -> List[TIMBarrelEntry]:
        logger.info(
            "Phase 1: collecting %s entries from InterPro", self.domain.display_name
        )

        # Re-use existing DB entries when available — domain families are
        # organism-independent, so there's no need to re-fetch on subsequent runs.
        with get_connection(self.db_path) as conn:
            existing = get_all_domain_entries(conn, table=self.domain.entries_table)
        if existing:
            logger.info(
                "Phase 1: found %d entries in DB — skipping InterPro fetch",
                len(existing),
            )
            return [TIMBarrelEntry(**row) for row in existing]

        entries = self.interpro.collect_domain_entries(
            annotation=self.domain.interpro_annotation,
            search=self.domain.interpro_search,
            extra_accessions=self.domain.extra_accessions,
        )
        with get_connection(self.db_path) as conn:
            upsert_domain_entries(conn, entries, table=self.domain.entries_table)
        return entries

    def _phase2_proteins(self, entries: List[TIMBarrelEntry]) -> List[Protein]:
        logger.info(
            "Phase 2: collecting %s proteins from InterPro", self.org.display_name
        )
        proteins = self.interpro.collect_proteins(
            entries,
            organism=self.org.display_name,
            taxon_id=self.org.taxon_id,
        )
        with get_connection(self.db_path) as conn:
            upsert_proteins(conn, proteins, table=self.protein_table)
        return proteins

    def _phase2b_deduplicate(self, proteins: List[Protein]) -> List[Protein]:
        """Mark redundant proteins in the DB and return only the canonical representatives."""
        with get_connection(self.db_path) as conn:
            n_redundant = deduplicate_proteins(
                conn,
                protein_table=self.protein_table,
                isoform_table=self.isoform_table,
            )
            redundant_ids = {
                row[0] for row in conn.execute(
                    f"SELECT uniprot_id FROM {self.protein_table} "
                    f"WHERE canonical_uniprot_id IS NOT NULL"
                ).fetchall()
            }
        canonical = [p for p in proteins if p.uniprot_id not in redundant_ids]
        logger.info(
            "Phase 2b: deduplication — %d redundant proteins marked; "
            "%d canonical proteins proceed to phase 3",
            n_redundant, len(canonical),
        )
        return canonical

    def _phase3_isoforms(
        self, proteins: List[Protein], report: CollectionReport
    ) -> List[Isoform]:
        total = len(proteins)
        logger.info(
            "Phase 3: collecting isoforms for %d proteins from UniProt", total
        )
        all_isoforms: List[Isoform] = []

        for idx, protein in enumerate(proteins, 1):
            if idx % 25 == 0 or idx == 1 or idx == total:
                logger.info(
                    "  %d/%d proteins — %d isoforms so far",
                    idx, total, len(all_isoforms),
                )
            isoforms = self.uniprot.collect_isoforms(protein)
            if not isoforms:
                report.failed_proteins.append(protein.uniprot_id)
                continue
            with get_connection(self.db_path) as conn:
                upsert_isoforms(conn, isoforms, table=self.isoform_table)
            all_isoforms.extend(isoforms)

        return all_isoforms
