"""Query engine for retrieving protein and isoform data from the database."""

import json
import logging
import sqlite3
from typing import Any, Dict, List, Optional

from ..config import DOMAINS, ORGANISMS, DomainConfig, OrganismConfig
from ..database.connection import get_connection
from ..database.storage import get_counts

logger = logging.getLogger(__name__)


class QueryEngine:
    def __init__(
        self,
        db_path: Optional[str] = None,
        domain: str = "tim_barrel",
        organism: str = "homo_sapiens",
    ):
        self.db_path = db_path
        self.domain: DomainConfig   = DOMAINS[domain]
        self.org:    OrganismConfig = ORGANISMS[organism]

        self.entries_table  = self.domain.entries_table
        self.protein_table  = self.org.protein_table(self.domain)
        self.isoform_table  = self.org.isoform_table(self.domain)
        self.affected_table = self.org.affected_isoforms_table(self.domain)

    # ------------------------------------------------------------------
    # Domain entries
    # ------------------------------------------------------------------

    def get_all_families(self) -> List[Dict[str, Any]]:
        with get_connection(self.db_path) as conn:
            rows = conn.execute(
                f"SELECT * FROM {self.entries_table} ORDER BY accession"
            ).fetchall()
            return [dict(r) for r in rows]

    # ------------------------------------------------------------------
    # Proteins
    # ------------------------------------------------------------------

    def get_all_proteins(self) -> List[Dict[str, Any]]:
        with get_connection(self.db_path) as conn:
            rows = conn.execute(
                f"SELECT * FROM {self.protein_table} ORDER BY uniprot_id"
            ).fetchall()
            return [dict(r) for r in rows]

    def get_proteins_by_family(self, accession: str) -> List[Dict[str, Any]]:
        with get_connection(self.db_path) as conn:
            rows = conn.execute(
                f"SELECT * FROM {self.protein_table} WHERE tim_barrel_accession = ? ORDER BY uniprot_id",
                (accession,),
            ).fetchall()
            return [dict(r) for r in rows]

    def get_protein(self, uniprot_id: str) -> Optional[Dict[str, Any]]:
        with get_connection(self.db_path) as conn:
            row = conn.execute(
                f"SELECT * FROM {self.protein_table} WHERE uniprot_id = ?", (uniprot_id,)
            ).fetchone()
            return dict(row) if row else None

    # ------------------------------------------------------------------
    # Isoforms
    # ------------------------------------------------------------------

    def get_isoforms_for_protein(self, uniprot_id: str) -> List[Dict[str, Any]]:
        with get_connection(self.db_path) as conn:
            rows = conn.execute(
                f"SELECT * FROM {self.isoform_table} WHERE uniprot_id = ? ORDER BY is_canonical DESC, isoform_id",
                (uniprot_id,),
            ).fetchall()
            return [_deserialize_isoform(dict(r)) for r in rows]

    def get_all_isoforms(self) -> List[Dict[str, Any]]:
        with get_connection(self.db_path) as conn:
            rows = conn.execute(
                f"SELECT * FROM {self.isoform_table} ORDER BY uniprot_id, is_canonical DESC"
            ).fetchall()
            return [_deserialize_isoform(dict(r)) for r in rows]

    def get_proteins_with_alternative_isoforms(self) -> List[Dict[str, Any]]:
        """Return proteins that have at least one non-canonical isoform."""
        with get_connection(self.db_path) as conn:
            rows = conn.execute(
                f"""
                SELECT p.*, COUNT(i.isoform_id) as isoform_count
                FROM {self.protein_table} p
                JOIN {self.isoform_table} i ON p.uniprot_id = i.uniprot_id
                WHERE i.is_canonical = 0
                GROUP BY p.uniprot_id
                ORDER BY isoform_count DESC
                """
            ).fetchall()
            return [dict(r) for r in rows]

    def get_isoforms_with_domain(self) -> List[Dict[str, Any]]:
        """Return all isoforms that have a resolved domain location."""
        with get_connection(self.db_path) as conn:
            rows = conn.execute(
                f"SELECT * FROM {self.isoform_table} WHERE tim_barrel_location IS NOT NULL"
            ).fetchall()
            return [_deserialize_isoform(dict(r)) for r in rows]

    def get_isoforms_with_splice_variants(self) -> List[Dict[str, Any]]:
        """Return non-canonical isoforms that have at least one splice variant."""
        with get_connection(self.db_path) as conn:
            rows = conn.execute(
                f"""
                SELECT * FROM {self.isoform_table}
                WHERE is_canonical = 0
                  AND splice_variants IS NOT NULL
                  AND splice_variants != '[]'
                ORDER BY uniprot_id
                """
            ).fetchall()
            return [_deserialize_isoform(dict(r)) for r in rows]

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------

    def summary(self) -> Dict[str, Any]:
        with get_connection(self.db_path) as conn:
            counts = get_counts(
                conn,
                entries_table=self.entries_table,
                protein_table=self.protein_table,
                isoform_table=self.isoform_table,
            )
            reviewed = conn.execute(
                f"SELECT COUNT(*) FROM {self.protein_table} WHERE reviewed = 1"
            ).fetchone()[0]
            avg_len = conn.execute(
                f"SELECT AVG(sequence_length) FROM {self.isoform_table} WHERE is_canonical = 1"
            ).fetchone()[0]
        return {
            **counts,
            "proteins": counts[self.protein_table],
            "isoforms": counts[self.isoform_table],
            "reviewed_proteins": reviewed,
            "avg_canonical_sequence_length": round(avg_len, 1) if avg_len else None,
        }


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _deserialize_isoform(row: Dict[str, Any]) -> Dict[str, Any]:
    """Parse JSON columns back to Python objects."""
    for col in ("exon_annotations", "splice_variants", "tim_barrel_location"):
        val = row.get(col)
        if isinstance(val, str):
            try:
                row[col] = json.loads(val)
            except json.JSONDecodeError:
                pass
    return row
