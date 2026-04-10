"""
Collect isoform data from UniProt for a list of proteins.

Extraction logic is ported from the ground-truth scripts that built the
original 407-protein database.  Key decisions:

- One Isoform row per isoform (canonical + alternatives).
- splice_variants stores the UniProt "Alternative sequence" features that
  apply to *this* isoform (positions + sequence changes).  This is the
  primary input for downstream AS-effect analysis.
- tim_barrel_location is populated only for the canonical isoform via
  InterPro; alternative isoforms inherit None and will be annotated during
  analysis once splice coordinates are mapped onto domain boundaries.
- exon_annotations is left None; it requires Ensembl coordinate mapping
  (a separate future collection step).
"""

import logging
from typing import Any, Dict, List, Optional

from ..api.interpro_client import InterProClient
from ..api.uniprot_client import UniProtClient
from ..models.entities import Isoform, Protein

logger = logging.getLogger(__name__)


class UniProtCollector:
    def __init__(
        self,
        uniprot: Optional[UniProtClient] = None,
        interpro: Optional[InterProClient] = None,
    ):
        self.uniprot = uniprot or UniProtClient()
        self.interpro = interpro or InterProClient()

    # ------------------------------------------------------------------
    # Public
    # ------------------------------------------------------------------

    def collect_isoforms(self, protein: Protein) -> List[Isoform]:
        """
        Fetch all isoforms for *protein* and return them as Isoform objects.

        Returns an empty list if UniProt has no data for this protein.
        """
        uid = protein.uniprot_id
        data = self.uniprot.get_protein(uid)
        if not data:
            logger.warning("No UniProt data for %s — skipping", uid)
            return []

        try:
            return self._extract_isoforms(data, protein)
        except Exception as e:
            logger.error("Extraction failed for %s: %s", uid, e)
            return []

    def collect_batch(
        self, proteins: List[Protein], log_every: int = 10
    ) -> List[Isoform]:
        """Collect isoforms for every protein; logs progress every *log_every* proteins."""
        all_isoforms: List[Isoform] = []
        total = len(proteins)
        for i, protein in enumerate(proteins, 1):
            isoforms = self.collect_isoforms(protein)
            all_isoforms.extend(isoforms)
            if i % log_every == 0 or i == total:
                logger.info(
                    "Progress: %d/%d proteins — %d isoforms collected so far",
                    i, total, len(all_isoforms),
                )
        return all_isoforms

    # ------------------------------------------------------------------
    # Extraction helpers
    # ------------------------------------------------------------------

    def _extract_isoforms(self, data: Dict[str, Any], protein: Protein) -> List[Isoform]:
        uid = protein.uniprot_id
        canonical_seq = data.get("sequence", {}).get("value", "")
        if not canonical_seq:
            logger.error("No canonical sequence for %s", uid)
            return []

        # Shared data extracted once per protein
        ensembl_gene_id = _extract_ensembl_gene_id(data)
        alphafold_id = _extract_alphafold_id(data)
        all_splice_features = _extract_all_splice_features(data)
        tim_barrel_loc = self._get_tim_barrel_location(uid, protein.tim_barrel_accession)

        # --- Canonical isoform ---
        canonical = Isoform(
            isoform_id=f"{uid}-1",
            uniprot_id=uid,
            is_canonical=True,
            sequence=canonical_seq,
            sequence_length=len(canonical_seq),
            splice_variants=[],   # canonical has no alternative sequences
            tim_barrel_location=tim_barrel_loc,
            ensembl_gene_id=ensembl_gene_id,
            alphafold_id=alphafold_id,
        )
        isoforms: List[Isoform] = [canonical]

        # --- Alternative isoforms ---
        alt_isoform_meta = _parse_alternative_products(data)
        for meta in alt_isoform_meta:
            isoform_id = meta["isoform_id"]
            if isoform_id == f"{uid}-1":
                continue  # already handled as canonical

            sequence = self.uniprot.get_isoform_sequence(isoform_id)
            if not sequence:
                logger.warning("No sequence for isoform %s — skipping", isoform_id)
                continue

            # Splice features whose VSP IDs are referenced by this isoform
            vsp_ids = set(meta.get("sequence_ids", []))
            splice_variants = [f for f in all_splice_features if f.get("featureId") in vsp_ids]

            isoforms.append(
                Isoform(
                    isoform_id=isoform_id,
                    uniprot_id=uid,
                    is_canonical=False,
                    sequence=sequence,
                    sequence_length=len(sequence),
                    splice_variants=splice_variants,
                    tim_barrel_location=None,  # computed during analysis
                    ensembl_gene_id=ensembl_gene_id,
                    alphafold_id=alphafold_id,
                )
            )

        return isoforms

    def _get_tim_barrel_location(
        self, uniprot_id: str, tim_barrel_accession: str
    ) -> Optional[Dict[str, Any]]:
        try:
            return self.interpro.get_domain_boundaries(uniprot_id, tim_barrel_accession)
        except Exception as e:
            logger.warning("Could not get TIM barrel location for %s: %s", uniprot_id, e)
            return None


# ---------------------------------------------------------------------------
# Pure extraction functions (stateless, easily unit-tested)
# ---------------------------------------------------------------------------

def _parse_alternative_products(data: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Parse the ALTERNATIVE PRODUCTS comment and return one dict per alternative isoform.

    Each dict: {isoform_id, name, sequence_status, sequence_ids}
    """
    for comment in data.get("comments", []):
        if comment.get("commentType") != "ALTERNATIVE PRODUCTS":
            continue
        results = []
        for iso in comment.get("isoforms", []):
            ids = iso.get("isoformIds", [])
            if not ids:
                continue
            results.append({
                "isoform_id": ids[0],
                "name": iso.get("name", {}).get("value", ""),
                "sequence_status": iso.get("isoformSequenceStatus", ""),
                "sequence_ids": iso.get("sequenceIds", []),
            })
        return results
    return []


def _extract_all_splice_features(data: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Return all features of type 'Alternative sequence' from the UniProt entry.

    Each feature describes a sequence change and which isoform(s) it applies to.
    Stored verbatim so downstream analysis has full fidelity.
    """
    return [
        {
            "featureId": f.get("featureId"),
            "location": f.get("location"),
            "description": f.get("description"),
            "evidences": f.get("evidences", []),
        }
        for f in data.get("features", [])
        if f.get("type") == "Alternative sequence"
    ]


def _extract_ensembl_gene_id(data: Dict[str, Any]) -> Optional[str]:
    """Return the first Ensembl gene ID from cross-references."""
    for ref in data.get("uniProtKBCrossReferences", []):
        if ref.get("database") == "Ensembl":
            return ref.get("id")
    return None


def _extract_alphafold_id(data: Dict[str, Any]) -> Optional[str]:
    """Return the AlphaFold DB accession from cross-references."""
    for ref in data.get("uniProtKBCrossReferences", []):
        if ref.get("database") == "AlphaFoldDB":
            return ref.get("id")
    return None
