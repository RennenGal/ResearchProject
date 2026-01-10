"""
Query engine for protein data retrieval and analysis.

This module provides comprehensive search and filtering capabilities for PFAM families,
TIM barrel features, and protein identifiers with result formatting.
"""

import logging
from typing import List, Dict, Any, Optional, Union
from dataclasses import dataclass
from sqlalchemy.orm import Session, joinedload
from sqlalchemy import and_, or_, func, text

from ..database.connection import get_db_session
from ..database.schema import TIMBarrelEntry, InterProProtein, Protein
from ..models.entities import TIMBarrelEntryModel, InterProProteinModel, ProteinModel

logger = logging.getLogger(__name__)


@dataclass
class QueryFilters:
    """Data class for query filtering parameters."""
    pfam_family: Optional[str] = None
    protein_id: Optional[str] = None
    organism: Optional[str] = None
    min_sequence_length: Optional[int] = None
    max_sequence_length: Optional[int] = None
    min_exon_count: Optional[int] = None
    max_exon_count: Optional[int] = None
    has_tim_barrel: Optional[bool] = None
    tim_barrel_confidence: Optional[float] = None


@dataclass
class QueryResult:
    """Data class for formatted query results."""
    pfam_families: List[Dict[str, Any]]
    proteins: List[Dict[str, Any]]
    isoforms: List[Dict[str, Any]]
    total_count: int
    query_metadata: Dict[str, Any]


class QueryEngine:
    """
    Query engine for retrieving and analyzing protein data.
    
    Provides search methods for PFAM families, TIM barrel features, and protein identifiers
    with comprehensive filtering and result formatting capabilities.
    """
    
    def __init__(self):
        """Initialize the query engine."""
        self.logger = logging.getLogger(__name__)
    
    def search_by_pfam_family(self, pfam_id: str, include_isoforms: bool = True) -> QueryResult:
        """
        Search for all proteins and isoforms belonging to a specific PFAM family.
        
        Args:
            pfam_id: PFAM family accession identifier
            include_isoforms: Whether to include detailed isoform data
            
        Returns:
            QueryResult containing PFAM family, proteins, and isoforms
            
        Requirements: 5.1
        """
        self.logger.info(f"Searching for PFAM family: {pfam_id}")
        
        with get_db_session() as session:
            # Get TIM barrel entry
            tim_barrel_entry = session.query(TIMBarrelEntry).filter(
                TIMBarrelEntry.accession == pfam_id
            ).first()
            
            if not pfam_family:
                return QueryResult(
                    pfam_families=[],
                    proteins=[],
                    isoforms=[],
                    total_count=0,
                    query_metadata={"pfam_id": pfam_id, "found": False}
                )
            
            # Get proteins in this family
            proteins_query = session.query(InterProProtein).filter(
                InterProProtein.pfam_accession == pfam_id
            )
            
            if include_isoforms:
                proteins_query = proteins_query.options(joinedload(InterProProtein.isoforms))
            
            proteins = proteins_query.all()
            
            # Format results
            pfam_data = self._format_pfam_family(pfam_family)
            protein_data = [self._format_interpro_protein(p) for p in proteins]
            
            isoform_data = []
            if include_isoforms:
                for protein in proteins:
                    isoform_data.extend([self._format_protein_isoform(iso) for iso in protein.isoforms])
            
            total_count = len(proteins) + len(isoform_data)
            
            return QueryResult(
                pfam_families=[pfam_data],
                proteins=protein_data,
                isoforms=isoform_data,
                total_count=total_count,
                query_metadata={
                    "pfam_id": pfam_id,
                    "found": True,
                    "protein_count": len(proteins),
                    "isoform_count": len(isoform_data)
                }
            )
    
    def search_by_tim_barrel_features(self, criteria: Dict[str, Any]) -> QueryResult:
        """
        Search for proteins by TIM barrel structural features.
        
        Args:
            criteria: Dictionary containing TIM barrel search criteria:
                - min_confidence: Minimum confidence score
                - max_confidence: Maximum confidence score
                - has_location: Whether TIM barrel location is present
                
        Returns:
            QueryResult containing matching proteins and isoforms
            
        Requirements: 5.2
        """
        self.logger.info(f"Searching by TIM barrel features: {criteria}")
        
        with get_db_session() as session:
            query = session.query(Protein).options(
                joinedload(Protein.parent_protein).joinedload(InterProProtein.pfam_family)
            )
            
            # Apply TIM barrel filtering
            conditions = []
            
            if criteria.get('has_location', True):
                conditions.append(Protein.tim_barrel_location.isnot(None))
            
            if 'min_confidence' in criteria:
                conditions.append(
                    text("JSON_EXTRACT(tim_barrel_location, '$.confidence') >= :min_conf")
                    .bindparam(min_conf=criteria['min_confidence'])
                )
            
            if 'max_confidence' in criteria:
                conditions.append(
                    text("JSON_EXTRACT(tim_barrel_location, '$.confidence') <= :max_conf")
                    .bindparam(max_conf=criteria['max_confidence'])
                )
            
            if conditions:
                query = query.filter(and_(*conditions))
            
            isoforms = query.all()
            
            # Group by PFAM families and proteins
            pfam_families = {}
            proteins = {}
            
            for isoform in isoforms:
                pfam_family = isoform.parent_protein.pfam_family
                if pfam_family.accession not in pfam_families:
                    pfam_families[pfam_family.accession] = pfam_family
                
                if isoform.parent_protein.uniprot_id not in proteins:
                    proteins[isoform.parent_protein.uniprot_id] = isoform.parent_protein
            
            # Format results
            pfam_data = [self._format_pfam_family(pf) for pf in pfam_families.values()]
            protein_data = [self._format_interpro_protein(p) for p in proteins.values()]
            isoform_data = [self._format_protein_isoform(iso) for iso in isoforms]
            
            return QueryResult(
                pfam_families=pfam_data,
                proteins=protein_data,
                isoforms=isoform_data,
                total_count=len(isoforms),
                query_metadata={
                    "criteria": criteria,
                    "pfam_family_count": len(pfam_families),
                    "protein_count": len(proteins),
                    "isoform_count": len(isoforms)
                }
            )
    
    def get_protein_isoforms(self, protein_id: str) -> QueryResult:
        """
        Retrieve all isoforms for a specific protein.
        
        Args:
            protein_id: UniProt protein identifier
            
        Returns:
            QueryResult containing protein and all its isoforms
            
        Requirements: 5.3
        """
        self.logger.info(f"Getting isoforms for protein: {protein_id}")
        
        with get_db_session() as session:
            # Get protein with isoforms
            protein = session.query(InterProProtein).options(
                joinedload(InterProProtein.isoforms),
                joinedload(InterProProtein.pfam_family)
            ).filter(InterProProtein.uniprot_id == protein_id).first()
            
            if not protein:
                return QueryResult(
                    pfam_families=[],
                    proteins=[],
                    isoforms=[],
                    total_count=0,
                    query_metadata={"protein_id": protein_id, "found": False}
                )
            
            # Format results
            pfam_data = self._format_pfam_family(protein.pfam_family)
            protein_data = self._format_interpro_protein(protein)
            isoform_data = [self._format_protein_isoform(iso) for iso in protein.isoforms]
            
            return QueryResult(
                pfam_families=[pfam_data],
                proteins=[protein_data],
                isoforms=isoform_data,
                total_count=len(protein.isoforms),
                query_metadata={
                    "protein_id": protein_id,
                    "found": True,
                    "isoform_count": len(protein.isoforms)
                }
            )
    
    def filter_proteins(self, filters: QueryFilters, limit: Optional[int] = None) -> QueryResult:
        """
        Filter proteins by various criteria.
        
        Args:
            filters: QueryFilters object containing filtering criteria
            limit: Maximum number of results to return
            
        Returns:
            QueryResult containing filtered proteins and isoforms
            
        Requirements: 5.4
        """
        self.logger.info(f"Filtering proteins with criteria: {filters}")
        
        with get_db_session() as session:
            query = session.query(Protein).options(
                joinedload(Protein.parent_protein).joinedload(InterProProtein.pfam_family)
            )
            
            # Apply filters
            conditions = []
            
            if filters.pfam_family:
                conditions.append(
                    Protein.parent_protein.has(
                        InterProProtein.pfam_accession == filters.pfam_family
                    )
                )
            
            if filters.protein_id:
                conditions.append(Protein.parent_protein_id == filters.protein_id)
            
            if filters.organism:
                conditions.append(Protein.organism == filters.organism)
            
            if filters.min_sequence_length:
                conditions.append(Protein.sequence_length >= filters.min_sequence_length)
            
            if filters.max_sequence_length:
                conditions.append(Protein.sequence_length <= filters.max_sequence_length)
            
            if filters.min_exon_count:
                conditions.append(Protein.exon_count >= filters.min_exon_count)
            
            if filters.max_exon_count:
                conditions.append(Protein.exon_count <= filters.max_exon_count)
            
            if filters.has_tim_barrel is not None:
                if filters.has_tim_barrel:
                    conditions.append(Protein.tim_barrel_location.isnot(None))
                else:
                    conditions.append(Protein.tim_barrel_location.is_(None))
            
            if filters.tim_barrel_confidence:
                conditions.append(
                    text("JSON_EXTRACT(tim_barrel_location, '$.confidence') >= :confidence")
                    .bindparam(confidence=filters.tim_barrel_confidence)
                )
            
            if conditions:
                query = query.filter(and_(*conditions))
            
            if limit:
                query = query.limit(limit)
            
            isoforms = query.all()
            
            # Group by PFAM families and proteins
            pfam_families = {}
            proteins = {}
            
            for isoform in isoforms:
                pfam_family = isoform.parent_protein.pfam_family
                if pfam_family.accession not in pfam_families:
                    pfam_families[pfam_family.accession] = pfam_family
                
                if isoform.parent_protein.uniprot_id not in proteins:
                    proteins[isoform.parent_protein.uniprot_id] = isoform.parent_protein
            
            # Format results
            pfam_data = [self._format_pfam_family(pf) for pf in pfam_families.values()]
            protein_data = [self._format_interpro_protein(p) for p in proteins.values()]
            isoform_data = [self._format_protein_isoform(iso) for iso in isoforms]
            
            return QueryResult(
                pfam_families=pfam_data,
                proteins=protein_data,
                isoforms=isoform_data,
                total_count=len(isoforms),
                query_metadata={
                    "filters": filters.__dict__,
                    "limit": limit,
                    "pfam_family_count": len(pfam_families),
                    "protein_count": len(proteins),
                    "isoform_count": len(isoforms)
                }
            )
    
    def get_summary_statistics(self) -> Dict[str, Any]:
        """
        Get summary statistics about the collected data.
        
        Returns:
            Dictionary containing data statistics
            
        Requirements: 5.5
        """
        self.logger.info("Getting summary statistics")
        
        with get_db_session() as session:
            # Count entities
            tim_barrel_count = session.query(TIMBarrelEntry).count()
            protein_count = session.query(InterProProtein).count()
            isoform_count = session.query(Protein).count()
            
            # Get sequence length statistics
            seq_stats = session.query(
                func.min(Protein.sequence_length).label('min_length'),
                func.max(Protein.sequence_length).label('max_length'),
                func.avg(Protein.sequence_length).label('avg_length')
            ).first()
            
            # Get exon count statistics
            exon_stats = session.query(
                func.min(Protein.exon_count).label('min_exons'),
                func.max(Protein.exon_count).label('max_exons'),
                func.avg(Protein.exon_count).label('avg_exons')
            ).filter(Protein.exon_count.isnot(None)).first()
            
            # Count proteins with TIM barrel locations
            tim_barrel_count = session.query(Protein).filter(
                Protein.tim_barrel_location.isnot(None)
            ).count()
            
            return {
                "pfam_families": pfam_count,
                "proteins": protein_count,
                "isoforms": isoform_count,
                "sequence_length": {
                    "min": seq_stats.min_length if seq_stats else None,
                    "max": seq_stats.max_length if seq_stats else None,
                    "avg": float(seq_stats.avg_length) if seq_stats and seq_stats.avg_length else None
                },
                "exon_count": {
                    "min": exon_stats.min_exons if exon_stats else None,
                    "max": exon_stats.max_exons if exon_stats else None,
                    "avg": float(exon_stats.avg_exons) if exon_stats and exon_stats.avg_exons else None
                },
                "tim_barrel_annotations": tim_barrel_count,
                "tim_barrel_coverage": (tim_barrel_count / isoform_count * 100) if isoform_count > 0 else 0
            }
    
    def _format_tim_barrel_entry(self, tim_barrel_entry: TIMBarrelEntry) -> Dict[str, Any]:
        """Format TIM barrel entry for display."""
        return {
            "accession": pfam_family.accession,
            "name": pfam_family.name,
            "description": pfam_family.description,
            "tim_barrel_annotation": pfam_family.tim_barrel_annotation,
            "created_at": pfam_family.created_at.isoformat() if pfam_family.created_at else None
        }
    
    def _format_interpro_protein(self, protein: InterProProtein) -> Dict[str, Any]:
        """Format InterPro protein for display."""
        return {
            "uniprot_id": protein.uniprot_id,
            "pfam_accession": protein.pfam_accession,
            "name": protein.name,
            "organism": protein.organism,
            "created_at": protein.created_at.isoformat() if protein.created_at else None
        }
    
    def _format_protein_isoform(self, isoform: Protein) -> Dict[str, Any]:
        """Format protein isoform for display with all required fields."""
        return {
            "isoform_id": isoform.isoform_id,
            "parent_protein_id": isoform.parent_protein_id,
            "sequence": isoform.sequence,
            "sequence_length": isoform.sequence_length,
            "exon_annotations": isoform.exon_annotations,
            "exon_count": isoform.exon_count,
            "tim_barrel_location": isoform.tim_barrel_location,
            "organism": isoform.organism,
            "name": isoform.name,
            "description": isoform.description,
            "created_at": isoform.created_at.isoformat() if isoform.created_at else None
        }