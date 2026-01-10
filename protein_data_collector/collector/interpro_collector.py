"""
InterPro data collection service for PFAM families and human proteins.

This module provides high-level data collection orchestration for InterPro
database queries, focusing on TIM barrel protein families and human proteins.
"""

import asyncio
import logging
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass, field
from datetime import datetime

from ..api.interpro_client import InterProAPIClient, get_tim_barrel_entries, get_human_proteins_for_tim_barrel_entries
from ..models.entities import TIMBarrelEntryModel, InterProProteinModel
from ..config import get_config
from ..retry import get_retry_controller
from ..errors import DataError, APIError


@dataclass
class CollectionStats:
    """Statistics for data collection operations."""
    pfam_families_found: int = 0
    pfam_families_processed: int = 0
    human_proteins_found: int = 0
    human_proteins_processed: int = 0
    validation_errors: int = 0
    api_errors: int = 0
    start_time: Optional[datetime] = None
    end_time: Optional[datetime] = None
    
    @property
    def duration_seconds(self) -> float:
        """Calculate collection duration in seconds."""
        if self.start_time and self.end_time:
            return (self.end_time - self.start_time).total_seconds()
        return 0.0
    
    @property
    def success_rate(self) -> float:
        """Calculate overall success rate."""
        total_operations = self.pfam_families_processed + self.human_proteins_processed
        if total_operations == 0:
            return 0.0
        successful_operations = total_operations - self.validation_errors - self.api_errors
        return successful_operations / total_operations


@dataclass
class CollectionResult:
    """Result of data collection operation."""
    tim_barrel_entries: List[TIMBarrelEntryModel] = field(default_factory=list)
    human_proteins: List[InterProProteinModel] = field(default_factory=list)
    stats: CollectionStats = field(default_factory=CollectionStats)
    errors: List[str] = field(default_factory=list)


class InterProCollector:
    """
    High-level data collector for InterPro database queries.
    
    Orchestrates the collection of PFAM families with TIM barrel annotations
    and human proteins belonging to those families.
    """
    
    def __init__(self, client: Optional[InterProAPIClient] = None):
        """
        Initialize InterPro collector.
        
        Args:
            client: InterPro API client, creates new one if not provided
        """
        self.client = client
        self.config = get_config()
        self.retry_controller = get_retry_controller()
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
    
    async def __aenter__(self):
        """Async context manager entry."""
        if self.client is None:
            self.client = InterProAPIClient()
        return self
    
    async def __aexit__(self, exc_type, exc_val, exc_tb):
        """Async context manager exit."""
        if self.client:
            await self.client.close()
    
    async def collect_tim_barrel_entries(
        self, 
        page_size: int = 200
    ) -> Tuple[List[TIMBarrelEntryModel], CollectionStats]:
        """
        Collect all TIM barrel entries (both PFAM families and InterPro entries).
        
        Args:
            page_size: Number of results per page for API queries
            
        Returns:
            Tuple of (collected entries, collection statistics)
            
        Raises:
            APIError: For API-related errors
            DataError: For data validation errors
        """
        stats = CollectionStats(start_time=datetime.now())
        
        self.logger.info("Starting TIM barrel entry collection (PFAM families and InterPro entries)")
        
        try:
            if self.client is None:
                self.client = InterProAPIClient()
            
            # Collect TIM barrel entries (both PFAM and InterPro)
            entries = await get_tim_barrel_entries(self.client, page_size)
            
            stats.pfam_families_found = len([e for e in entries if e.is_pfam])
            stats.interpro_entries_found = len([e for e in entries if e.is_interpro])
            stats.pfam_families_processed = len(families)
            
            self.logger.info(
                "Completed PFAM family collection",
                extra={
                    "families_found": stats.pfam_families_found,
                    "families_processed": stats.pfam_families_processed
                }
            )
            
            return entries, stats
            
        except Exception as e:
            stats.api_errors += 1
            self.logger.error(
                "Failed to collect PFAM families: %s",
                str(e),
                extra={"error_type": type(e).__name__}
            )
            raise
        finally:
            stats.end_time = datetime.now()
    
    async def collect_human_proteins_for_entries(
        self, 
        tim_barrel_entries: List[TIMBarrelEntryModel],
        page_size: int = 200,
        organism: str = "Homo sapiens"
    ) -> Tuple[List[InterProProteinModel], CollectionStats]:
        """
        Collect human proteins for given PFAM families.
        
        Args:
            pfam_families: List of PFAM families to query
            page_size: Number of results per page for API queries
            organism: Target organism (default: Homo sapiens)
            
        Returns:
            Tuple of (collected proteins, collection statistics)
            
        Raises:
            APIError: For API-related errors
            DataError: For data validation errors
        """
        stats = CollectionStats(start_time=datetime.now())
        
        self.logger.info(
            "Starting human protein collection for %d PFAM families",
            len(pfam_families),
            extra={
                "pfam_families_count": len(pfam_families),
                "target_organism": organism
            }
        )
        
        try:
            if self.client is None:
                self.client = InterProAPIClient()
            
            # Collect human proteins for all families
            proteins = await get_human_proteins_for_pfam_families(
                pfam_families, 
                self.client, 
                page_size
            )
            
            # Filter proteins to ensure they are from the target organism
            filtered_proteins = []
            for protein in proteins:
                if protein.organism.lower() == organism.lower():
                    filtered_proteins.append(protein)
                else:
                    self.logger.warning(
                        "Protein %s has unexpected organism: %s (expected: %s)",
                        protein.uniprot_id,
                        protein.organism,
                        organism,
                        extra={
                            "protein_id": protein.uniprot_id,
                            "actual_organism": protein.organism,
                            "expected_organism": organism
                        }
                    )
                    stats.validation_errors += 1
            
            stats.human_proteins_found = len(proteins)
            stats.human_proteins_processed = len(filtered_proteins)
            
            self.logger.info(
                "Completed human protein collection",
                extra={
                    "proteins_found": stats.human_proteins_found,
                    "proteins_processed": stats.human_proteins_processed,
                    "validation_errors": stats.validation_errors,
                    "target_organism": organism
                }
            )
            
            return filtered_proteins, stats
            
        except Exception as e:
            stats.api_errors += 1
            self.logger.error(
                "Failed to collect human proteins: %s",
                str(e),
                extra={"error_type": type(e).__name__}
            )
            raise
        finally:
            stats.end_time = datetime.now()
    
    async def collect_tim_barrel_data(
        self, 
        page_size: int = 200,
        organism: str = "Homo sapiens"
    ) -> CollectionResult:
        """
        Perform complete TIM barrel data collection workflow.
        
        Collects PFAM families with TIM barrel annotations, then collects
        human proteins belonging to those families.
        
        Args:
            page_size: Number of results per page for API queries
            organism: Target organism (default: Homo sapiens)
            
        Returns:
            CollectionResult with families, proteins, and statistics
        """
        result = CollectionResult()
        overall_stats = CollectionStats(start_time=datetime.now())
        
        self.logger.info(
            "Starting complete TIM barrel data collection workflow",
            extra={
                "page_size": page_size,
                "target_organism": organism
            }
        )
        
        try:
            # Phase 1: Collect PFAM families
            self.logger.info("Phase 1: Collecting PFAM families with TIM barrel annotations")
            families, family_stats = await self.collect_tim_barrel_pfam_families(page_size)
            result.pfam_families = families
            
            # Update overall stats
            overall_stats.pfam_families_found = family_stats.pfam_families_found
            overall_stats.pfam_families_processed = family_stats.pfam_families_processed
            overall_stats.api_errors += family_stats.api_errors
            overall_stats.validation_errors += family_stats.validation_errors
            
            if not families:
                self.logger.warning("No PFAM families found with TIM barrel annotations")
                result.errors.append("No PFAM families found with TIM barrel annotations")
                return result
            
            # Phase 2: Collect human proteins
            self.logger.info("Phase 2: Collecting human proteins for PFAM families")
            proteins, protein_stats = await self.collect_human_proteins_for_entries(
                entries, page_size, organism
            )
            result.human_proteins = proteins
            
            # Update overall stats
            overall_stats.human_proteins_found = protein_stats.human_proteins_found
            overall_stats.human_proteins_processed = protein_stats.human_proteins_processed
            overall_stats.api_errors += protein_stats.api_errors
            overall_stats.validation_errors += protein_stats.validation_errors
            
            self.logger.info(
                "Completed TIM barrel data collection workflow",
                extra={
                    "pfam_families": len(result.pfam_families),
                    "human_proteins": len(result.human_proteins),
                    "success_rate": overall_stats.success_rate,
                    "duration_seconds": overall_stats.duration_seconds
                }
            )
            
        except Exception as e:
            error_msg = f"TIM barrel data collection failed: {str(e)}"
            result.errors.append(error_msg)
            self.logger.error(
                error_msg,
                extra={"error_type": type(e).__name__}
            )
            raise
        finally:
            overall_stats.end_time = datetime.now()
            result.stats = overall_stats
        
        return result
    
    async def validate_human_organism_filtering(
        self, 
        proteins: List[InterProProteinModel],
        expected_organism: str = "Homo sapiens"
    ) -> Tuple[List[InterProProteinModel], int]:
        """
        Validate that all proteins are from the expected organism.
        
        Args:
            proteins: List of proteins to validate
            expected_organism: Expected organism name
            
        Returns:
            Tuple of (valid proteins, error count)
        """
        valid_proteins = []
        error_count = 0
        
        for protein in proteins:
            if protein.organism.lower() == expected_organism.lower():
                valid_proteins.append(protein)
            else:
                error_count += 1
                self.logger.warning(
                    "Protein %s has unexpected organism: %s (expected: %s)",
                    protein.uniprot_id,
                    protein.organism,
                    expected_organism,
                    extra={
                        "protein_id": protein.uniprot_id,
                        "actual_organism": protein.organism,
                        "expected_organism": expected_organism
                    }
                )
        
        self.logger.info(
            "Organism filtering validation completed",
            extra={
                "total_proteins": len(proteins),
                "valid_proteins": len(valid_proteins),
                "organism_errors": error_count,
                "expected_organism": expected_organism
            }
        )
        
        return valid_proteins, error_count
    
    def get_pfam_family_summary(self, families: List[PfamFamilyModel]) -> Dict[str, Any]:
        """
        Generate summary statistics for collected PFAM families.
        
        Args:
            families: List of PFAM families
            
        Returns:
            Summary statistics dictionary
        """
        if not families:
            return {
                "total_families": 0,
                "accessions": [],
                "names": [],
                "tim_barrel_annotations": []
            }
        
        return {
            "total_families": len(families),
            "accessions": [family.accession for family in families],
            "names": [family.name for family in families],
            "tim_barrel_annotations": [family.tim_barrel_annotation for family in families],
            "unique_interpro_ids": list(set(
                family.interpro_id for family in families 
                if family.interpro_id
            ))
        }
    
    def get_protein_summary(self, proteins: List[InterProProteinModel]) -> Dict[str, Any]:
        """
        Generate summary statistics for collected proteins.
        
        Args:
            proteins: List of proteins
            
        Returns:
            Summary statistics dictionary
        """
        if not proteins:
            return {
                "total_proteins": 0,
                "unique_proteins": 0,
                "pfam_families": [],
                "organisms": []
            }
        
        unique_proteins = set(protein.uniprot_id for protein in proteins)
        pfam_families = set(protein.pfam_accession for protein in proteins)
        organisms = set(protein.organism for protein in proteins)
        
        return {
            "total_proteins": len(proteins),
            "unique_proteins": len(unique_proteins),
            "pfam_families": list(pfam_families),
            "organisms": list(organisms),
            "proteins_per_family": {
                family: len([p for p in proteins if p.pfam_accession == family])
                for family in pfam_families
            }
        }


# Convenience functions for common operations
async def collect_tim_barrel_families_and_proteins(
    page_size: int = 200,
    organism: str = "Homo sapiens"
) -> CollectionResult:
    """
    Convenience function to collect TIM barrel families and human proteins.
    
    Args:
        page_size: Number of results per page for API queries
        organism: Target organism (default: Homo sapiens)
        
    Returns:
        CollectionResult with families, proteins, and statistics
    """
    async with InterProCollector() as collector:
        return await collector.collect_tim_barrel_data(page_size, organism)


async def validate_organism_filtering(
    proteins: List[InterProProteinModel],
    expected_organism: str = "Homo sapiens"
) -> Tuple[List[InterProProteinModel], int]:
    """
    Convenience function to validate organism filtering.
    
    Args:
        proteins: List of proteins to validate
        expected_organism: Expected organism name
        
    Returns:
        Tuple of (valid proteins, error count)
    """
    async with InterProCollector() as collector:
        return await collector.validate_human_organism_filtering(proteins, expected_organism)