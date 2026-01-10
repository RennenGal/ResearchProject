"""
UniProt protein isoform data collection service.

This module provides functionality to collect detailed isoform data for proteins
identified from InterPro, including sequence, exon annotations, TIM barrel locations,
and metadata using the unified UniProt API client.
"""

import asyncio
import logging
from typing import List, Dict, Any, Optional, Set
from dataclasses import dataclass, field
from datetime import datetime

from ..api.uniprot_client import UnifiedUniProtClient
from ..models.entities import InterProProteinModel, ProteinModel
from ..config import get_config
from ..errors import APIError, NetworkError, DataError, ValidationError, ErrorContext, create_error_context


@dataclass
class IsoformCollectionReport:
    """Report of protein isoform collection operation."""
    total_proteins_processed: int = 0
    total_isoforms_collected: int = 0
    successful_proteins: int = 0
    failed_proteins: int = 0
    validation_errors: List[str] = field(default_factory=list)
    api_errors: List[str] = field(default_factory=list)
    collection_duration: Optional[float] = None
    start_time: Optional[datetime] = None
    end_time: Optional[datetime] = None
    
    @property
    def success_rate(self) -> float:
        """Calculate success rate as percentage."""
        if self.total_proteins_processed == 0:
            return 0.0
        return (self.successful_proteins / self.total_proteins_processed) * 100.0
    
    @property
    def average_isoforms_per_protein(self) -> float:
        """Calculate average isoforms per successful protein."""
        if self.successful_proteins == 0:
            return 0.0
        return self.total_isoforms_collected / self.successful_proteins


class UniProtIsoformCollector:
    """
    Service for collecting protein isoform data from UniProt.
    
    Retrieves all isoforms for each InterPro protein, collecting sequence,
    exon annotations, TIM barrel locations, and metadata.
    """
    
    def __init__(self, client: Optional[UnifiedUniProtClient] = None):
        """
        Initialize UniProt isoform collector.
        
        Args:
            client: Unified UniProt client, creates new one if not provided
        """
        self.client = client
        self.config = get_config()
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        self._processed_proteins: Set[str] = set()
    
    async def __aenter__(self):
        """Async context manager entry."""
        if self.client is None:
            self.client = UnifiedUniProtClient()
            await self.client.__aenter__()
        return self
    
    async def __aexit__(self, exc_type, exc_val, exc_tb):
        """Async context manager exit."""
        if self.client is not None:
            await self.client.__aexit__(exc_type, exc_val, exc_tb)
    
    async def collect_protein_isoforms(self, interpro_protein: InterProProteinModel) -> List[ProteinModel]:
        """
        Collect all isoforms for a single InterPro protein.
        
        Args:
            interpro_protein: InterPro protein model to collect isoforms for
            
        Returns:
            List of ProteinModel instances for all isoforms
            
        Raises:
            APIError: For API-specific errors
            NetworkError: For network-related errors
            ValidationError: For data validation errors
        """
        uniprot_id = interpro_protein.uniprot_id
        
        # Skip if already processed
        if uniprot_id in self._processed_proteins:
            self.logger.debug(
                "Skipping already processed protein %s",
                uniprot_id,
                extra={"uniprot_id": uniprot_id}
            )
            return []
        
        self.logger.info(
            "Collecting isoforms for protein %s",
            uniprot_id,
            extra={
                "uniprot_id": uniprot_id,
                "pfam_accession": interpro_protein.tim_barrel_accession,
                "protein_name": interpro_protein.name
            }
        )
        
        try:
            # Get isoform data from UniProt
            isoforms_data = await self.client.get_protein_isoforms(uniprot_id)
            
            if not isoforms_data:
                self.logger.warning(
                    "No isoforms found for protein %s",
                    uniprot_id,
                    extra={"uniprot_id": uniprot_id}
                )
                return []
            
            # Parse and validate isoform data
            proteins = []
            method = await self.client._determine_access_method()
            
            for isoform_data in isoforms_data:
                try:
                    # Parse isoform data into ProteinModel
                    protein = self.client.parse_protein_isoform_data(
                        isoform_data, 
                        uniprot_id, 
                        method
                    )
                    
                    # Additional processing for exon count calculation
                    protein = self._process_exon_annotations(protein)
                    
                    # Validate the protein model
                    self._validate_protein_model(protein)
                    
                    proteins.append(protein)
                    
                    self.logger.debug(
                        "Successfully processed isoform %s for protein %s",
                        protein.isoform_id,
                        uniprot_id,
                        extra={
                            "uniprot_id": uniprot_id,
                            "isoform_id": protein.isoform_id,
                            "sequence_length": protein.sequence_length,
                            "exon_count": protein.exon_count
                        }
                    )
                    
                except (DataError, ValidationError) as e:
                    self.logger.warning(
                        "Failed to process isoform data for protein %s: %s",
                        uniprot_id,
                        str(e),
                        extra={
                            "uniprot_id": uniprot_id,
                            "isoform_data": isoform_data,
                            "error_type": type(e).__name__
                        }
                    )
                    continue
            
            # Mark protein as processed
            self._processed_proteins.add(uniprot_id)
            
            self.logger.info(
                "Collected %d isoforms for protein %s",
                len(proteins),
                uniprot_id,
                extra={
                    "uniprot_id": uniprot_id,
                    "isoform_count": len(proteins),
                    "pfam_accession": interpro_protein.tim_barrel_accession
                }
            )
            
            return proteins
            
        except Exception as e:
            # Create error context
            context = create_error_context(
                operation="collect_protein_isoforms",
                database="UniProt",
                entity_id=uniprot_id,
                entity_type="protein",
                additional_data={
                    "pfam_accession": interpro_protein.tim_barrel_accession,
                    "protein_name": interpro_protein.name
                }
            )
            
            self.logger.error(
                "Failed to collect isoforms for protein %s: %s",
                uniprot_id,
                str(e),
                extra={
                    "uniprot_id": uniprot_id,
                    "error_type": type(e).__name__,
                    "pfam_accession": interpro_protein.tim_barrel_accession
                }
            )
            
            # Re-raise with context
            if isinstance(e, (APIError, NetworkError, DataError, ValidationError)):
                raise
            else:
                raise APIError(f"Unexpected error collecting isoforms for {uniprot_id}: {str(e)}")
    
    def _process_exon_annotations(self, protein: ProteinModel) -> ProteinModel:
        """
        Process and calculate exon count from exon annotations.
        
        Args:
            protein: ProteinModel to process
            
        Returns:
            Updated ProteinModel with calculated exon count
        """
        if not protein.exon_annotations:
            return protein
        
        # Calculate exon count from annotations
        exon_count = 0
        
        if isinstance(protein.exon_annotations, dict):
            # Handle different annotation formats
            if "exons" in protein.exon_annotations:
                exons = protein.exon_annotations["exons"]
                if isinstance(exons, list):
                    exon_count = len(exons)
            elif "exon_count" in protein.exon_annotations:
                exon_count = protein.exon_annotations.get("exon_count", 0)
            elif "features" in protein.exon_annotations:
                # Count exon features
                features = protein.exon_annotations["features"]
                if isinstance(features, list):
                    exon_count = len([f for f in features if f.get("type") == "exon"])
        
        # Update protein model with calculated exon count
        if exon_count > 0 and protein.exon_count != exon_count:
            # Create new protein model with updated exon count
            protein_dict = protein.model_dump()
            protein_dict["exon_count"] = exon_count
            return ProteinModel(**protein_dict)
        
        return protein
    
    def _validate_protein_model(self, protein: ProteinModel) -> None:
        """
        Validate protein model data integrity.
        
        Args:
            protein: ProteinModel to validate
            
        Raises:
            ValidationError: If validation fails
        """
        # Basic validation is handled by Pydantic, but we can add additional checks
        
        # Validate sequence length consistency
        if protein.sequence and len(protein.sequence) != protein.sequence_length:
            raise ValidationError(
                f"Sequence length mismatch for {protein.isoform_id}: "
                f"expected {protein.sequence_length}, got {len(protein.sequence)}"
            )
        
        # Validate exon count consistency
        if protein.exon_count is not None and protein.exon_annotations:
            if isinstance(protein.exon_annotations, dict) and "exons" in protein.exon_annotations:
                actual_exon_count = len(protein.exon_annotations["exons"])
                if protein.exon_count != actual_exon_count:
                    raise ValidationError(
                        f"Exon count mismatch for {protein.isoform_id}: "
                        f"expected {protein.exon_count}, got {actual_exon_count}"
                    )
        
        # Validate TIM barrel location coordinates (only if location has meaningful data)
        if protein.tim_barrel_location and protein.sequence_length:
            location = protein.tim_barrel_location
            if isinstance(location, dict) and location:  # Check if dict is not empty
                start = location.get("start", 0)
                end = location.get("end", 0)
                # Only validate if we have actual coordinate values
                if start > 0 and end > 0:
                    if start > protein.sequence_length or end > protein.sequence_length:
                        raise ValidationError(
                            f"TIM barrel location out of bounds for {protein.isoform_id}: "
                            f"location {start}-{end} exceeds sequence length {protein.sequence_length}"
                        )
                    if start >= end:
                        raise ValidationError(
                            f"Invalid TIM barrel location for {protein.isoform_id}: "
                            f"start {start} >= end {end}"
                        )
    
    async def collect_isoforms_batch(
        self, 
        interpro_proteins: List[InterProProteinModel],
        batch_size: Optional[int] = None
    ) -> List[ProteinModel]:
        """
        Collect isoforms for multiple InterPro proteins in batches.
        
        Args:
            interpro_proteins: List of InterPro proteins to process
            batch_size: Number of proteins to process concurrently
            
        Returns:
            List of all collected ProteinModel instances
        """
        if not interpro_proteins:
            return []
        
        batch_size = batch_size or self.config.collection.batch_size
        all_proteins = []
        
        self.logger.info(
            "Starting batch collection of isoforms for %d proteins",
            len(interpro_proteins),
            extra={
                "total_proteins": len(interpro_proteins),
                "batch_size": batch_size
            }
        )
        
        # Process proteins in batches
        for i in range(0, len(interpro_proteins), batch_size):
            batch = interpro_proteins[i:i + batch_size]
            batch_num = (i // batch_size) + 1
            total_batches = (len(interpro_proteins) + batch_size - 1) // batch_size
            
            self.logger.info(
                "Processing batch %d/%d with %d proteins",
                batch_num,
                total_batches,
                len(batch),
                extra={
                    "batch_number": batch_num,
                    "total_batches": total_batches,
                    "batch_size": len(batch)
                }
            )
            
            # Create tasks for concurrent processing
            tasks = [
                self.collect_protein_isoforms(protein) 
                for protein in batch
            ]
            
            # Execute batch concurrently
            batch_results = await asyncio.gather(*tasks, return_exceptions=True)
            
            # Collect results and handle exceptions
            for protein, result in zip(batch, batch_results):
                if isinstance(result, Exception):
                    self.logger.error(
                        "Failed to process protein %s in batch: %s",
                        protein.uniprot_id,
                        str(result),
                        extra={
                            "uniprot_id": protein.uniprot_id,
                            "batch_number": batch_num,
                            "error_type": type(result).__name__
                        }
                    )
                else:
                    all_proteins.extend(result)
        
        self.logger.info(
            "Completed batch collection: %d total isoforms collected from %d proteins",
            len(all_proteins),
            len(interpro_proteins),
            extra={
                "total_proteins": len(interpro_proteins),
                "total_isoforms": len(all_proteins),
                "average_isoforms_per_protein": len(all_proteins) / len(interpro_proteins) if interpro_proteins else 0
            }
        )
        
        return all_proteins
    
    async def collect_isoforms_with_report(
        self, 
        interpro_proteins: List[InterProProteinModel],
        batch_size: Optional[int] = None
    ) -> tuple[List[ProteinModel], IsoformCollectionReport]:
        """
        Collect isoforms with detailed reporting.
        
        Args:
            interpro_proteins: List of InterPro proteins to process
            batch_size: Number of proteins to process concurrently
            
        Returns:
            Tuple of (collected proteins, collection report)
        """
        report = IsoformCollectionReport()
        report.start_time = datetime.now()
        report.total_proteins_processed = len(interpro_proteins)
        
        try:
            proteins = await self.collect_isoforms_batch(interpro_proteins, batch_size)
            
            # Calculate report statistics
            report.total_isoforms_collected = len(proteins)
            report.successful_proteins = len(self._processed_proteins)
            report.failed_proteins = report.total_proteins_processed - report.successful_proteins
            
            report.end_time = datetime.now()
            report.collection_duration = (report.end_time - report.start_time).total_seconds()
            
            self.logger.info(
                "Collection completed successfully",
                extra={
                    "total_proteins": report.total_proteins_processed,
                    "successful_proteins": report.successful_proteins,
                    "failed_proteins": report.failed_proteins,
                    "total_isoforms": report.total_isoforms_collected,
                    "success_rate": report.success_rate,
                    "duration_seconds": report.collection_duration
                }
            )
            
            return proteins, report
            
        except Exception as e:
            report.end_time = datetime.now()
            report.collection_duration = (report.end_time - report.start_time).total_seconds()
            report.api_errors.append(str(e))
            
            self.logger.error(
                "Collection failed with error: %s",
                str(e),
                extra={
                    "error_type": type(e).__name__,
                    "duration_seconds": report.collection_duration
                }
            )
            
            raise
    
    def get_collection_statistics(self) -> Dict[str, Any]:
        """
        Get current collection statistics.
        
        Returns:
            Dictionary with collection statistics
        """
        return {
            "processed_proteins_count": len(self._processed_proteins),
            "processed_protein_ids": list(self._processed_proteins),
            "client_status": self.client.get_access_method_status() if self.client else None
        }
    
    def reset_collection_state(self) -> None:
        """Reset collection state for fresh collection."""
        self._processed_proteins.clear()
        self.logger.info("Collection state reset")


# Convenience functions for common operations
async def collect_all_isoforms(
    interpro_proteins: List[InterProProteinModel],
    batch_size: Optional[int] = None,
    client: Optional[UnifiedUniProtClient] = None
) -> List[ProteinModel]:
    """
    Collect all isoforms for a list of InterPro proteins.
    
    Args:
        interpro_proteins: List of InterPro proteins to process
        batch_size: Number of proteins to process concurrently
        client: Unified client, creates new one if not provided
        
    Returns:
        List of all collected ProteinModel instances
    """
    async with UniProtIsoformCollector(client) as collector:
        return await collector.collect_isoforms_batch(interpro_proteins, batch_size)


async def collect_isoforms_with_detailed_report(
    interpro_proteins: List[InterProProteinModel],
    batch_size: Optional[int] = None,
    client: Optional[UnifiedUniProtClient] = None
) -> tuple[List[ProteinModel], IsoformCollectionReport]:
    """
    Collect isoforms with detailed reporting.
    
    Args:
        interpro_proteins: List of InterPro proteins to process
        batch_size: Number of proteins to process concurrently
        client: Unified client, creates new one if not provided
        
    Returns:
        Tuple of (collected proteins, collection report)
    """
    async with UniProtIsoformCollector(client) as collector:
        return await collector.collect_isoforms_with_report(interpro_proteins, batch_size)


async def collect_single_protein_isoforms(
    interpro_protein: InterProProteinModel,
    client: Optional[UnifiedUniProtClient] = None
) -> List[ProteinModel]:
    """
    Collect isoforms for a single InterPro protein.
    
    Args:
        interpro_protein: InterPro protein to process
        client: Unified client, creates new one if not provided
        
    Returns:
        List of collected ProteinModel instances
    """
    async with UniProtIsoformCollector(client) as collector:
        return await collector.collect_protein_isoforms(interpro_protein)