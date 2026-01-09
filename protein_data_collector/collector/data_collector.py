"""
Main data collection orchestration service.

This module provides the primary DataCollector class that orchestrates the complete
three-phase collection workflow: PFAM families → InterPro proteins → UniProt isoforms.
Includes progress tracking, status reporting, and collection resume functionality.
"""

import asyncio
import json
import logging
from datetime import datetime, timedelta
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass, field
from pathlib import Path

from ..models.entities import PfamFamilyModel, InterProProteinModel, ProteinModel
from ..database.connection import get_database_manager, get_db_transaction
from ..database.schema import PfamFamily, InterProProtein, Protein
from ..database.storage import DatabaseStorage, StorageResult
from .interpro_collector import InterProCollector, CollectionResult as InterProResult
from .uniprot_collector import UniProtIsoformCollector, IsoformCollectionReport
from ..config import get_config
from ..errors import DataError, APIError, ValidationError


@dataclass
class CollectionProgress:
    """Tracks progress of data collection operation."""
    phase: str = "not_started"  # not_started, pfam_families, interpro_proteins, uniprot_isoforms, storage, completed
    pfam_families_collected: int = 0
    interpro_proteins_collected: int = 0
    uniprot_isoforms_collected: int = 0
    pfam_families_stored: int = 0
    interpro_proteins_stored: int = 0
    uniprot_isoforms_stored: int = 0
    start_time: Optional[datetime] = None
    last_checkpoint: Optional[datetime] = None
    errors: List[str] = field(default_factory=list)
    
    @property
    def total_entities_collected(self) -> int:
        """Total entities collected across all phases."""
        return (self.pfam_families_collected + 
                self.interpro_proteins_collected + 
                self.uniprot_isoforms_collected)
    
    @property
    def total_entities_stored(self) -> int:
        """Total entities stored in database."""
        return (self.pfam_families_stored + 
                self.interpro_proteins_stored + 
                self.uniprot_isoforms_stored)
    
    @property
    def duration_seconds(self) -> float:
        """Collection duration in seconds."""
        if self.start_time:
            end_time = self.last_checkpoint or datetime.now()
            return (end_time - self.start_time).total_seconds()
        return 0.0
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "phase": self.phase,
            "pfam_families_collected": self.pfam_families_collected,
            "interpro_proteins_collected": self.interpro_proteins_collected,
            "uniprot_isoforms_collected": self.uniprot_isoforms_collected,
            "pfam_families_stored": self.pfam_families_stored,
            "interpro_proteins_stored": self.interpro_proteins_stored,
            "uniprot_isoforms_stored": self.uniprot_isoforms_stored,
            "start_time": self.start_time.isoformat() if self.start_time else None,
            "last_checkpoint": self.last_checkpoint.isoformat() if self.last_checkpoint else None,
            "errors": self.errors,
            "total_entities_collected": self.total_entities_collected,
            "total_entities_stored": self.total_entities_stored,
            "duration_seconds": self.duration_seconds
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'CollectionProgress':
        """Create from dictionary."""
        progress = cls()
        progress.phase = data.get("phase", "not_started")
        progress.pfam_families_collected = data.get("pfam_families_collected", 0)
        progress.interpro_proteins_collected = data.get("interpro_proteins_collected", 0)
        progress.uniprot_isoforms_collected = data.get("uniprot_isoforms_collected", 0)
        progress.pfam_families_stored = data.get("pfam_families_stored", 0)
        progress.interpro_proteins_stored = data.get("interpro_proteins_stored", 0)
        progress.uniprot_isoforms_stored = data.get("uniprot_isoforms_stored", 0)
        progress.errors = data.get("errors", [])
        
        if data.get("start_time"):
            progress.start_time = datetime.fromisoformat(data["start_time"])
        if data.get("last_checkpoint"):
            progress.last_checkpoint = datetime.fromisoformat(data["last_checkpoint"])
        
        return progress


@dataclass
class CollectionReport:
    """Complete collection operation report."""
    progress: CollectionProgress = field(default_factory=CollectionProgress)
    pfam_families: List[PfamFamilyModel] = field(default_factory=list)
    interpro_proteins: List[InterProProteinModel] = field(default_factory=list)
    uniprot_isoforms: List[ProteinModel] = field(default_factory=list)
    validation_errors: List[str] = field(default_factory=list)
    api_errors: List[str] = field(default_factory=list)
    storage_errors: List[str] = field(default_factory=list)
    
    @property
    def success_rate(self) -> float:
        """Calculate overall success rate."""
        total_operations = (len(self.pfam_families) + 
                          len(self.interpro_proteins) + 
                          len(self.uniprot_isoforms))
        if total_operations == 0:
            return 0.0
        
        total_errors = (len(self.validation_errors) + 
                       len(self.api_errors) + 
                       len(self.storage_errors))
        successful_operations = total_operations - total_errors
        return (successful_operations / total_operations) * 100.0
    
    def to_summary_dict(self) -> Dict[str, Any]:
        """Generate summary dictionary."""
        return {
            "collection_phase": self.progress.phase,
            "duration_seconds": self.progress.duration_seconds,
            "pfam_families_collected": len(self.pfam_families),
            "interpro_proteins_collected": len(self.interpro_proteins),
            "uniprot_isoforms_collected": len(self.uniprot_isoforms),
            "entities_stored": self.progress.total_entities_stored,
            "success_rate": self.success_rate,
            "validation_errors": len(self.validation_errors),
            "api_errors": len(self.api_errors),
            "storage_errors": len(self.storage_errors),
            "total_errors": len(self.validation_errors) + len(self.api_errors) + len(self.storage_errors)
        }


class DataCollector:
    """
    Main data collection orchestration service.
    
    Coordinates the three-phase collection workflow:
    1. PFAM families with TIM barrel annotations from InterPro
    2. Human proteins belonging to those families from InterPro
    3. Protein isoforms with detailed data from UniProt
    
    Provides progress tracking, status reporting, and resume functionality.
    """
    
    def __init__(self, 
                 progress_file: Optional[str] = None,
                 batch_size: Optional[int] = None):
        """
        Initialize data collector.
        
        Args:
            progress_file: Path to progress file for resume functionality
            batch_size: Batch size for processing operations
        """
        self.config = get_config()
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        self.db_manager = get_database_manager()
        
        # Configuration
        self.batch_size = batch_size or self.config.collection.batch_size
        self.progress_file = Path(progress_file) if progress_file else None
        
        # Collection state
        self.progress = CollectionProgress()
        
        # Initialize collectors immediately for direct method access
        self.interpro_collector = InterProCollector()
        self.uniprot_collector = UniProtIsoformCollector()
        self.storage = DatabaseStorage(self.batch_size)
        
        # Resume functionality
        self._load_progress()
    
    def _load_progress(self) -> None:
        """Load progress from file if it exists."""
        if self.progress_file and self.progress_file.exists():
            try:
                with open(self.progress_file, 'r') as f:
                    progress_data = json.load(f)
                self.progress = CollectionProgress.from_dict(progress_data)
                
                self.logger.info(
                    "Loaded collection progress from file",
                    extra={
                        "progress_file": str(self.progress_file),
                        "current_phase": self.progress.phase,
                        "entities_collected": self.progress.total_entities_collected,
                        "entities_stored": self.progress.total_entities_stored
                    }
                )
            except Exception as e:
                self.logger.warning(
                    "Failed to load progress file, starting fresh: %s",
                    str(e),
                    extra={"progress_file": str(self.progress_file)}
                )
                self.progress = CollectionProgress()
    
    def _save_progress(self) -> None:
        """Save current progress to file."""
        if self.progress_file:
            try:
                self.progress.last_checkpoint = datetime.now()
                
                # Ensure directory exists
                self.progress_file.parent.mkdir(parents=True, exist_ok=True)
                
                with open(self.progress_file, 'w') as f:
                    json.dump(self.progress.to_dict(), f, indent=2)
                
                self.logger.debug(
                    "Saved collection progress to file",
                    extra={
                        "progress_file": str(self.progress_file),
                        "current_phase": self.progress.phase
                    }
                )
            except Exception as e:
                self.logger.error(
                    "Failed to save progress file: %s",
                    str(e),
                    extra={"progress_file": str(self.progress_file)}
                )
    
    async def __aenter__(self):
        """Async context manager entry."""
        # Collectors are already initialized, just ensure they're ready
        await self.interpro_collector.__aenter__()
        await self.uniprot_collector.__aenter__()
        
        return self
    
    async def __aexit__(self, exc_type, exc_val, exc_tb):
        """Async context manager exit."""
        if self.interpro_collector:
            await self.interpro_collector.__aexit__(exc_type, exc_val, exc_tb)
        if self.uniprot_collector:
            await self.uniprot_collector.__aexit__(exc_type, exc_val, exc_tb)
    
    async def collect_pfam_families(self, page_size: int = 200) -> List[PfamFamilyModel]:
        """
        Phase 1: Collect PFAM families with TIM barrel annotations.
        
        Args:
            page_size: Number of results per page for API queries
            
        Returns:
            List of collected PFAM families
        """
        if self.progress.phase in ["pfam_families", "interpro_proteins", "uniprot_isoforms", "storage", "completed"]:
            self.logger.info("PFAM families phase already completed, skipping")
            return []
        
        self.logger.info("Starting Phase 1: PFAM family collection")
        self.progress.phase = "pfam_families"
        self.progress.start_time = self.progress.start_time or datetime.now()
        self._save_progress()
        
        try:
            families, stats = await self.interpro_collector.collect_tim_barrel_pfam_families(page_size)
            
            self.progress.pfam_families_collected = len(families)
            self.logger.info(
                "Phase 1 completed: collected %d PFAM families",
                len(families),
                extra={
                    "families_collected": len(families),
                    "duration_seconds": stats.duration_seconds,
                    "success_rate": stats.success_rate
                }
            )
            
            self._save_progress()
            return families
            
        except Exception as e:
            error_msg = f"Phase 1 failed: {str(e)}"
            self.progress.errors.append(error_msg)
            self.logger.error(error_msg, extra={"error_type": type(e).__name__})
            self._save_progress()
            raise
    
    async def collect_interpro_proteins(self, 
                                      pfam_families: List[PfamFamilyModel],
                                      page_size: int = 200) -> List[InterProProteinModel]:
        """
        Phase 2: Collect human proteins for PFAM families.
        
        Args:
            pfam_families: PFAM families to collect proteins for
            page_size: Number of results per page for API queries
            
        Returns:
            List of collected InterPro proteins
        """
        if self.progress.phase in ["interpro_proteins", "uniprot_isoforms", "storage", "completed"]:
            self.logger.info("InterPro proteins phase already completed, skipping")
            return []
        
        self.logger.info("Starting Phase 2: InterPro protein collection")
        self.progress.phase = "interpro_proteins"
        self._save_progress()
        
        try:
            proteins, stats = await self.interpro_collector.collect_human_proteins_for_families(
                pfam_families, page_size
            )
            
            self.progress.interpro_proteins_collected = len(proteins)
            self.logger.info(
                "Phase 2 completed: collected %d InterPro proteins",
                len(proteins),
                extra={
                    "proteins_collected": len(proteins),
                    "duration_seconds": stats.duration_seconds,
                    "success_rate": stats.success_rate
                }
            )
            
            self._save_progress()
            return proteins
            
        except Exception as e:
            error_msg = f"Phase 2 failed: {str(e)}"
            self.progress.errors.append(error_msg)
            self.logger.error(error_msg, extra={"error_type": type(e).__name__})
            self._save_progress()
            raise
    
    async def collect_human_proteins(self, 
                                   pfam_families: List[PfamFamilyModel],
                                   page_size: int = 200) -> List[InterProProteinModel]:
        """
        Alias for collect_interpro_proteins for backward compatibility.
        
        Args:
            pfam_families: PFAM families to collect proteins for
            page_size: Number of results per page for API queries
            
        Returns:
            List of collected InterPro proteins
        """
        return await self.collect_interpro_proteins(pfam_families, page_size)
    
    async def collect_uniprot_isoforms(self, 
                                     interpro_proteins: List[InterProProteinModel]) -> List[ProteinModel]:
        """
        Phase 3: Collect UniProt isoforms for InterPro proteins.
        
        Args:
            interpro_proteins: InterPro proteins to collect isoforms for
            
        Returns:
            List of collected protein isoforms
        """
        if self.progress.phase in ["uniprot_isoforms", "storage", "completed"]:
            self.logger.info("UniProt isoforms phase already completed, skipping")
            return []
        
        self.logger.info("Starting Phase 3: UniProt isoform collection")
        self.progress.phase = "uniprot_isoforms"
        self._save_progress()
        
        try:
            isoforms, report = await self.uniprot_collector.collect_isoforms_with_report(
                interpro_proteins, self.batch_size
            )
            
            self.progress.uniprot_isoforms_collected = len(isoforms)
            self.logger.info(
                "Phase 3 completed: collected %d UniProt isoforms",
                len(isoforms),
                extra={
                    "isoforms_collected": len(isoforms),
                    "successful_proteins": report.successful_proteins,
                    "failed_proteins": report.failed_proteins,
                    "success_rate": report.success_rate,
                    "duration_seconds": report.collection_duration
                }
            )
            
            self._save_progress()
            return isoforms
            
        except Exception as e:
            error_msg = f"Phase 3 failed: {str(e)}"
            self.progress.errors.append(error_msg)
            self.logger.error(error_msg, extra={"error_type": type(e).__name__})
            self._save_progress()
            raise
    
    async def collect_protein_isoforms(self, 
                                     interpro_proteins: List[InterProProteinModel]) -> List[ProteinModel]:
        """
        Alias for collect_uniprot_isoforms for backward compatibility.
        
        Args:
            interpro_proteins: InterPro proteins to collect isoforms for
            
        Returns:
            List of collected protein isoforms
        """
        return await self.collect_uniprot_isoforms(interpro_proteins)
    
    async def run_full_collection(self, 
                                page_size: int = 200,
                                store_data: bool = True) -> CollectionReport:
        """
        Execute the complete three-phase collection workflow.
        
        Args:
            page_size: Number of results per page for API queries
            store_data: Whether to store collected data in database
            
        Returns:
            Complete collection report
        """
        report = CollectionReport()
        report.progress = self.progress
        
        self.logger.info(
            "Starting complete TIM barrel data collection workflow",
            extra={
                "page_size": page_size,
                "store_data": store_data,
                "resume_from_phase": self.progress.phase
            }
        )
        
        try:
            # Phase 1: PFAM families
            pfam_families = await self.collect_pfam_families(page_size)
            report.pfam_families = pfam_families
            
            if not pfam_families and self.progress.pfam_families_collected == 0:
                error_msg = "No PFAM families found with TIM barrel annotations"
                report.api_errors.append(error_msg)
                self.logger.warning(error_msg)
                return report
            
            # Phase 2: InterPro proteins
            interpro_proteins = await self.collect_interpro_proteins(pfam_families, page_size)
            report.interpro_proteins = interpro_proteins
            
            if not interpro_proteins and self.progress.interpro_proteins_collected == 0:
                error_msg = "No human proteins found for PFAM families"
                report.api_errors.append(error_msg)
                self.logger.warning(error_msg)
                return report
            
            # Phase 3: UniProt isoforms
            uniprot_isoforms = await self.collect_uniprot_isoforms(interpro_proteins)
            report.uniprot_isoforms = uniprot_isoforms
            
            # Phase 4: Data storage (if requested)
            if store_data:
                await self.store_collected_data(
                    pfam_families, interpro_proteins, uniprot_isoforms, report
                )
            
            # Mark as completed
            self.progress.phase = "completed"
            self._save_progress()
            
            self.logger.info(
                "Complete collection workflow finished successfully",
                extra={
                    "pfam_families": len(report.pfam_families),
                    "interpro_proteins": len(report.interpro_proteins),
                    "uniprot_isoforms": len(report.uniprot_isoforms),
                    "success_rate": report.success_rate,
                    "duration_seconds": self.progress.duration_seconds
                }
            )
            
        except Exception as e:
            error_msg = f"Collection workflow failed: {str(e)}"
            report.api_errors.append(error_msg)
            self.logger.error(error_msg, extra={"error_type": type(e).__name__})
            raise
        
        return report
    
    async def store_collected_data(self,
                                 pfam_families: List[PfamFamilyModel],
                                 interpro_proteins: List[InterProProteinModel],
                                 uniprot_isoforms: List[ProteinModel],
                                 report: CollectionReport) -> None:
        """
        Store collected data in database with validation and duplicate prevention.
        
        Args:
            pfam_families: PFAM families to store
            interpro_proteins: InterPro proteins to store
            uniprot_isoforms: UniProt isoforms to store
            report: Collection report to update with storage results
        """
        if self.progress.phase in ["storage", "completed"]:
            self.logger.info("Data storage phase already completed, skipping")
            return
        
        self.logger.info("Starting data storage phase")
        self.progress.phase = "storage"
        self._save_progress()
        
        try:
            # Store PFAM families
            if pfam_families:
                pfam_result = self.storage.store_pfam_families(pfam_families)
                self.progress.pfam_families_stored = pfam_result.stats.successfully_stored
                if pfam_result.errors:
                    report.storage_errors.extend(pfam_result.errors)
                self.logger.info(f"Stored {pfam_result.stats.successfully_stored} PFAM families")
            
            # Store InterPro proteins
            if interpro_proteins:
                protein_result = self.storage.store_interpro_proteins(interpro_proteins)
                self.progress.interpro_proteins_stored = protein_result.stats.successfully_stored
                if protein_result.errors:
                    report.storage_errors.extend(protein_result.errors)
                self.logger.info(f"Stored {protein_result.stats.successfully_stored} InterPro proteins")
            
            # Store UniProt isoforms
            if uniprot_isoforms:
                isoform_result = self.storage.store_uniprot_isoforms(uniprot_isoforms)
                self.progress.uniprot_isoforms_stored = isoform_result.stats.successfully_stored
                if isoform_result.errors:
                    report.storage_errors.extend(isoform_result.errors)
                self.logger.info(f"Stored {isoform_result.stats.successfully_stored} UniProt isoforms")
            
            self.logger.info(
                "Data storage completed successfully",
                extra={
                    "pfam_families_stored": self.progress.pfam_families_stored,
                    "interpro_proteins_stored": self.progress.interpro_proteins_stored,
                    "uniprot_isoforms_stored": self.progress.uniprot_isoforms_stored,
                    "total_entities_stored": self.progress.total_entities_stored
                }
            )
            
        except Exception as e:
            error_msg = f"Data storage failed: {str(e)}"
            report.storage_errors.append(error_msg)
            self.logger.error(error_msg, extra={"error_type": type(e).__name__})
            raise
        finally:
            self._save_progress()
    
    def get_collection_status(self) -> Dict[str, Any]:
        """Get current collection status and progress."""
        return {
            "current_phase": self.progress.phase,
            "progress": self.progress.to_dict(),
            "batch_size": self.batch_size,
            "progress_file": str(self.progress_file) if self.progress_file else None,
            "database_connected": self.db_manager.test_connection()
        }
    
    def reset_collection(self) -> None:
        """Reset collection state for fresh start."""
        self.progress = CollectionProgress()
        if self.progress_file and self.progress_file.exists():
            self.progress_file.unlink()
        
        if self.uniprot_collector:
            self.uniprot_collector.reset_collection_state()
        
        self.logger.info("Collection state reset")


# Convenience functions for common operations
async def run_complete_collection(
    progress_file: Optional[str] = None,
    page_size: int = 200,
    batch_size: Optional[int] = None,
    store_data: bool = True
) -> CollectionReport:
    """
    Run complete TIM barrel data collection workflow.
    
    Args:
        progress_file: Path to progress file for resume functionality
        page_size: Number of results per page for API queries
        batch_size: Batch size for processing operations
        store_data: Whether to store collected data in database
        
    Returns:
        Complete collection report
    """
    async with DataCollector(progress_file, batch_size) as collector:
        return await collector.run_full_collection(page_size, store_data)


async def resume_collection(progress_file: str) -> CollectionReport:
    """
    Resume interrupted collection from progress file.
    
    Args:
        progress_file: Path to existing progress file
        
    Returns:
        Complete collection report
    """
    if not Path(progress_file).exists():
        raise FileNotFoundError(f"Progress file not found: {progress_file}")
    
    async with DataCollector(progress_file) as collector:
        return await collector.run_full_collection()