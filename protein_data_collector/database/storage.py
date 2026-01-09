"""
Database storage operations with validation and duplicate prevention.

This module provides high-level database operations for storing collected protein data
with comprehensive validation, duplicate prevention, and batch processing capabilities.
"""

import logging
from typing import List, Dict, Any, Optional, Tuple, Set
from dataclasses import dataclass, field
from datetime import datetime
from sqlalchemy.orm import Session
from sqlalchemy.exc import IntegrityError, SQLAlchemyError
from sqlalchemy import and_, or_

from .connection import get_db_transaction, get_db_session
from .schema import TIMBarrelEntry, InterProProtein, Protein
from ..models.entities import TIMBarrelEntryModel, InterProProteinModel, ProteinModel
from ..errors import DataError, ValidationError


@dataclass
class StorageStats:
    """Statistics for database storage operations."""
    total_attempted: int = 0
    successfully_stored: int = 0
    duplicates_skipped: int = 0
    validation_errors: int = 0
    database_errors: int = 0
    start_time: Optional[datetime] = None
    end_time: Optional[datetime] = None
    
    @property
    def success_rate(self) -> float:
        """Calculate success rate as percentage."""
        if self.total_attempted == 0:
            return 0.0
        return (self.successfully_stored / self.total_attempted) * 100.0
    
    @property
    def duration_seconds(self) -> float:
        """Calculate operation duration in seconds."""
        if self.start_time and self.end_time:
            return (self.end_time - self.start_time).total_seconds()
        return 0.0


@dataclass
class StorageResult:
    """Result of database storage operation."""
    stats: StorageStats = field(default_factory=StorageStats)
    errors: List[str] = field(default_factory=list)
    stored_entities: List[str] = field(default_factory=list)  # IDs of successfully stored entities
    duplicate_entities: List[str] = field(default_factory=list)  # IDs of duplicate entities
    
    def add_error(self, error: str) -> None:
        """Add an error message to the result."""
        self.errors.append(error)
    
    def add_stored_entity(self, entity_id: str) -> None:
        """Add a successfully stored entity ID."""
        self.stored_entities.append(entity_id)
    
    def add_duplicate_entity(self, entity_id: str) -> None:
        """Add a duplicate entity ID."""
        self.duplicate_entities.append(entity_id)


class DatabaseStorage:
    """
    High-level database storage service with validation and duplicate prevention.
    
    Provides batch processing capabilities and comprehensive error handling
    for storing protein data entities.
    """
    
    def __init__(self, batch_size: int = 100):
        """
        Initialize database storage service.
        
        Args:
            batch_size: Number of entities to process in each batch
        """
        self.batch_size = batch_size
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
    
    def store_tim_barrel_entries(self, entries: List[TIMBarrelEntryModel]) -> StorageResult:
        """
        Store TIM barrel entries with duplicate prevention and validation.
        
        Args:
            entries: List of TIM barrel entry models to store
            
        Returns:
            StorageResult with statistics and error information
        """
        result = StorageResult()
        result.stats.start_time = datetime.now()
        result.stats.total_attempted = len(entries)
        
        if not entries:
            result.stats.end_time = datetime.now()
            return result
        
        self.logger.info(f"Starting storage of {len(entries)} TIM barrel entries")
        
        # Process in batches
        for i in range(0, len(entries), self.batch_size):
            batch = entries[i:i + self.batch_size]
            batch_result = self._store_tim_barrel_entries_batch(batch)
            
            # Merge batch results
            result.stats.successfully_stored += batch_result.stats.successfully_stored
            result.stats.duplicates_skipped += batch_result.stats.duplicates_skipped
            result.stats.validation_errors += batch_result.stats.validation_errors
            result.stats.database_errors += batch_result.stats.database_errors
            result.errors.extend(batch_result.errors)
            result.stored_entities.extend(batch_result.stored_entities)
            result.duplicate_entities.extend(batch_result.duplicate_entities)
        
        result.stats.end_time = datetime.now()
        
        self.logger.info(
            f"TIM barrel entry storage completed: {result.stats.successfully_stored} stored, "
            f"{result.stats.duplicates_skipped} duplicates, {len(result.errors)} errors"
        )
        
        return result
    
    def _store_tim_barrel_entries_batch(self, entries: List[TIMBarrelEntryModel]) -> StorageResult:
        """Store a batch of TIM barrel entries."""
        result = StorageResult()
        
        try:
            with get_db_transaction() as session:
                # Get existing accessions to check for duplicates
                existing_accessions = self._get_existing_tim_barrel_accessions(session, entries)
                
                for entry in entries:
                    try:
                        # Check for duplicates
                        if entry.accession in existing_accessions:
                            result.stats.duplicates_skipped += 1
                            result.add_duplicate_entity(entry.accession)
                            self.logger.debug(f"Skipping duplicate TIM barrel entry: {entry.accession}")
                            continue
                        
                        # Validate entry data
                        self._validate_tim_barrel_entry(entry)
                        
                        # Create database entity
                        db_entry = TIMBarrelEntry(
                            accession=entry.accession,
                            entry_type=entry.entry_type,
                            name=entry.name,
                            description=entry.description,
                            interpro_type=entry.interpro_type,
                            tim_barrel_annotation=entry.tim_barrel_annotation,
                            member_databases=entry.member_databases,
                            interpro_id=entry.interpro_id
                        )
                        
                        session.add(db_entry)
                        result.stats.successfully_stored += 1
                        result.add_stored_entity(entry.accession)
                        
                        # Add to existing set to prevent duplicates within batch
                        existing_accessions.add(entry.accession)
                        
                    except ValidationError as e:
                        result.stats.validation_errors += 1
                        error_msg = f"Validation error for TIM barrel entry {entry.accession}: {str(e)}"
                        result.add_error(error_msg)
                        self.logger.warning(error_msg)
                        
                    except Exception as e:
                        result.stats.database_errors += 1
                        error_msg = f"Database error for TIM barrel entry {entry.accession}: {str(e)}"
                        result.add_error(error_msg)
                        self.logger.error(error_msg)
                        
        except Exception as e:
            result.stats.database_errors += len(entries)
            error_msg = f"Batch storage failed for TIM barrel entries: {str(e)}"
            result.add_error(error_msg)
            self.logger.error(error_msg)
        
        return result
    
    def store_interpro_proteins(self, proteins: List[InterProProteinModel]) -> StorageResult:
        """
        Store InterPro proteins with duplicate prevention and validation.
        
        Args:
            proteins: List of InterPro protein models to store
            
        Returns:
            StorageResult with statistics and error information
        """
        result = StorageResult()
        result.stats.start_time = datetime.now()
        result.stats.total_attempted = len(proteins)
        
        if not proteins:
            result.stats.end_time = datetime.now()
            return result
        
        self.logger.info(f"Starting storage of {len(proteins)} InterPro proteins")
        
        # Process in batches
        for i in range(0, len(proteins), self.batch_size):
            batch = proteins[i:i + self.batch_size]
            batch_result = self._store_interpro_proteins_batch(batch)
            
            # Merge batch results
            result.stats.successfully_stored += batch_result.stats.successfully_stored
            result.stats.duplicates_skipped += batch_result.stats.duplicates_skipped
            result.stats.validation_errors += batch_result.stats.validation_errors
            result.stats.database_errors += batch_result.stats.database_errors
            result.errors.extend(batch_result.errors)
            result.stored_entities.extend(batch_result.stored_entities)
            result.duplicate_entities.extend(batch_result.duplicate_entities)
        
        result.stats.end_time = datetime.now()
        
        self.logger.info(
            f"InterPro protein storage completed: {result.stats.successfully_stored} stored, "
            f"{result.stats.duplicates_skipped} duplicates, {len(result.errors)} errors"
        )
        
        return result
    
    def _store_interpro_proteins_batch(self, proteins: List[InterProProteinModel]) -> StorageResult:
        """Store a batch of InterPro proteins."""
        result = StorageResult()
        
        try:
            with get_db_transaction() as session:
                # Get existing protein IDs and validate PFAM references
                existing_protein_ids = self._get_existing_protein_ids(session, proteins)
                valid_pfam_accessions = self._get_valid_pfam_accessions(session, proteins)
                
                for protein in proteins:
                    try:
                        # Check for duplicates
                        if protein.uniprot_id in existing_protein_ids:
                            result.stats.duplicates_skipped += 1
                            result.add_duplicate_entity(protein.uniprot_id)
                            self.logger.debug(f"Skipping duplicate InterPro protein: {protein.uniprot_id}")
                            continue
                        
                        # Validate protein data and references
                        self._validate_interpro_protein(protein, valid_pfam_accessions)
                        
                        # Create database entity
                        db_protein = InterProProtein(
                            uniprot_id=protein.uniprot_id,
                            pfam_accession=protein.pfam_accession,
                            name=protein.name,
                            organism=protein.organism
                        )
                        
                        session.add(db_protein)
                        result.stats.successfully_stored += 1
                        result.add_stored_entity(protein.uniprot_id)
                        
                        # Add to existing set to prevent duplicates within batch
                        existing_protein_ids.add(protein.uniprot_id)
                        
                    except ValidationError as e:
                        result.stats.validation_errors += 1
                        error_msg = f"Validation error for InterPro protein {protein.uniprot_id}: {str(e)}"
                        result.add_error(error_msg)
                        self.logger.warning(error_msg)
                        
                    except Exception as e:
                        result.stats.database_errors += 1
                        error_msg = f"Database error for InterPro protein {protein.uniprot_id}: {str(e)}"
                        result.add_error(error_msg)
                        self.logger.error(error_msg)
                        
        except Exception as e:
            result.stats.database_errors += len(proteins)
            error_msg = f"Batch storage failed for InterPro proteins: {str(e)}"
            result.add_error(error_msg)
            self.logger.error(error_msg)
        
        return result
    
    def store_protein_isoforms(self, isoforms: List[ProteinModel]) -> StorageResult:
        """
        Store protein isoforms with duplicate prevention and validation.
        
        Args:
            isoforms: List of protein isoform models to store
            
        Returns:
            StorageResult with statistics and error information
        """
        result = StorageResult()
        result.stats.start_time = datetime.now()
        result.stats.total_attempted = len(isoforms)
        
        if not isoforms:
            result.stats.end_time = datetime.now()
            return result
        
        self.logger.info(f"Starting storage of {len(isoforms)} protein isoforms")
        
        # Process in batches
        for i in range(0, len(isoforms), self.batch_size):
            batch = isoforms[i:i + self.batch_size]
            batch_result = self._store_protein_isoforms_batch(batch)
            
            # Merge batch results
            result.stats.successfully_stored += batch_result.stats.successfully_stored
            result.stats.duplicates_skipped += batch_result.stats.duplicates_skipped
            result.stats.validation_errors += batch_result.stats.validation_errors
            result.stats.database_errors += batch_result.stats.database_errors
            result.errors.extend(batch_result.errors)
            result.stored_entities.extend(batch_result.stored_entities)
            result.duplicate_entities.extend(batch_result.duplicate_entities)
        
        result.stats.end_time = datetime.now()
        
        self.logger.info(
            f"Protein isoform storage completed: {result.stats.successfully_stored} stored, "
            f"{result.stats.duplicates_skipped} duplicates, {len(result.errors)} errors"
        )
        
        return result
    
    def _store_protein_isoforms_batch(self, isoforms: List[ProteinModel]) -> StorageResult:
        """Store a batch of protein isoforms."""
        result = StorageResult()
        
        try:
            with get_db_transaction() as session:
                # Get existing isoform IDs and validate parent protein references
                existing_isoform_ids = self._get_existing_isoform_ids(session, isoforms)
                valid_parent_protein_ids = self._get_valid_parent_protein_ids(session, isoforms)
                
                for isoform in isoforms:
                    try:
                        # Check for duplicates
                        if isoform.isoform_id in existing_isoform_ids:
                            result.stats.duplicates_skipped += 1
                            result.add_duplicate_entity(isoform.isoform_id)
                            self.logger.debug(f"Skipping duplicate protein isoform: {isoform.isoform_id}")
                            continue
                        
                        # Validate isoform data and references
                        self._validate_protein_isoform(isoform, valid_parent_protein_ids)
                        
                        # Create database entity
                        db_isoform = Protein(
                            isoform_id=isoform.isoform_id,
                            parent_protein_id=isoform.parent_protein_id,
                            sequence=isoform.sequence,
                            sequence_length=isoform.sequence_length,
                            exon_annotations=isoform.exon_annotations,
                            exon_count=isoform.exon_count,
                            tim_barrel_location=isoform.tim_barrel_location,
                            organism=isoform.organism,
                            name=isoform.name,
                            description=isoform.description
                        )
                        
                        session.add(db_isoform)
                        result.stats.successfully_stored += 1
                        result.add_stored_entity(isoform.isoform_id)
                        
                        # Add to existing set to prevent duplicates within batch
                        existing_isoform_ids.add(isoform.isoform_id)
                        
                    except ValidationError as e:
                        result.stats.validation_errors += 1
                        error_msg = f"Validation error for protein isoform {isoform.isoform_id}: {str(e)}"
                        result.add_error(error_msg)
                        self.logger.warning(error_msg)
                        
                    except Exception as e:
                        result.stats.database_errors += 1
                        error_msg = f"Database error for protein isoform {isoform.isoform_id}: {str(e)}"
                        result.add_error(error_msg)
                        self.logger.error(error_msg)
                        
        except Exception as e:
            result.stats.database_errors += len(isoforms)
            error_msg = f"Batch storage failed for protein isoforms: {str(e)}"
            result.add_error(error_msg)
            self.logger.error(error_msg)
        
        return result
    
    def _get_existing_tim_barrel_accessions(self, session: Session, entries: List[TIMBarrelEntryModel]) -> Set[str]:
        """Get existing TIM barrel entry accessions from database."""
        accessions = [entry.accession for entry in entries]
        existing = session.query(TIMBarrelEntry.accession).filter(TIMBarrelEntry.accession.in_(accessions)).all()
        return {acc[0] for acc in existing}
    
    def _get_existing_protein_ids(self, session: Session, proteins: List[InterProProteinModel]) -> Set[str]:
        """Get existing protein IDs from database."""
        protein_ids = [protein.uniprot_id for protein in proteins]
        existing = session.query(InterProProtein.uniprot_id).filter(InterProProtein.uniprot_id.in_(protein_ids)).all()
        return {pid[0] for pid in existing}
    
    def _get_existing_isoform_ids(self, session: Session, isoforms: List[ProteinModel]) -> Set[str]:
        """Get existing isoform IDs from database."""
        isoform_ids = [isoform.isoform_id for isoform in isoforms]
        existing = session.query(Protein.isoform_id).filter(Protein.isoform_id.in_(isoform_ids)).all()
        return {iid[0] for iid in existing}
    
    def _get_valid_pfam_accessions(self, session: Session, proteins: List[InterProProteinModel]) -> Set[str]:
        """Get valid TIM barrel entry accessions that exist in database."""
        tim_barrel_accessions = list(set(protein.tim_barrel_accession for protein in proteins))
        existing = session.query(TIMBarrelEntry.accession).filter(TIMBarrelEntry.accession.in_(tim_barrel_accessions)).all()
        return {acc[0] for acc in existing}
    
    def _get_valid_parent_protein_ids(self, session: Session, isoforms: List[ProteinModel]) -> Set[str]:
        """Get valid parent protein IDs that exist in database."""
        parent_ids = list(set(isoform.parent_protein_id for isoform in isoforms))
        existing = session.query(InterProProtein.uniprot_id).filter(InterProProtein.uniprot_id.in_(parent_ids)).all()
        return {pid[0] for pid in existing}
    
    def _validate_tim_barrel_entry(self, entry: TIMBarrelEntryModel) -> None:
        """Validate TIM barrel entry data."""
        if not entry.accession or not entry.accession.strip():
            raise ValidationError("TIM barrel entry accession cannot be empty")
        
        if not entry.name or not entry.name.strip():
            raise ValidationError("TIM barrel entry name cannot be empty")
        
        if not entry.tim_barrel_annotation or not entry.tim_barrel_annotation.strip():
            raise ValidationError("TIM barrel annotation cannot be empty")
        
        if entry.entry_type not in ['pfam', 'interpro']:
            raise ValidationError("Entry type must be 'pfam' or 'interpro'")
    
    def _validate_interpro_protein(self, protein: InterProProteinModel, valid_tim_barrel_accessions: Set[str]) -> None:
        """Validate InterPro protein data and references."""
        if not protein.uniprot_id or not protein.uniprot_id.strip():
            raise ValidationError("UniProt ID cannot be empty")
        
        if not protein.tim_barrel_accession or not protein.tim_barrel_accession.strip():
            raise ValidationError("TIM barrel accession cannot be empty")
        
        if protein.tim_barrel_accession not in valid_tim_barrel_accessions:
            raise ValidationError(f"TIM barrel accession {protein.tim_barrel_accession} does not exist in database")
        
        if not protein.organism or not protein.organism.strip():
            raise ValidationError("Organism cannot be empty")
        
        # Validate organism is human
        if protein.organism.lower() not in ["homo sapiens", "human"]:
            raise ValidationError(f"Expected human organism, got: {protein.organism}")
    
    def _validate_protein_isoform(self, isoform: ProteinModel, valid_parent_protein_ids: Set[str]) -> None:
        """Validate protein isoform data and references."""
        if not isoform.isoform_id or not isoform.isoform_id.strip():
            raise ValidationError("Isoform ID cannot be empty")
        
        if not isoform.parent_protein_id or not isoform.parent_protein_id.strip():
            raise ValidationError("Parent protein ID cannot be empty")
        
        if isoform.parent_protein_id not in valid_parent_protein_ids:
            raise ValidationError(f"Parent protein ID {isoform.parent_protein_id} does not exist in database")
        
        if not isoform.sequence or not isoform.sequence.strip():
            raise ValidationError("Protein sequence cannot be empty")
        
        if isoform.sequence_length <= 0:
            raise ValidationError("Sequence length must be positive")
        
        if len(isoform.sequence) != isoform.sequence_length:
            raise ValidationError(f"Sequence length mismatch: expected {isoform.sequence_length}, got {len(isoform.sequence)}")
        
        # Validate amino acid sequence
        valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
        invalid_chars = set(isoform.sequence.upper()) - valid_amino_acids
        if invalid_chars:
            raise ValidationError(f"Invalid amino acid characters in sequence: {invalid_chars}")
        
        # Validate TIM barrel location if present
        if isoform.tim_barrel_location:
            self._validate_tim_barrel_location(isoform.tim_barrel_location, isoform.sequence_length)
        
        # Validate exon count consistency
        if isoform.exon_count is not None and isoform.exon_annotations:
            if isinstance(isoform.exon_annotations, dict) and "exons" in isoform.exon_annotations:
                actual_exon_count = len(isoform.exon_annotations["exons"])
                if isoform.exon_count != actual_exon_count:
                    raise ValidationError(f"Exon count mismatch: expected {isoform.exon_count}, got {actual_exon_count}")
    
    def _validate_tim_barrel_location(self, location: Dict[str, Any], sequence_length: int) -> None:
        """Validate TIM barrel location coordinates."""
        if not isinstance(location, dict):
            raise ValidationError("TIM barrel location must be a dictionary")
        
        start = location.get("start")
        end = location.get("end")
        
        if start is not None and end is not None:
            if not isinstance(start, int) or not isinstance(end, int):
                raise ValidationError("TIM barrel start and end must be integers")
            
            if start < 1 or end < 1:
                raise ValidationError("TIM barrel coordinates must be positive")
            
            if start >= end:
                raise ValidationError("TIM barrel start must be less than end")
            
            if start > sequence_length or end > sequence_length:
                raise ValidationError(f"TIM barrel location ({start}-{end}) exceeds sequence length ({sequence_length})")
    
    def get_storage_statistics(self) -> Dict[str, Any]:
        """Get current storage statistics from database."""
        try:
            with get_db_session() as session:
                tim_barrel_count = session.query(TIMBarrelEntry).count()
                protein_count = session.query(InterProProtein).count()
                isoform_count = session.query(Protein).count()
                
                return {
                    "tim_barrel_entries": tim_barrel_count,
                    "interpro_proteins": protein_count,
                    "protein_isoforms": isoform_count,
                    "total_entities": tim_barrel_count + protein_count + isoform_count
                }
        except Exception as e:
            self.logger.error(f"Failed to get storage statistics: {str(e)}")
            return {
                "tim_barrel_entries": 0,
                "interpro_proteins": 0,
                "protein_isoforms": 0,
                "total_entities": 0,
                "error": str(e)
            }


# Convenience functions for common operations
def store_all_entities(
    tim_barrel_entries: List[TIMBarrelEntryModel],
    interpro_proteins: List[InterProProteinModel],
    protein_isoforms: List[ProteinModel],
    batch_size: int = 100
) -> Tuple[StorageResult, StorageResult, StorageResult]:
    """
    Store all entity types in the correct order with proper dependencies.
    
    Args:
        tim_barrel_entries: TIM barrel entries to store
        interpro_proteins: InterPro proteins to store
        protein_isoforms: Protein isoforms to store
        batch_size: Batch size for processing
        
    Returns:
        Tuple of (pfam_result, protein_result, isoform_result)
    """
    storage = DatabaseStorage(batch_size)
    
    # Store in dependency order: TIM barrel entries → proteins → isoforms
    tim_barrel_result = storage.store_tim_barrel_entries(tim_barrel_entries)
    protein_result = storage.store_interpro_proteins(interpro_proteins)
    isoform_result = storage.store_protein_isoforms(protein_isoforms)
    
    return pfam_result, protein_result, isoform_result


def get_database_statistics() -> Dict[str, Any]:
    """Get current database statistics."""
    storage = DatabaseStorage()
    return storage.get_storage_statistics()