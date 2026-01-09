"""
Data validation utilities for protein data.

This module provides validators for protein sequences, TIM barrel coordinates,
and other biological data with comprehensive error reporting.
"""

import re
from typing import Dict, Any, List, Optional, Union
from dataclasses import dataclass
from enum import Enum


class ValidationErrorType(Enum):
    """Types of validation errors."""
    INVALID_SEQUENCE = "invalid_sequence"
    INVALID_COORDINATES = "invalid_coordinates"
    MISSING_REQUIRED_FIELD = "missing_required_field"
    INVALID_FORMAT = "invalid_format"
    OUT_OF_BOUNDS = "out_of_bounds"
    INCONSISTENT_DATA = "inconsistent_data"


@dataclass
class ValidationError:
    """Represents a validation error with context."""
    error_type: ValidationErrorType
    field_name: str
    message: str
    value: Any = None
    context: Optional[Dict[str, Any]] = None


@dataclass
class ValidationResult:
    """Result of a validation operation."""
    is_valid: bool
    errors: List[ValidationError]
    warnings: List[str] = None
    
    @property
    def error_message(self) -> str:
        """Get formatted error message."""
        if not self.errors:
            return ""
        return "; ".join([error.message for error in self.errors])
    
    def add_error(self, error: ValidationError) -> None:
        """Add a validation error."""
        self.errors.append(error)
        self.is_valid = False
    
    def add_warning(self, warning: str) -> None:
        """Add a validation warning."""
        if self.warnings is None:
            self.warnings = []
        self.warnings.append(warning)


class ProteinSequenceValidator:
    """Validator for protein amino acid sequences."""
    
    # Standard 20 amino acids
    VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
    
    # Extended amino acids including ambiguous codes
    EXTENDED_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWYXBZJUO*")
    
    def __init__(self, allow_extended: bool = False):
        """
        Initialize validator.
        
        Args:
            allow_extended: Whether to allow extended amino acid codes (X, B, Z, etc.)
        """
        self.allow_extended = allow_extended
        self.valid_chars = self.EXTENDED_AMINO_ACIDS if allow_extended else self.VALID_AMINO_ACIDS
    
    def validate(self, sequence: str) -> ValidationResult:
        """
        Validate a protein sequence.
        
        Args:
            sequence: Protein sequence to validate
            
        Returns:
            ValidationResult with validation status and errors
        """
        result = ValidationResult(is_valid=True, errors=[])
        
        if not sequence:
            result.add_error(ValidationError(
                error_type=ValidationErrorType.MISSING_REQUIRED_FIELD,
                field_name="sequence",
                message="Protein sequence cannot be empty",
                value=sequence
            ))
            return result
        
        # Clean sequence (remove whitespace)
        clean_sequence = sequence.strip().upper()
        
        if not clean_sequence:
            result.add_error(ValidationError(
                error_type=ValidationErrorType.INVALID_SEQUENCE,
                field_name="sequence",
                message="Protein sequence cannot be whitespace-only",
                value=sequence
            ))
            return result
        
        # Check for invalid characters
        invalid_chars = set(clean_sequence) - self.valid_chars
        if invalid_chars:
            result.add_error(ValidationError(
                error_type=ValidationErrorType.INVALID_SEQUENCE,
                field_name="sequence",
                message=f"Invalid amino acid characters found: {', '.join(sorted(invalid_chars))}",
                value=sequence,
                context={"invalid_characters": list(invalid_chars)}
            ))
        
        # Check sequence length
        if len(clean_sequence) < 10:
            result.add_warning("Protein sequence is unusually short (< 10 amino acids)")
        elif len(clean_sequence) > 50000:
            result.add_warning("Protein sequence is unusually long (> 50,000 amino acids)")
        
        # Check for unusual patterns
        if '*' in clean_sequence and clean_sequence.index('*') < len(clean_sequence) - 1:
            result.add_warning("Stop codon (*) found in middle of sequence")
        
        return result
    
    def clean_sequence(self, sequence: str) -> str:
        """
        Clean and normalize a protein sequence.
        
        Args:
            sequence: Raw protein sequence
            
        Returns:
            Cleaned sequence in uppercase without whitespace
        """
        return sequence.strip().upper()


class TIMBarrelLocationValidator:
    """Validator for TIM barrel location coordinates."""
    
    def validate(self, location: Dict[str, Any], sequence_length: int) -> ValidationResult:
        """
        Validate TIM barrel location coordinates.
        
        Args:
            location: Dictionary containing location information
            sequence_length: Length of the protein sequence
            
        Returns:
            ValidationResult with validation status and errors
        """
        result = ValidationResult(is_valid=True, errors=[])
        
        if not isinstance(location, dict):
            result.add_error(ValidationError(
                error_type=ValidationErrorType.INVALID_FORMAT,
                field_name="tim_barrel_location",
                message="TIM barrel location must be a dictionary",
                value=location
            ))
            return result
        
        # Check for required fields
        required_fields = ['start', 'end']
        for field in required_fields:
            if field not in location:
                result.add_error(ValidationError(
                    error_type=ValidationErrorType.MISSING_REQUIRED_FIELD,
                    field_name=f"tim_barrel_location.{field}",
                    message=f"Missing required field: {field}",
                    value=location
                ))
        
        if result.errors:
            return result
        
        start = location.get('start')
        end = location.get('end')
        
        # Validate coordinate types
        if not isinstance(start, int) or not isinstance(end, int):
            result.add_error(ValidationError(
                error_type=ValidationErrorType.INVALID_FORMAT,
                field_name="tim_barrel_location",
                message="Start and end coordinates must be integers",
                value=location
            ))
            return result
        
        # Validate coordinate values
        if start < 1:
            result.add_error(ValidationError(
                error_type=ValidationErrorType.OUT_OF_BOUNDS,
                field_name="tim_barrel_location.start",
                message="Start coordinate must be >= 1",
                value=start
            ))
        
        if end > sequence_length:
            result.add_error(ValidationError(
                error_type=ValidationErrorType.OUT_OF_BOUNDS,
                field_name="tim_barrel_location.end",
                message=f"End coordinate {end} exceeds sequence length {sequence_length}",
                value=end,
                context={"sequence_length": sequence_length}
            ))
        
        if start >= end:
            result.add_error(ValidationError(
                error_type=ValidationErrorType.INCONSISTENT_DATA,
                field_name="tim_barrel_location",
                message=f"Start coordinate {start} must be less than end coordinate {end}",
                value=location
            ))
        
        # Validate optional confidence score
        if 'confidence' in location:
            confidence = location['confidence']
            if not isinstance(confidence, (int, float)) or not (0.0 <= confidence <= 1.0):
                result.add_error(ValidationError(
                    error_type=ValidationErrorType.INVALID_FORMAT,
                    field_name="tim_barrel_location.confidence",
                    message="Confidence score must be a number between 0.0 and 1.0",
                    value=confidence
                ))
        
        # Check for reasonable TIM barrel length
        if result.is_valid:
            barrel_length = end - start + 1
            if barrel_length < 200:
                result.add_warning(f"TIM barrel region is unusually short ({barrel_length} residues)")
            elif barrel_length > 400:
                result.add_warning(f"TIM barrel region is unusually long ({barrel_length} residues)")
        
        return result


class DataValidator:
    """General data validator for protein entities."""
    
    def __init__(self):
        """Initialize validator with sub-validators."""
        self.sequence_validator = ProteinSequenceValidator()
        self.location_validator = TIMBarrelLocationValidator()
    
    def validate_pfam_family(self, data: Dict[str, Any]) -> ValidationResult:
        """Validate PFAM family data."""
        result = ValidationResult(is_valid=True, errors=[])
        
        # Check required fields
        required_fields = ['accession', 'name', 'tim_barrel_annotation']
        for field in required_fields:
            if field not in data or not data[field]:
                result.add_error(ValidationError(
                    error_type=ValidationErrorType.MISSING_REQUIRED_FIELD,
                    field_name=field,
                    message=f"Missing required field: {field}",
                    value=data.get(field)
                ))
        
        # Validate accession format
        if 'accession' in data and data['accession']:
            accession = data['accession'].strip()
            if not re.match(r'^[A-Z0-9]+$', accession):
                result.add_error(ValidationError(
                    error_type=ValidationErrorType.INVALID_FORMAT,
                    field_name="accession",
                    message="PFAM accession must contain only uppercase letters and numbers",
                    value=accession
                ))
        
        return result
    
    def validate_interpro_protein(self, data: Dict[str, Any]) -> ValidationResult:
        """Validate InterPro protein data."""
        result = ValidationResult(is_valid=True, errors=[])
        
        # Check required fields
        required_fields = ['uniprot_id', 'pfam_accession']
        for field in required_fields:
            if field not in data or not data[field]:
                result.add_error(ValidationError(
                    error_type=ValidationErrorType.MISSING_REQUIRED_FIELD,
                    field_name=field,
                    message=f"Missing required field: {field}",
                    value=data.get(field)
                ))
        
        # Validate UniProt ID format
        if 'uniprot_id' in data and data['uniprot_id']:
            uniprot_id = data['uniprot_id'].strip()
            if not re.match(r'^[A-Z0-9]+$', uniprot_id):
                result.add_error(ValidationError(
                    error_type=ValidationErrorType.INVALID_FORMAT,
                    field_name="uniprot_id",
                    message="UniProt ID must contain only uppercase letters and numbers",
                    value=uniprot_id
                ))
        
        return result
    
    def validate_protein_isoform(self, data: Dict[str, Any]) -> ValidationResult:
        """Validate protein isoform data."""
        result = ValidationResult(is_valid=True, errors=[])
        
        # Check required fields
        required_fields = ['isoform_id', 'parent_protein_id', 'sequence', 'sequence_length']
        for field in required_fields:
            if field not in data or data[field] is None:
                result.add_error(ValidationError(
                    error_type=ValidationErrorType.MISSING_REQUIRED_FIELD,
                    field_name=field,
                    message=f"Missing required field: {field}",
                    value=data.get(field)
                ))
        
        if result.errors:
            return result
        
        # Validate sequence
        sequence_result = self.sequence_validator.validate(data['sequence'])
        result.errors.extend(sequence_result.errors)
        if sequence_result.warnings:
            if result.warnings is None:
                result.warnings = []
            result.warnings.extend(sequence_result.warnings)
        
        # Validate sequence length consistency
        actual_length = len(data['sequence'].strip())
        declared_length = data['sequence_length']
        if actual_length != declared_length:
            result.add_error(ValidationError(
                error_type=ValidationErrorType.INCONSISTENT_DATA,
                field_name="sequence_length",
                message=f"Declared sequence length {declared_length} does not match actual length {actual_length}",
                value=declared_length,
                context={"actual_length": actual_length}
            ))
        
        # Validate TIM barrel location if present
        if 'tim_barrel_location' in data and data['tim_barrel_location']:
            location_result = self.location_validator.validate(
                data['tim_barrel_location'], 
                data['sequence_length']
            )
            result.errors.extend(location_result.errors)
            if location_result.warnings:
                if result.warnings is None:
                    result.warnings = []
                result.warnings.extend(location_result.warnings)
        
        # Update validity based on errors
        result.is_valid = len(result.errors) == 0
        
        return result