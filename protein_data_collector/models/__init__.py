"""
Data models and validation for the Protein Data Collector system.

This package provides Pydantic models for data validation and serialization
of protein data entities including PFAM families, InterPro proteins, and protein isoforms.
"""

from .entities import TIMBarrelEntryModel, InterProProteinModel, ProteinModel
from .validation import (
    ProteinSequenceValidator,
    TIMBarrelLocationValidator,
    DataValidator,
    ValidationResult,
    ValidationError
)

__all__ = [
    # Entity models
    "TIMBarrelEntryModel",
    "InterProProteinModel", 
    "ProteinModel",
    
    # Validation
    "ProteinSequenceValidator",
    "TIMBarrelLocationValidator",
    "DataValidator",
    "ValidationResult",
    "ValidationError",
]