"""
Data collection orchestration and workflow management.

This package contains the main data collector service that coordinates
the three-phase collection process: PFAM families → proteins → isoforms.
"""

from .interpro_collector import (
    InterProCollector,
    CollectionResult,
    CollectionStats,
    collect_tim_barrel_families_and_proteins,
    validate_organism_filtering
)
from .data_collector import (
    DataCollector,
    CollectionProgress,
    CollectionReport,
    run_complete_collection,
    resume_collection
)

__all__ = [
    'InterProCollector',
    'CollectionResult', 
    'CollectionStats',
    'collect_tim_barrel_families_and_proteins',
    'validate_organism_filtering',
    'DataCollector',
    'CollectionProgress',
    'CollectionReport',
    'run_complete_collection',
    'resume_collection'
]