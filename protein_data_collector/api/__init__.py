"""
API integration layer for external database access.

This package contains clients for InterPro and UniProt REST API integration,
along with retry logic and comprehensive rate limiting functionality.
"""

from .interpro_client import (
    InterProAPIClient,
    get_tim_barrel_entries,
    get_human_proteins_for_tim_barrel_entries
)

from .uniprot_client import (
    UnifiedUniProtClient,
    UniProtRESTClient
)

__all__ = [
    'InterProAPIClient',
    'UnifiedUniProtClient',
    'UniProtRESTClient',
    'get_tim_barrel_entries',
    'get_human_proteins_for_tim_barrel_entries'
]