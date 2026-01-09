"""
API integration layer for external database access.

This package contains clients for InterPro, UniProt, and MCP server integration,
along with retry logic and comprehensive rate limiting functionality.
"""

from .interpro_client import (
    InterProAPIClient,
    get_tim_barrel_entries,
    get_human_proteins_for_tim_barrel_entries
)

from .uniprot_client import (
    UnifiedUniProtClient,
    get_protein_with_isoforms,
    get_proteins_batch
)

__all__ = [
    'InterProAPIClient',
    'UnifiedUniProtClient',
    'get_tim_barrel_entries',
    'get_human_proteins_for_tim_barrel_entries',
    'get_protein_with_isoforms',
    'get_proteins_batch'
]