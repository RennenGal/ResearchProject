"""
UniProt REST API client for protein data collection.

This module provides a simplified interface for accessing UniProt data through
direct REST API calls with comprehensive rate limiting and retry logic.
"""

import asyncio
import json
import logging
from typing import Dict, List, Optional, Any, Union
from urllib.parse import urljoin, urlencode
from datetime import datetime
import httpx

from ..config import APIConfig, get_config
from ..retry import RetryController, get_retry_controller
from ..rate_limiter import RateLimitConfig, get_rate_limit_manager
from ..cache import get_global_cache, CachedAPIClient
from ..models.entities import ProteinModel
from ..errors import APIError, NetworkError, DataError, ConfigurationError, ErrorContext, create_error_context


class UniProtRESTClient:
    """Client for UniProt REST API with rate limiting and caching."""
    
    def __init__(self, api_config: APIConfig, retry_controller: RetryController):
        """
        Initialize UniProt REST API client.
        
        Args:
            api_config: API configuration
            retry_controller: Retry controller for operations
        """
        self.config = api_config
        self.retry_controller = retry_controller
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        
        # Set up rate limiting for UniProt REST API
        system_config = get_config()
        rate_limit_config = RateLimitConfig(
            requests_per_second=system_config.rate_limiting.uniprot_requests_per_second,
            burst_limit=system_config.rate_limiting.uniprot_burst_limit,
            burst_window_seconds=system_config.rate_limiting.uniprot_burst_window_seconds,
            violation_initial_delay=system_config.rate_limiting.violation_initial_delay,
            violation_backoff_multiplier=system_config.rate_limiting.violation_backoff_multiplier,
            violation_max_delay=system_config.rate_limiting.violation_max_delay,
            soft_limit_threshold=system_config.rate_limiting.soft_limit_threshold,
            enable_monitoring=system_config.rate_limiting.enable_monitoring,
            enable_reporting=system_config.rate_limiting.enable_reporting
        )
        
        # Create rate limiter for UniProt REST API
        rate_limit_manager = get_rate_limit_manager()
        self.rate_limiter = rate_limit_manager.create_limiter("UniProt_REST", rate_limit_config)
        
        # Set up HTTP client
        self.client = httpx.AsyncClient(
            timeout=httpx.Timeout(self.config.request_timeout),
            limits=httpx.Limits(max_connections=100, max_keepalive_connections=20)
        )
    
    async def __aenter__(self):
        """Async context manager entry."""
        return self
    
    async def __aexit__(self, exc_type, exc_val, exc_tb):
        """Async context manager exit."""
        await self.client.aclose()
    
    async def _make_request(self, endpoint: str, params: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        Make rate-limited request to UniProt REST API.
        
        Args:
            endpoint: API endpoint
            params: Query parameters
            
        Returns:
            API response data
            
        Raises:
            APIError: If API request fails
            NetworkError: If network request fails
            DataError: If response data is invalid
        """
        url = urljoin(self.config.uniprot_base_url, endpoint)
        
        # Apply rate limiting
        await self.rate_limiter.acquire()
        
        try:
            # Make request with retry logic
            async def make_request():
                response = await self.client.get(url, params=params or {})
                response.raise_for_status()
                return response.json()
            
            result = await self.retry_controller.execute_with_retry_async(
                make_request,
                database="UniProt_REST",
                operation_name=f"request_{endpoint}"
            )
            
            return result
            
        except httpx.HTTPStatusError as e:
            if e.response.status_code == 404:
                return {}  # Not found is not an error for our use case
            elif e.response.status_code == 429:
                raise APIError(
                    f"Rate limit exceeded for UniProt REST API: {e.response.status_code}",
                    status_code=e.response.status_code,
                    context=create_error_context(
                        operation="uniprot_rest_request",
                        endpoint=endpoint,
                        params=params,
                        response_status=e.response.status_code
                    )
                )
            else:
                raise APIError(
                    f"UniProt REST API error: {e.response.status_code} - {e.response.text}",
                    status_code=e.response.status_code,
                    context=create_error_context(
                        operation="uniprot_rest_request",
                        endpoint=endpoint,
                        params=params,
                        response_status=e.response.status_code
                    )
                )
                
        except httpx.RequestError as e:
            raise NetworkError(
                f"Network error accessing UniProt REST API: {str(e)}",
                context=create_error_context(
                    operation="uniprot_rest_request",
                    endpoint=endpoint,
                    params=params,
                    error=str(e)
                )
            )
            
        except json.JSONDecodeError as e:
            raise DataError(f"Invalid JSON response from UniProt REST API: {str(e)}")
        except Exception as e:
            if isinstance(e, (APIError, DataError)):
                raise
            raise APIError(f"Unexpected error in UniProt REST API request: {str(e)}")
    
    async def get_protein_data(self, uniprot_id: str) -> Optional[Dict[str, Any]]:
        """
        Get protein data from UniProt REST API.
        
        Args:
            uniprot_id: UniProt accession ID
            
        Returns:
            Protein data dictionary or None if not found
        """
        try:
            self.logger.debug(f"Fetching protein data for {uniprot_id} from UniProt REST API")
            
            # Request comprehensive protein data
            endpoint = f"uniprotkb/{uniprot_id}"
            params = {
                "format": "json",
                "fields": ",".join([
                    "accession", "id", "protein_name", "organism_name", "organism_id",
                    "sequence", "length", "mass", "reviewed", "protein_existence",
                    "annotation_score", "keywords", "go", "cc_alternative_products",
                    "ft_domain", "ft_region", "ft_site", "xref_interpro", "xref_pfam",
                    "xref_smart", "xref_cdd", "xref_ensembl", "xref_refseq", "xref_embl",
                    "xref_pdb", "cc_catalytic_activity", "cc_pathway", "cc_function"
                ])
            }
            
            result = await self._make_request(endpoint, params)
            
            if not result:
                self.logger.warning(f"No data found for protein {uniprot_id}")
                return None
            
            self.logger.debug(f"Successfully retrieved protein data for {uniprot_id}")
            return result
            
        except Exception as e:
            self.logger.error(f"Failed to get protein data for {uniprot_id}: {str(e)}")
            raise
    
    async def get_protein_isoforms(self, uniprot_id: str) -> List[Dict[str, Any]]:
        """
        Get protein isoform data from UniProt REST API.
        
        Args:
            uniprot_id: UniProt accession ID
            
        Returns:
            List of isoform data dictionaries
        """
        try:
            self.logger.debug(f"Fetching isoform data for {uniprot_id} from UniProt REST API")
            
            # Get protein data which includes isoform information
            protein_data = await self.get_protein_data(uniprot_id)
            
            if not protein_data:
                return []
            
            # Extract isoform information from alternative products
            isoforms = []
            
            # Add canonical isoform
            if "sequence" in protein_data:
                canonical_isoform = {
                    "id": f"{uniprot_id}-1",
                    "name": "Canonical",
                    "sequence": protein_data["sequence"]["value"],
                    "length": protein_data["sequence"]["length"],
                    "is_canonical": True,
                    "description": "Canonical sequence"
                }
                isoforms.append(canonical_isoform)
            
            # Extract alternative isoforms from comments
            comments = protein_data.get("comments", [])
            for comment in comments:
                if comment.get("commentType") == "ALTERNATIVE PRODUCTS":
                    alt_isoforms = comment.get("isoforms", [])
                    for isoform in alt_isoforms:
                        isoform_ids = isoform.get("isoformIds", [])
                        if isoform_ids:
                            isoform_id = isoform_ids[0]
                            # Skip canonical (already added)
                            if isoform_id.endswith("-1"):
                                continue
                            
                            isoform_data = {
                                "id": isoform_id,
                                "name": isoform.get("name", {}).get("value", "Alternative isoform"),
                                "sequence": "",  # Would need separate API call
                                "length": 0,     # Would need separate API call
                                "is_canonical": False,
                                "description": f"Alternative isoform - {isoform.get('name', {}).get('value', 'unnamed')}"
                            }
                            isoforms.append(isoform_data)
            
            self.logger.debug(f"Found {len(isoforms)} isoforms for {uniprot_id}")
            return isoforms
            
        except Exception as e:
            self.logger.error(f"Failed to get isoform data for {uniprot_id}: {str(e)}")
            raise
    
    async def search_proteins(self, query: str, limit: int = 100) -> List[Dict[str, Any]]:
        """
        Search for proteins using UniProt REST API.
        
        Args:
            query: Search query
            limit: Maximum number of results
            
        Returns:
            List of protein data dictionaries
        """
        try:
            self.logger.debug(f"Searching proteins with query: {query}")
            
            endpoint = "uniprotkb/search"
            params = {
                "query": query,
                "format": "json",
                "size": min(limit, 500),  # UniProt API limit
                "fields": "accession,id,protein_name,organism_name,reviewed,length"
            }
            
            result = await self._make_request(endpoint, params)
            
            proteins = result.get("results", [])
            self.logger.debug(f"Found {len(proteins)} proteins for query: {query}")
            
            return proteins
            
        except Exception as e:
            self.logger.error(f"Failed to search proteins with query '{query}': {str(e)}")
            raise
    
    async def get_status(self) -> Dict[str, Any]:
        """
        Get client status information.
        
        Returns:
            Dictionary with status information
        """
        # Test API availability
        try:
            test_result = await self._make_request("uniprotkb/P04637", {"format": "json", "fields": "accession"})
            api_available = test_result is not None
        except Exception:
            api_available = False
        
        # Get rate limiting statistics
        rate_limiting_stats = self.rate_limiter.get_stats()
        
        return {
            "api_available": api_available,
            "base_url": self.config.uniprot_base_url,
            "timeout": self.config.request_timeout,
            "rate_limiting": rate_limiting_stats,
            "client_type": "rest_only"
        }


class UnifiedUniProtClient:
    """
    Simplified UniProt client that uses only REST API.
    
    This class provides a unified interface for UniProt data access,
    using only the REST API for simplicity and reliability.
    """
    
    def __init__(self, api_config: Optional[APIConfig] = None, retry_controller: Optional[RetryController] = None):
        """
        Initialize unified UniProt client.
        
        Args:
            api_config: API configuration (uses default if None)
            retry_controller: Retry controller (uses default if None)
        """
        # Use default configurations if not provided
        if api_config is None:
            system_config = get_config()
            api_config = system_config.api
        
        if retry_controller is None:
            retry_controller = get_retry_controller()
        
        self.api_config = api_config
        self.retry_controller = retry_controller
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        
        # Initialize REST client only
        self.rest_client = UniProtRESTClient(self.api_config, self.retry_controller)
        
        # Set up caching
        self.cache = get_global_cache()
        self.cached_client = CachedAPIClient(self.cache)
    
    async def __aenter__(self):
        """Async context manager entry."""
        await self.rest_client.__aenter__()
        return self
    
    async def __aexit__(self, exc_type, exc_val, exc_tb):
        """Async context manager exit."""
        await self.rest_client.__aexit__(exc_type, exc_val, exc_tb)
    
    async def get_protein_data(self, uniprot_id: str) -> Optional[Dict[str, Any]]:
        """
        Get comprehensive protein data.
        
        Args:
            uniprot_id: UniProt accession ID
            
        Returns:
            Protein data dictionary or None if not found
        """
        try:
            self.logger.info(f"Fetching protein data for {uniprot_id}")
            
            # Use cached client for better performance
            cache_key = f"protein_data:{uniprot_id}"
            
            async def fetch_data():
                return await self.rest_client.get_protein_data(uniprot_id)
            
            result = await self.cached_client.get_or_fetch(cache_key, fetch_data, ttl_hours=24)
            
            if result:
                self.logger.info(f"Successfully retrieved protein data for {uniprot_id}")
            else:
                self.logger.warning(f"No protein data found for {uniprot_id}")
            
            return result
            
        except Exception as e:
            self.logger.error(f"Failed to get protein data for {uniprot_id}: {str(e)}")
            raise
    
    async def get_protein_isoforms(self, uniprot_id: str) -> List[Dict[str, Any]]:
        """
        Get protein isoform data.
        
        Args:
            uniprot_id: UniProt accession ID
            
        Returns:
            List of isoform data dictionaries
        """
        try:
            self.logger.info(f"Fetching isoform data for {uniprot_id}")
            
            # Use cached client for better performance
            cache_key = f"protein_isoforms:{uniprot_id}"
            
            async def fetch_isoforms():
                return await self.rest_client.get_protein_isoforms(uniprot_id)
            
            result = await self.cached_client.get_or_fetch(cache_key, fetch_isoforms, ttl_hours=24)
            
            self.logger.info(f"Retrieved {len(result)} isoforms for {uniprot_id}")
            return result
            
        except Exception as e:
            self.logger.error(f"Failed to get isoform data for {uniprot_id}: {str(e)}")
            raise
    
    async def search_proteins(self, query: str, limit: int = 100) -> List[Dict[str, Any]]:
        """
        Search for proteins.
        
        Args:
            query: Search query
            limit: Maximum number of results
            
        Returns:
            List of protein data dictionaries
        """
        return await self.rest_client.search_proteins(query, limit)
    
    def transform_to_protein_model(self, raw_data: Dict[str, Any], parent_protein_id: str) -> ProteinModel:
        """
        Transform raw UniProt data to ProteinModel.
        
        Args:
            raw_data: Raw protein data from REST API
            parent_protein_id: UniProt ID of the parent protein
            
        Returns:
            ProteinModel instance
        """
        try:
            # Extract basic information
            uniprot_id = raw_data.get("primaryAccession", parent_protein_id)
            accession = raw_data.get("primaryAccession", parent_protein_id)
            name = raw_data.get("uniProtkbId", "")
            
            # Extract protein name
            protein_desc = raw_data.get("proteinDescription", {})
            recommended_name = protein_desc.get("recommendedName", {})
            protein_name = ""
            if recommended_name:
                full_name = recommended_name.get("fullName", {})
                protein_name = full_name.get("value", "")
            
            # Extract organism information
            organism_info = raw_data.get("organism", {})
            organism_name = organism_info.get("scientificName", "")
            organism_id = organism_info.get("taxonId", 0)
            
            # Extract sequence information
            sequence_info = raw_data.get("sequence", {})
            sequence = sequence_info.get("value", "")
            sequence_length = sequence_info.get("length", 0)
            
            # Extract quality indicators
            reviewed = raw_data.get("entryType", "").lower() == "uniprotkb reviewed (swiss-prot)"
            protein_existence = raw_data.get("proteinExistence", "")
            annotation_score = raw_data.get("annotationScore", 0)
            
            # Create ProteinModel
            protein_model = ProteinModel(
                uniprot_id=uniprot_id,
                accession=accession,
                name=name,
                protein_name=protein_name,
                organism_name=organism_name,
                organism_id=organism_id,
                sequence=sequence,
                sequence_length=sequence_length,
                reviewed=reviewed,
                protein_existence=protein_existence,
                annotation_score=annotation_score,
                
                # Set defaults for complex fields (would need more processing)
                alternative_products="{}",
                isoforms="[]",
                isoform_count=1,
                features="[]",
                active_sites="[]",
                binding_sites="[]",
                domains="[]",
                tim_barrel_features="{}",
                secondary_structure="{}",
                interpro_references="[]",
                pfam_references="[]",
                smart_references="[]",
                cdd_references="[]",
                ensembl_references="[]",
                refseq_references="[]",
                embl_references="[]",
                pdb_references="[]",
                comments="[]",
                keywords="[]",
                go_references="[]",
                external_references="{}",
                
                # Metadata
                data_source="uniprot_rest",
                collection_method="rest_api_only",
                last_updated=datetime.now().isoformat(),
                created_at=datetime.now().isoformat()
            )
            
            return protein_model
            
        except Exception as e:
            self.logger.error(f"Failed to transform protein data for {parent_protein_id}: {str(e)}")
            raise DataError(f"Failed to transform protein data: {str(e)}")
    
    async def get_access_method_status(self) -> Dict[str, Any]:
        """
        Get access method status information.
        
        Returns:
            Dictionary with status information
        """
        rest_status = await self.rest_client.get_status()
        
        return {
            "current_method": "rest",
            "rest_available": rest_status["api_available"],
            "rest_client_stats": rest_status,
            "simplified_client": True,
            "rest_only": True
        }
    
    async def get_client_statistics(self) -> Dict[str, Any]:
        """
        Get comprehensive client statistics.
        
        Returns:
            Dictionary with statistics
        """
        rest_status = await self.rest_client.get_status()
        access_method_status = await self.get_access_method_status()
        
        return {
            "client_type": "simplified_rest_only",
            "access_method": access_method_status,
            "rest_client_stats": rest_status,
            "cache_stats": await self.cached_client.get_statistics() if hasattr(self.cached_client, 'get_statistics') else {}
        }