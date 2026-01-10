"""
InterPro API client for PFAM family and protein data collection.

This module provides HTTP client functionality for InterPro REST API endpoints
with comprehensive rate limiting, retry logic, and error handling.
"""

import asyncio
import logging
import time
from typing import Dict, List, Optional, Any, Union
from urllib.parse import urljoin, urlencode
import httpx

from ..config import APIConfig, get_config
from ..retry import RetryController, get_retry_controller
from ..rate_limiter import RateLimitConfig, get_rate_limit_manager, create_api_rate_limiter
from ..cache import get_global_cache, CachedAPIClient
from ..models.entities import TIMBarrelEntryModel, InterProProteinModel
from ..errors import APIError, NetworkError, DataError


class InterProAPIClient:
    """
    HTTP client for InterPro REST API endpoints.
    
    Provides methods for PFAM family queries and protein discovery
    with comprehensive rate limiting and retry logic.
    """
    
    def __init__(self, config: Optional[APIConfig] = None, retry_controller: Optional[RetryController] = None):
        """
        Initialize InterPro API client.
        
        Args:
            config: API configuration, uses global config if not provided
            retry_controller: Retry controller, uses global if not provided
        """
        system_config = get_config()
        self.config = config or system_config.api
        self.retry_controller = retry_controller or get_retry_controller()
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        
        # Set up comprehensive rate limiting
        rate_limit_config = RateLimitConfig(
            requests_per_second=system_config.rate_limiting.interpro_requests_per_second,
            burst_limit=system_config.rate_limiting.interpro_burst_limit,
            burst_window_seconds=system_config.rate_limiting.interpro_burst_window_seconds,
            violation_initial_delay=system_config.rate_limiting.violation_initial_delay,
            violation_backoff_multiplier=system_config.rate_limiting.violation_backoff_multiplier,
            violation_max_delay=system_config.rate_limiting.violation_max_delay,
            soft_limit_threshold=system_config.rate_limiting.soft_limit_threshold,
            enable_monitoring=system_config.rate_limiting.enable_monitoring,
            enable_reporting=system_config.rate_limiting.enable_reporting
        )
        
        # Create rate limiter for InterPro API
        rate_limit_manager = get_rate_limit_manager()
        self.rate_limiter = rate_limit_manager.create_limiter("InterPro", rate_limit_config)
        
        # Set up response caching
        self.cache = get_global_cache()
        self.cached_client = CachedAPIClient(self.cache)
        
        # HTTP client configuration
        self.client = httpx.AsyncClient(
            timeout=httpx.Timeout(
                connect=self.config.connection_timeout,
                read=self.config.request_timeout,
                write=self.config.request_timeout,
                pool=self.config.request_timeout
            ),
            headers={
                "User-Agent": "ProteinDataCollector/1.0",
                "Accept": "application/json"
            }
        )
    
    async def __aenter__(self):
        """Async context manager entry."""
        return self
    
    async def __aexit__(self, exc_type, exc_val, exc_tb):
        """Async context manager exit."""
        await self.client.aclose()
    
    def _build_url(self, endpoint: str, params: Optional[Dict[str, Any]] = None) -> str:
        """
        Build complete URL for InterPro API endpoint.
        
        Args:
            endpoint: API endpoint path
            params: Query parameters
            
        Returns:
            Complete URL with parameters
        """
        url = urljoin(self.config.interpro_base_url, endpoint.lstrip('/'))
        if params:
            # Filter out None values and convert to strings
            clean_params = {k: str(v) for k, v in params.items() if v is not None}
            if clean_params:
                url += '?' + urlencode(clean_params)
        return url
    
    async def _make_request(self, endpoint: str, params: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        Make HTTP request to InterPro API with rate limiting, caching, and error handling.
        
        Args:
            endpoint: API endpoint path
            params: Query parameters
            
        Returns:
            JSON response data
            
        Raises:
            APIError: For API-specific errors
            NetworkError: For network-related errors
            DataError: For invalid response data
        """
        # Use cached request if caching is enabled
        async def request_func():
            return await self._make_direct_request(endpoint, params)
        
        return await self.cached_client.cached_request(
            api_name="InterPro",
            endpoint=endpoint,
            params=params,
            request_func=request_func
        )
    
    async def _make_direct_request(self, endpoint: str, params: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        Make direct HTTP request to InterPro API without caching.
        
        Args:
            endpoint: API endpoint path
            params: Query parameters
            
        Returns:
            JSON response data
            
        Raises:
            APIError: For API-specific errors
            NetworkError: For network-related errors
            DataError: For invalid response data
        """
        url = self._build_url(endpoint, params)
        
        # Apply comprehensive rate limiting
        delay = await self.rate_limiter.acquire()
        
        self.logger.debug(
            "Making InterPro API request",
            extra={
                "url": url,
                "endpoint": endpoint,
                "params": params,
                "rate_limit_delay": delay
            }
        )
        
        try:
            response = await self.client.get(url)
            
            # Log response details
            self.logger.debug(
                "InterPro API response received",
                extra={
                    "status_code": response.status_code,
                    "response_size": len(response.content),
                    "url": url
                }
            )
            
            # Handle HTTP errors
            if response.status_code == 204:
                # No content - return empty response
                return {"results": [], "count": 0, "next": None}
            elif response.status_code == 429:
                raise APIError(f"Rate limit exceeded for InterPro API: {url}", response.status_code)
            elif response.status_code >= 500:
                raise APIError(f"InterPro API server error: {response.status_code}", response.status_code)
            elif response.status_code >= 400:
                raise APIError(f"InterPro API client error: {response.status_code}", response.status_code)
            
            # Parse JSON response
            try:
                data = response.json()
            except Exception as e:
                # Handle empty response body
                if len(response.content) == 0:
                    return {"results": [], "count": 0, "next": None}
                raise DataError(f"Invalid JSON response from InterPro API: {str(e)}")
            
            return data
            
        except httpx.TimeoutException as e:
            raise NetworkError(f"Timeout connecting to InterPro API: {str(e)}")
        except httpx.ConnectError as e:
            raise NetworkError(f"Connection error to InterPro API: {str(e)}")
        except httpx.HTTPError as e:
            raise NetworkError(f"HTTP error from InterPro API: {str(e)}")
    
    async def get_pfam_families_with_tim_barrel(self, page_size: int = 200) -> List[Dict[str, Any]]:
        """
        Retrieve PFAM families with TIM barrel annotations from InterPro.
        
        Args:
            page_size: Number of results per page
            
        Returns:
            List of PFAM family data dictionaries
            
        Raises:
            APIError: For API-specific errors
            NetworkError: For network-related errors
        """
        async def operation():
            families = []
            next_url = None
            page = 1
            
            # Search for PFAM entries with TIM barrel annotations
            # Using InterPro's search functionality for structural annotations
            params = {
                'page_size': page_size,
                'search': 'TIM barrel'
            }
            
            while True:
                self.logger.info(
                    "Fetching PFAM families page %d",
                    page,
                    extra={"page": page, "page_size": page_size}
                )
                
                endpoint = 'entry/pfam/'
                if next_url:
                    # Use the next URL provided by the API
                    response_data = await self._make_request(next_url.replace(self.config.interpro_base_url, ''))
                else:
                    response_data = await self._make_request(endpoint, params)
                
                # Extract results from response
                if 'results' in response_data:
                    page_families = response_data['results']
                    families.extend(page_families)
                    
                    self.logger.info(
                        "Retrieved %d PFAM families from page %d",
                        len(page_families),
                        page,
                        extra={
                            "page": page,
                            "families_in_page": len(page_families),
                            "total_families": len(families)
                        }
                    )
                
                # Check for next page
                next_url = response_data.get('next')
                if not next_url:
                    break
                
                page += 1
            
            self.logger.info(
                "Completed PFAM family collection",
                extra={
                    "total_families": len(families),
                    "total_pages": page
                }
            )
            
            return families
        
        return await self.retry_controller.execute_with_retry_async(
            operation,
            database="InterPro",
            operation_name="get_pfam_families_with_tim_barrel"
        )
    
    async def get_pfam_families_with_tim_barrel_search(self, search_term: str, page_size: int = 200) -> List[Dict[str, Any]]:
        """
        Retrieve PFAM families using a specific search term from InterPro.
        
        Args:
            search_term: Search term to use
            page_size: Number of results per page
            
        Returns:
            List of PFAM family data dictionaries
            
        Raises:
            APIError: For API-specific errors
            NetworkError: For network-related errors
        """
        async def operation():
            families = []
            next_url = None
            page = 1
            
            # Search for PFAM entries with the given search term
            params = {
                'page_size': page_size,
                'search': search_term
            }
            
            while True:
                self.logger.debug(
                    f"Fetching PFAM families page {page} for search '{search_term}'",
                    extra={"page": page, "page_size": page_size, "search_term": search_term}
                )
                
                endpoint = 'entry/pfam/'
                if next_url:
                    # Use the next URL provided by the API
                    response_data = await self._make_request(next_url.replace(self.config.interpro_base_url, ''))
                else:
                    response_data = await self._make_request(endpoint, params)
                
                # Extract results from response
                if 'results' in response_data:
                    page_families = response_data['results']
                    families.extend(page_families)
                    
                    self.logger.debug(
                        f"Retrieved {len(page_families)} PFAM families from page {page} for '{search_term}'",
                        extra={
                            "page": page,
                            "families_in_page": len(page_families),
                            "total_families": len(families),
                            "search_term": search_term
                        }
                    )
                
                # Check for next page
                next_url = response_data.get('next')
                if not next_url:
                    break
                
                page += 1
            
            self.logger.debug(
                f"Completed PFAM family collection for '{search_term}'",
                extra={
                    "total_families": len(families),
                    "total_pages": page,
                    "search_term": search_term
                }
            )
            
            return families
        
        return await self.retry_controller.execute_with_retry_async(
            operation,
            database="InterPro",
            operation_name=f"get_pfam_families_search_{search_term.replace(' ', '_')}"
        )
    
    async def get_interpro_entries_with_tim_barrel_search(self, search_term: str, page_size: int = 200) -> List[Dict[str, Any]]:
        """
        Retrieve InterPro entries (IPR) using a specific search term from InterPro.
        
        Args:
            search_term: Search term to use
            page_size: Number of results per page
            
        Returns:
            List of InterPro entry data dictionaries
            
        Raises:
            APIError: For API-specific errors
            NetworkError: For network-related errors
        """
        async def operation():
            entries = []
            next_url = None
            page = 1
            
            # Search for InterPro entries with the given search term
            params = {
                'page_size': page_size,
                'search': search_term
            }
            
            while True:
                self.logger.debug(
                    f"Fetching InterPro entries page {page} for search '{search_term}'",
                    extra={"page": page, "page_size": page_size, "search_term": search_term}
                )
                
                endpoint = 'entry/interpro/'
                if next_url:
                    # Use the next URL provided by the API
                    response_data = await self._make_request(next_url.replace(self.config.interpro_base_url, ''))
                else:
                    response_data = await self._make_request(endpoint, params)
                
                # Extract results from response
                if 'results' in response_data:
                    page_entries = response_data['results']
                    entries.extend(page_entries)
                    
                    self.logger.debug(
                        f"Retrieved {len(page_entries)} InterPro entries from page {page} for '{search_term}'",
                        extra={
                            "page": page,
                            "entries_in_page": len(page_entries),
                            "total_entries": len(entries),
                            "search_term": search_term
                        }
                    )
                
                # Check for next page
                next_url = response_data.get('next')
                if not next_url:
                    break
                
                page += 1
            
            self.logger.debug(
                f"Completed InterPro entry collection for '{search_term}'",
                extra={
                    "total_entries": len(entries),
                    "total_pages": page,
                    "search_term": search_term
                }
            )
            
            return entries
        
        return await self.retry_controller.execute_with_retry_async(
            operation,
            database="InterPro",
            operation_name=f"get_interpro_entries_search_{search_term.replace(' ', '_')}"
        )

    async def get_pfam_families_from_interpro_entry(self, interpro_accession: str) -> List[Dict[str, Any]]:
        """
        Get PFAM families associated with a specific InterPro entry.
        
        Args:
            interpro_accession: InterPro entry accession identifier (e.g., IPR013785)
            
        Returns:
            List of PFAM family data dictionaries
            
        Raises:
            APIError: For API-specific errors
            NetworkError: For network-related errors
        """
        async def operation():
            endpoint = f'entry/interpro/{interpro_accession}/'
            
            self.logger.debug(
                "Fetching InterPro entry details",
                extra={"interpro_accession": interpro_accession}
            )
            
            entry_data = await self._make_request(endpoint)
            
            # Extract PFAM families from member databases
            pfam_families = []
            member_databases = entry_data.get('member_databases', {})
            
            if 'pfam' in member_databases:
                pfam_entries = member_databases['pfam']
                for pfam_entry in pfam_entries:
                    # Get detailed information for each PFAM family
                    pfam_accession = pfam_entry.get('accession')
                    if pfam_accession:
                        try:
                            pfam_details = await self.get_pfam_family_details(pfam_accession)
                            pfam_families.append(pfam_details)
                        except Exception as e:
                            self.logger.warning(
                                f"Failed to get details for PFAM {pfam_accession} from InterPro {interpro_accession}: {str(e)}"
                            )
            
            self.logger.debug(
                f"Found {len(pfam_families)} PFAM families in InterPro entry {interpro_accession}",
                extra={
                    "interpro_accession": interpro_accession,
                    "pfam_count": len(pfam_families)
                }
            )
            
            return pfam_families
        
        return await self.retry_controller.execute_with_retry_async(
            operation,
            database="InterPro",
            operation_name=f"get_pfam_families_from_interpro_{interpro_accession}"
        )

    async def get_pfam_family_details(self, pfam_accession: str) -> Dict[str, Any]:
        """
        Get detailed information for a specific PFAM family.
        
        Args:
            pfam_accession: PFAM family accession identifier
            
        Returns:
            PFAM family details dictionary
            
        Raises:
            APIError: For API-specific errors
            NetworkError: For network-related errors
        """
        async def operation():
            endpoint = f'entry/pfam/{pfam_accession}/'
            
            self.logger.debug(
                "Fetching PFAM family details",
                extra={"pfam_accession": pfam_accession}
            )
            
            return await self._make_request(endpoint)
        
        return await self.retry_controller.execute_with_retry_async(
            operation,
            database="InterPro",
            operation_name=f"get_pfam_family_details_{pfam_accession}"
        )
    
    async def get_proteins_in_pfam_family(
        self, 
        pfam_accession: str, 
        organism: str = "Homo sapiens",
        page_size: int = 200
    ) -> List[Dict[str, Any]]:
        """
        Get human proteins belonging to a specific PFAM family.
        
        Args:
            pfam_accession: PFAM family accession identifier
            organism: Target organism (default: Homo sapiens)
            page_size: Number of results per page
            
        Returns:
            List of protein data dictionaries
            
        Raises:
            APIError: For API-specific errors
            NetworkError: For network-related errors
        """
        async def operation():
            proteins = []
            next_url = None
            page = 1
            
            # Map organism to taxonomy ID
            organism_tax_id = "9606" if organism == "Homo sapiens" else organism
            
            # Query parameters
            params = {
                'page_size': page_size
            }
            
            while True:
                self.logger.info(
                    "Fetching proteins for PFAM %s, page %d",
                    pfam_accession,
                    page,
                    extra={
                        "pfam_accession": pfam_accession,
                        "organism": organism,
                        "taxonomy_id": organism_tax_id,
                        "page": page,
                        "page_size": page_size
                    }
                )
                
                # Use correct endpoint: taxonomy first, then entry
                endpoint = f'protein/UniProt/taxonomy/uniprot/{organism_tax_id}/entry/pfam/{pfam_accession}/'
                if next_url:
                    # Use the next URL provided by the API
                    response_data = await self._make_request(next_url.replace(self.config.interpro_base_url, ''))
                else:
                    response_data = await self._make_request(endpoint, params)
                
                # Extract results from response
                if 'results' in response_data:
                    page_proteins = response_data['results']
                    proteins.extend(page_proteins)
                    
                    self.logger.info(
                        "Retrieved %d proteins from page %d for PFAM %s",
                        len(page_proteins),
                        page,
                        pfam_accession,
                        extra={
                            "pfam_accession": pfam_accession,
                            "page": page,
                            "proteins_in_page": len(page_proteins),
                            "total_proteins": len(proteins)
                        }
                    )
                
                # Check for next page
                next_url = response_data.get('next')
                if not next_url:
                    break
                
                page += 1
            
            self.logger.info(
                "Completed protein collection for PFAM %s",
                pfam_accession,
                extra={
                    "pfam_accession": pfam_accession,
                    "total_proteins": len(proteins),
                    "total_pages": page,
                    "organism": organism
                }
            )
            
            return proteins
        
        return await self.retry_controller.execute_with_retry_async(
            operation,
            database="InterPro",
            operation_name=f"get_proteins_in_pfam_family_{pfam_accession}"
        )
    
    def parse_interpro_entry_data(self, entry_data: Dict[str, Any]) -> TIMBarrelEntryModel:
        """
        Parse InterPro API response data into TIMBarrelEntryModel.
        
        Args:
            entry_data: Raw InterPro entry data from InterPro API
            
        Returns:
            Validated TIMBarrelEntryModel instance
            
        Raises:
            DataError: For invalid or incomplete data
        """
        try:
            # Extract metadata from InterPro response
            metadata = entry_data.get('metadata', {})
            if not metadata:
                raise DataError("Missing metadata in InterPro entry data")
            
            # Extract required fields
            accession = metadata.get('accession')
            if not accession:
                raise DataError("Missing InterPro accession in entry data")
            
            name = metadata.get('name', '')
            if not name:
                raise DataError("Missing InterPro name in entry data")
            
            # Extract optional fields
            description = metadata.get('description', name)
            interpro_type = metadata.get('type', 'Unknown')
            
            # Create TIM barrel annotation from the name and description
            tim_barrel_annotation = f"TIM barrel InterPro entry: {name}"
            if description and description != name:
                tim_barrel_annotation += f" - {description}"
            
            # Extract member databases
            member_databases = entry_data.get('member_databases', {})
            
            return TIMBarrelEntryModel(
                accession=accession,
                entry_type='interpro',
                name=name,
                description=description,
                interpro_type=interpro_type,
                tim_barrel_annotation=tim_barrel_annotation,
                member_databases=member_databases
            )
            
        except Exception as e:
            raise DataError(f"Failed to parse InterPro entry data: {str(e)}")

    def parse_pfam_family_data(self, family_data: Dict[str, Any]) -> TIMBarrelEntryModel:
        """
        Parse InterPro API response data into TIMBarrelEntryModel.
        
        Args:
            family_data: Raw family data from InterPro API
            
        Returns:
            Validated TIMBarrelEntryModel instance
            
        Raises:
            DataError: For invalid or incomplete data
        """
        try:
            # Extract metadata from InterPro response
            metadata = family_data.get('metadata', {})
            if not metadata:
                raise DataError("Missing metadata in family data")
            
            # Extract required fields
            accession = metadata.get('accession')
            if not accession:
                raise DataError("Missing PFAM accession in family data")
            
            name = metadata.get('name', '')
            if not name:
                raise DataError("Missing PFAM name in family data")
            
            # Use name as description since InterPro doesn't provide separate description
            description = name
            
            # Create TIM barrel annotation from the name (since we searched for TIM barrel)
            tim_barrel_annotation = f"TIM barrel family: {name}"
            
            # Get InterPro ID from integrated field
            interpro_id = metadata.get('integrated')
            
            return TIMBarrelEntryModel(
                accession=accession,
                entry_type='pfam',
                name=name,
                description=description,
                tim_barrel_annotation=tim_barrel_annotation,
                interpro_id=interpro_id
            )
            
        except Exception as e:
            raise DataError(f"Failed to parse PFAM family data: {str(e)}")
    
    def parse_protein_data(self, protein_data: Dict[str, Any], tim_barrel_accession: str) -> InterProProteinModel:
        """
        Parse InterPro API response data into InterProProteinModel.
        
        Args:
            protein_data: Raw protein data from InterPro API
            tim_barrel_accession: Associated TIM barrel entry accession (PFAM or InterPro)
            
        Returns:
            Validated InterProProteinModel instance
            
        Raises:
            DataError: For invalid or incomplete data
        """
        try:
            # Extract UniProt ID from metadata
            uniprot_id = protein_data.get('metadata', {}).get('accession')
            if not uniprot_id:
                raise DataError("Missing UniProt ID in protein data")
            
            # Extract protein name and organism
            name = protein_data.get('metadata', {}).get('name', '')
            
            # Extract organism information
            organism = "Homo sapiens"  # Default, as we filter by this
            source_organism = protein_data.get('metadata', {}).get('source_organism', {})
            if isinstance(source_organism, dict):
                organism_name = source_organism.get('fullName', source_organism.get('scientificName', ''))
                if organism_name:
                    organism = organism_name
            
            # Extract additional metadata
            gene_info = protein_data.get('metadata', {}).get('gene', '')
            gene_name = ''
            if isinstance(gene_info, dict):
                gene_name = gene_info.get('name', '')
            elif isinstance(gene_info, str):
                gene_name = gene_info
            
            basic_metadata = {
                'source_database': protein_data.get('metadata', {}).get('source_database', ''),
                'length': protein_data.get('metadata', {}).get('length'),
                'gene_name': gene_name,
                'protein_existence': protein_data.get('metadata', {}).get('protein_existence')
            }
            
            # Remove None values from metadata
            basic_metadata = {k: v for k, v in basic_metadata.items() if v is not None}
            
            return InterProProteinModel(
                uniprot_id=uniprot_id,
                tim_barrel_accession=tim_barrel_accession,
                name=name,
                organism=organism,
                basic_metadata=basic_metadata
            )
            
        except Exception as e:
            raise DataError(f"Failed to parse protein data: {str(e)}")
    
    async def get_proteins_in_interpro_entry(
        self, 
        interpro_accession: str, 
        organism: str = "Homo sapiens",
        page_size: int = 200
    ) -> List[Dict[str, Any]]:
        """
        Get human proteins belonging to a specific InterPro entry.
        
        Args:
            interpro_accession: InterPro entry accession identifier
            organism: Target organism (default: Homo sapiens)
            page_size: Number of results per page
            
        Returns:
            List of protein data dictionaries
            
        Raises:
            APIError: For API-specific errors
            NetworkError: For network-related errors
        """
        async def operation():
            proteins = []
            next_url = None
            page = 1
            
            # Map organism to taxonomy ID
            organism_tax_id = "9606" if organism == "Homo sapiens" else organism
            
            # Query parameters
            params = {
                'page_size': page_size
            }
            
            while True:
                self.logger.info(
                    "Fetching proteins for InterPro %s, page %d",
                    interpro_accession,
                    page,
                    extra={
                        "interpro_accession": interpro_accession,
                        "organism": organism,
                        "taxonomy_id": organism_tax_id,
                        "page": page,
                        "page_size": page_size
                    }
                )
                
                # Use correct endpoint: taxonomy first, then entry
                endpoint = f'protein/UniProt/taxonomy/uniprot/{organism_tax_id}/entry/interpro/{interpro_accession}/'
                if next_url:
                    # Use the next URL provided by the API
                    response_data = await self._make_request(next_url.replace(self.config.interpro_base_url, ''))
                else:
                    response_data = await self._make_request(endpoint, params)
                
                # Extract results from response
                if 'results' in response_data:
                    page_proteins = response_data['results']
                    proteins.extend(page_proteins)
                    
                    self.logger.info(
                        "Retrieved %d proteins from page %d for InterPro %s",
                        len(page_proteins),
                        page,
                        interpro_accession,
                        extra={
                            "interpro_accession": interpro_accession,
                            "page": page,
                            "proteins_in_page": len(page_proteins),
                            "total_proteins": len(proteins)
                        }
                    )
                
                # Check for next page
                next_url = response_data.get('next')
                if not next_url:
                    break
                
                page += 1
            
            self.logger.info(
                "Completed protein collection for InterPro %s",
                interpro_accession,
                extra={
                    "interpro_accession": interpro_accession,
                    "total_proteins": len(proteins),
                    "total_pages": page,
                    "organism": organism
                }
            )
            
            return proteins
        
        return await self.retry_controller.execute_with_retry_async(
            operation,
            database="InterPro",
            operation_name=f"get_proteins_in_interpro_entry_{interpro_accession}"
        )

    async def get_performance_metrics(self) -> Dict[str, Any]:
        """
        Get comprehensive performance metrics for InterPro API client.
        
        Returns:
            Dictionary with performance metrics
        """
        cache_metrics = self.cache.get_metrics()
        rate_limit_stats = self.rate_limiter.get_stats()
        
        return {
            "api_name": "InterPro",
            "cache_metrics": cache_metrics,
            "rate_limiting": rate_limit_stats,
            "client_config": {
                "base_url": self.config.interpro_base_url,
                "connection_timeout": self.config.connection_timeout,
                "request_timeout": self.config.request_timeout
            }
        }
    
    async def close(self) -> None:
        """Close the HTTP client and cleanup resources."""
        await self.client.aclose()
        if hasattr(self.cache, 'close'):
            await self.cache.close()


# Convenience functions for common operations
async def get_tim_barrel_entries(
    client: Optional[InterProAPIClient] = None,
    page_size: int = 200
) -> List[TIMBarrelEntryModel]:
    """
    Get all TIM barrel entries (both PFAM families and InterPro entries).
    
    Args:
        client: InterPro API client, creates new one if not provided
        page_size: Number of results per page
        
    Returns:
        List of validated TIMBarrelEntryModel instances
    """
    if client is None:
        async with InterProAPIClient() as client:
            return await get_tim_barrel_entries(client, page_size)
    
    # Get both PFAM families and InterPro entries
    pfam_data = await client.get_pfam_families_with_tim_barrel(page_size)
    interpro_data = await client.get_interpro_entries_with_tim_barrel_search("TIM barrel", page_size)
    
    entries = []
    
    # Parse PFAM families
    for family_data in pfam_data:
        try:
            entry = client.parse_pfam_family_data(family_data)
            entries.append(entry)
        except DataError as e:
            client.logger.warning(
                "Skipping invalid PFAM family data: %s",
                str(e),
                extra={"family_data": family_data}
            )
    
    # Parse InterPro entries
    for entry_data in interpro_data:
        try:
            entry = client.parse_interpro_entry_data(entry_data)
            entries.append(entry)
        except DataError as e:
            client.logger.warning(
                "Skipping invalid InterPro entry data: %s",
                str(e),
                extra={"entry_data": entry_data}
            )
    
    return entries


async def get_human_proteins_for_tim_barrel_entries(
    tim_barrel_entries: List[TIMBarrelEntryModel],
    client: Optional[InterProAPIClient] = None,
    page_size: int = 200
) -> List[InterProProteinModel]:
    """
    Get all human proteins for a list of TIM barrel entries.
    
    Args:
        tim_barrel_entries: List of TIM barrel entries to query
        client: InterPro API client, creates new one if not provided
        page_size: Number of results per page
        
    Returns:
        List of validated InterProProteinModel instances
    """
    if client is None:
        async with InterProAPIClient() as client:
            return await get_human_proteins_for_tim_barrel_entries(tim_barrel_entries, client, page_size)
    
    all_proteins = []
    
    for entry in tim_barrel_entries:
        # Only get proteins for PFAM entries (InterPro entries don't have direct protein associations)
        if entry.is_pfam:
            try:
                proteins_data = await client.get_proteins_in_pfam_family(
                    entry.accession, 
                    organism="Homo sapiens",
                    page_size=page_size
                )
                
                for protein_data in proteins_data:
                    try:
                        protein = client.parse_protein_data(protein_data, entry.accession)
                        all_proteins.append(protein)
                    except DataError as e:
                        client.logger.warning(
                            "Skipping invalid protein data for entry %s: %s",
                            entry.accession,
                            str(e),
                            extra={"protein_data": protein_data, "entry_accession": entry.accession}
                        )
                        
            except Exception as e:
                client.logger.error(
                    "Failed to get proteins for TIM barrel entry %s: %s",
                    entry.accession,
                    str(e),
                    extra={"entry_accession": entry.accession}
                )
    
    return all_proteins