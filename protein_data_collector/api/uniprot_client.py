"""
Unified UniProt API client with MCP server integration and REST API fallback.

This module provides a unified interface for accessing UniProt data through
either MCP servers (preferred) or direct REST API calls (fallback), abstracting
the differences between the two approaches with comprehensive rate limiting.
"""

import asyncio
import json
import logging
import subprocess
import tempfile
from typing import Dict, List, Optional, Any, Union
from urllib.parse import urljoin, urlencode
import httpx

from ..config import APIConfig, MCPConfig, get_config
from ..retry import RetryController, get_retry_controller
from ..rate_limiter import RateLimitConfig, get_rate_limit_manager
from ..cache import get_global_cache, CachedAPIClient
from ..models.entities import ProteinModel
from ..errors import APIError, NetworkError, DataError, ConfigurationError, ErrorContext, create_error_context


class MCPServerClient:
    """Client for interacting with UniProt MCP servers with rate limiting."""
    
    def __init__(self, mcp_config: MCPConfig, retry_controller: RetryController):
        """
        Initialize MCP server client.
        
        Args:
            mcp_config: MCP server configuration
            retry_controller: Retry controller for operations
        """
        self.config = mcp_config
        self.retry_controller = retry_controller
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        self._server_available = None  # Cache server availability status
        
        # Set up rate limiting for MCP server
        system_config = get_config()
        rate_limit_config = RateLimitConfig(
            requests_per_second=system_config.rate_limiting.mcp_requests_per_second,
            burst_limit=system_config.rate_limiting.mcp_burst_limit,
            burst_window_seconds=system_config.rate_limiting.mcp_burst_window_seconds,
            violation_initial_delay=system_config.rate_limiting.violation_initial_delay,
            violation_backoff_multiplier=system_config.rate_limiting.violation_backoff_multiplier,
            violation_max_delay=system_config.rate_limiting.violation_max_delay,
            soft_limit_threshold=system_config.rate_limiting.soft_limit_threshold,
            enable_monitoring=system_config.rate_limiting.enable_monitoring,
            enable_reporting=system_config.rate_limiting.enable_reporting
        )
        
        # Create rate limiter for MCP server
        rate_limit_manager = get_rate_limit_manager()
        self.rate_limiter = rate_limit_manager.create_limiter("UniProt_MCP", rate_limit_config)
    
    async def is_server_available(self) -> bool:
        """
        Check if MCP server is available and functional.
        
        Returns:
            True if server is available, False otherwise
        """
        if not self.config.enabled:
            return False
        
        # Return cached result if available
        if self._server_available is not None:
            return self._server_available
        
        try:
            # Test server availability with a simple query
            result = await self._call_mcp_server("list_tools", {})
            self._server_available = result is not None
            
            if self._server_available:
                self.logger.info("MCP server is available and functional")
            else:
                self.logger.warning("MCP server is not responding properly")
                
            return self._server_available
            
        except Exception as e:
            self.logger.warning(
                "MCP server availability check failed: %s",
                str(e),
                extra={"server_path": self.config.server_path}
            )
            self._server_available = False
            return False
    
    async def _call_mcp_server(self, tool_name: str, parameters: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """
        Call MCP server with specified tool and parameters.
        
        Args:
            tool_name: Name of the MCP tool to call
            parameters: Parameters for the tool
            
        Returns:
            MCP server response data or None if failed
            
        Raises:
            ConfigurationError: If MCP server is not properly configured
            APIError: If MCP server call fails
        """
        if not self.config.enabled:
            raise ConfigurationError("MCP server is disabled in configuration")
        
        # Apply rate limiting
        delay = await self.rate_limiter.acquire()
        
        # Prepare MCP request
        mcp_request = {
            "jsonrpc": "2.0",
            "id": 1,
            "method": "tools/call",
            "params": {
                "name": tool_name,
                "arguments": parameters
            }
        }
        
        self.logger.debug(
            "Making MCP server request",
            extra={
                "tool_name": tool_name,
                "parameters": parameters,
                "rate_limit_delay": delay
            }
        )
        
        try:
            # Create temporary file for MCP communication
            with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as temp_file:
                json.dump(mcp_request, temp_file)
                temp_file_path = temp_file.name
            
            # Call MCP server using subprocess
            process = await asyncio.create_subprocess_exec(
                self.config.server_path,
                stdin=asyncio.subprocess.PIPE,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE
            )
            
            # Send request and get response
            request_data = json.dumps(mcp_request).encode('utf-8')
            stdout, stderr = await process.communicate(input=request_data)
            
            if process.returncode != 0:
                error_msg = stderr.decode('utf-8') if stderr else "Unknown MCP server error"
                raise APIError(f"MCP server call failed: {error_msg}")
            
            # Parse response
            response_text = stdout.decode('utf-8')
            response_data = json.loads(response_text)
            
            # Check for JSON-RPC errors
            if "error" in response_data:
                error_info = response_data["error"]
                raise APIError(f"MCP server error: {error_info.get('message', 'Unknown error')}")
            
            # Return result
            return response_data.get("result")
            
        except json.JSONDecodeError as e:
            raise DataError(f"Invalid JSON response from MCP server: {str(e)}")
        except Exception as e:
            if isinstance(e, (APIError, DataError)):
                raise
            raise APIError(f"MCP server communication failed: {str(e)}")
    
    async def get_protein_data(self, uniprot_id: str) -> Optional[Dict[str, Any]]:
        """
        Get protein data from UniProt via MCP server.
        
        Args:
            uniprot_id: UniProt protein identifier
            
        Returns:
            Protein data dictionary or None if not found
        """
        if not await self.is_server_available():
            return None
        
        try:
            # Call MCP server to get protein information
            result = await self._call_mcp_server("get_protein", {"accession": uniprot_id})
            
            if result and "protein_data" in result:
                self.logger.debug(
                    "Retrieved protein data via MCP server",
                    extra={"uniprot_id": uniprot_id}
                )
                return result["protein_data"]
            
            return None
            
        except Exception as e:
            self.logger.warning(
                "Failed to get protein data via MCP server: %s",
                str(e),
                extra={"uniprot_id": uniprot_id}
            )
            return None
    
    async def get_protein_isoforms(self, uniprot_id: str) -> List[Dict[str, Any]]:
        """
        Get all isoforms for a protein via MCP server.
        
        Args:
            uniprot_id: UniProt protein identifier
            
        Returns:
            List of isoform data dictionaries
        """
        if not await self.is_server_available():
            return []
        
        try:
            # Call MCP server to get isoform information
            result = await self._call_mcp_server("get_isoforms", {"accession": uniprot_id})
            
            if result and "isoforms" in result:
                isoforms = result["isoforms"]
                self.logger.debug(
                    "Retrieved %d isoforms via MCP server",
                    len(isoforms),
                    extra={"uniprot_id": uniprot_id, "isoform_count": len(isoforms)}
                )
                return isoforms
            
            return []
            
        except Exception as e:
            self.logger.warning(
                "Failed to get protein isoforms via MCP server: %s",
                str(e),
                extra={"uniprot_id": uniprot_id}
            )
            return []


class UniProtRESTClient:
    """Client for direct UniProt REST API access with comprehensive rate limiting."""
    
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
        
        # Set up comprehensive rate limiting
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
        Build complete URL for UniProt REST API endpoint.
        
        Args:
            endpoint: API endpoint path
            params: Query parameters
            
        Returns:
            Complete URL with parameters
        """
        url = urljoin(self.config.uniprot_base_url, endpoint.lstrip('/'))
        if params:
            # Filter out None values and convert to strings
            clean_params = {k: str(v) for k, v in params.items() if v is not None}
            if clean_params:
                url += '?' + urlencode(clean_params)
        return url
    
    async def _make_request(self, endpoint: str, params: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        Make HTTP request to UniProt REST API with caching.
        
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
            api_name="UniProt_REST",
            endpoint=endpoint,
            params=params,
            request_func=request_func
        )
    
    async def _make_direct_request(self, endpoint: str, params: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        Make direct HTTP request to UniProt REST API without caching.
        
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
            "Making UniProt REST API request",
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
                "UniProt REST API response received",
                extra={
                    "status_code": response.status_code,
                    "response_size": len(response.content),
                    "url": url
                }
            )
            
            # Handle HTTP errors
            if response.status_code == 429:
                raise APIError(f"Rate limit exceeded for UniProt API: {url}", response.status_code)
            elif response.status_code >= 500:
                raise APIError(f"UniProt API server error: {response.status_code}", response.status_code)
            elif response.status_code == 404:
                # Return empty result for not found
                return {}
            elif response.status_code >= 400:
                raise APIError(f"UniProt API client error: {response.status_code}", response.status_code)
            
            # Parse JSON response
            try:
                data = response.json()
            except Exception as e:
                raise DataError(f"Invalid JSON response from UniProt API: {str(e)}")
            
            return data
            
        except httpx.TimeoutException as e:
            raise NetworkError(f"Timeout connecting to UniProt API: {str(e)}")
        except httpx.ConnectError as e:
            raise NetworkError(f"Connection error to UniProt API: {str(e)}")
        except httpx.HTTPError as e:
            raise NetworkError(f"HTTP error from UniProt API: {str(e)}")
    
    async def get_protein_data(self, uniprot_id: str) -> Optional[Dict[str, Any]]:
        """
        Get protein data from UniProt REST API.
        
        Args:
            uniprot_id: UniProt protein identifier
            
        Returns:
            Protein data dictionary or None if not found
        """
        try:
            # Request specific fields for protein data
            fields = [
                "accession", "id", "protein_name", "gene_names", "organism_name",
                "sequence", "length", "mass", "cc_alternative_products",
                "ft_region", "ft_domain", "xref_interpro"
            ]
            
            params = {
                "fields": ",".join(fields),
                "format": "json"
            }
            
            endpoint = f"uniprotkb/{uniprot_id}"
            data = await self._make_request(endpoint, params)
            
            if data:
                self.logger.debug(
                    "Retrieved protein data via REST API",
                    extra={"uniprot_id": uniprot_id}
                )
                return data
            
            return None
            
        except APIError as e:
            if e.status_code == 404:
                return None
            raise
    
    async def get_protein_isoforms(self, uniprot_id: str) -> List[Dict[str, Any]]:
        """
        Get all isoforms for a protein via UniProt REST API.
        
        Args:
            uniprot_id: UniProt protein identifier
            
        Returns:
            List of isoform data dictionaries
        """
        try:
            # Get protein data which includes isoform information
            protein_data = await self.get_protein_data(uniprot_id)
            
            if not protein_data:
                return []
            
            isoforms = []
            
            # Extract canonical sequence as first isoform
            canonical_isoform = {
                "isoform_id": f"{uniprot_id}-1",
                "sequence": protein_data.get("sequence", {}).get("value", ""),
                "length": protein_data.get("sequence", {}).get("length", 0),
                "is_canonical": True,
                "description": "Canonical sequence"
            }
            isoforms.append(canonical_isoform)
            
            # Extract alternative isoforms from cc_alternative_products
            alt_products = protein_data.get("comments", [])
            for comment in alt_products:
                if comment.get("commentType") == "ALTERNATIVE_PRODUCTS":
                    for isoform in comment.get("isoforms", []):
                        isoform_data = {
                            "isoform_id": isoform.get("name", {}).get("value", ""),
                            "sequence": isoform.get("sequence", {}).get("value", ""),
                            "length": isoform.get("sequence", {}).get("length", 0),
                            "is_canonical": False,
                            "description": isoform.get("note", {}).get("texts", [{}])[0].get("value", "")
                        }
                        isoforms.append(isoform_data)
            
            self.logger.debug(
                "Retrieved %d isoforms via REST API",
                len(isoforms),
                extra={"uniprot_id": uniprot_id, "isoform_count": len(isoforms)}
            )
            
            return isoforms
            
        except Exception as e:
            self.logger.warning(
                "Failed to get protein isoforms via REST API: %s",
                str(e),
                extra={"uniprot_id": uniprot_id}
            )
            return []
    
    async def close(self) -> None:
        """Close the HTTP client."""
        await self.client.aclose()


class UnifiedUniProtClient:
    """
    Unified client for UniProt data access with MCP server integration and REST API fallback.
    
    This client provides a single interface that abstracts the differences between
    MCP server calls and direct REST API calls, automatically falling back to
    REST API when MCP servers are unavailable.
    """
    
    def __init__(
        self,
        api_config: Optional[APIConfig] = None,
        mcp_config: Optional[MCPConfig] = None,
        retry_controller: Optional[RetryController] = None
    ):
        """
        Initialize unified UniProt client.
        
        Args:
            api_config: API configuration, uses global config if not provided
            mcp_config: MCP configuration, uses global config if not provided
            retry_controller: Retry controller, uses global if not provided
        """
        config = get_config()
        self.api_config = api_config or config.api
        self.mcp_config = mcp_config or config.mcp
        self.retry_controller = retry_controller or get_retry_controller()
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        
        # Initialize clients
        self.mcp_client = MCPServerClient(self.mcp_config, self.retry_controller)
        self.rest_client = UniProtRESTClient(self.api_config, self.retry_controller)
        
        # Track which method is being used
        self._using_mcp = None
    
    async def __aenter__(self):
        """Async context manager entry."""
        await self.rest_client.__aenter__()
        return self
    
    async def __aexit__(self, exc_type, exc_val, exc_tb):
        """Async context manager exit."""
        await self.rest_client.__aexit__(exc_type, exc_val, exc_tb)
    
    async def _determine_access_method(self) -> str:
        """
        Determine which access method to use (MCP or REST).
        
        Returns:
            "mcp" if MCP server is available, "rest" otherwise
        """
        if self._using_mcp is not None:
            return "mcp" if self._using_mcp else "rest"
        
        # Check MCP server availability
        if self.mcp_config.enabled and await self.mcp_client.is_server_available():
            self._using_mcp = True
            self.logger.info("Using MCP server for UniProt data access")
            return "mcp"
        else:
            self._using_mcp = False
            if self.mcp_config.fallback_to_rest:
                self.logger.info("MCP server unavailable, falling back to REST API")
                return "rest"
            else:
                raise ConfigurationError("MCP server unavailable and fallback disabled")
    
    async def get_protein_data(self, uniprot_id: str) -> Optional[Dict[str, Any]]:
        """
        Get protein data using the best available method.
        
        Args:
            uniprot_id: UniProt protein identifier
            
        Returns:
            Protein data dictionary or None if not found
            
        Raises:
            APIError: For API-specific errors
            NetworkError: For network-related errors
            ConfigurationError: If no access method is available
        """
        async def operation():
            method = await self._determine_access_method()
            
            if method == "mcp":
                # Try MCP server first
                try:
                    result = await self.mcp_client.get_protein_data(uniprot_id)
                    if result is not None:
                        return result
                    
                    # If MCP returns None but fallback is enabled, try REST
                    if self.mcp_config.fallback_to_rest:
                        self.logger.debug(
                            "MCP server returned no data, trying REST API fallback",
                            extra={"uniprot_id": uniprot_id}
                        )
                        return await self.rest_client.get_protein_data(uniprot_id)
                    
                    return None
                    
                except Exception as e:
                    # If MCP fails and fallback is enabled, try REST
                    if self.mcp_config.fallback_to_rest:
                        self.logger.warning(
                            "MCP server failed, falling back to REST API: %s",
                            str(e),
                            extra={"uniprot_id": uniprot_id}
                        )
                        return await self.rest_client.get_protein_data(uniprot_id)
                    raise
            else:
                # Use REST API directly
                return await self.rest_client.get_protein_data(uniprot_id)
        
        return await self.retry_controller.execute_with_retry_async(
            operation,
            database="UniProt",
            operation_name=f"get_protein_data_{uniprot_id}"
        )
    
    async def get_protein_isoforms(self, uniprot_id: str) -> List[Dict[str, Any]]:
        """
        Get all isoforms for a protein using the best available method.
        
        Args:
            uniprot_id: UniProt protein identifier
            
        Returns:
            List of isoform data dictionaries
            
        Raises:
            APIError: For API-specific errors
            NetworkError: For network-related errors
            ConfigurationError: If no access method is available
        """
        async def operation():
            method = await self._determine_access_method()
            
            if method == "mcp":
                # Try MCP server first
                try:
                    result = await self.mcp_client.get_protein_isoforms(uniprot_id)
                    if result:
                        return result
                    
                    # If MCP returns empty but fallback is enabled, try REST
                    if self.mcp_config.fallback_to_rest:
                        self.logger.debug(
                            "MCP server returned no isoforms, trying REST API fallback",
                            extra={"uniprot_id": uniprot_id}
                        )
                        return await self.rest_client.get_protein_isoforms(uniprot_id)
                    
                    return []
                    
                except Exception as e:
                    # If MCP fails and fallback is enabled, try REST
                    if self.mcp_config.fallback_to_rest:
                        self.logger.warning(
                            "MCP server failed for isoforms, falling back to REST API: %s",
                            str(e),
                            extra={"uniprot_id": uniprot_id}
                        )
                        return await self.rest_client.get_protein_isoforms(uniprot_id)
                    raise
            else:
                # Use REST API directly
                return await self.rest_client.get_protein_isoforms(uniprot_id)
        
        return await self.retry_controller.execute_with_retry_async(
            operation,
            database="UniProt",
            operation_name=f"get_protein_isoforms_{uniprot_id}"
        )
    
    def parse_protein_isoform_data(
        self,
        isoform_data: Dict[str, Any],
        parent_protein_id: str,
        method: str = "unknown"
    ) -> ProteinModel:
        """
        Parse isoform data into ProteinModel regardless of source (MCP or REST).
        
        Args:
            isoform_data: Raw isoform data from MCP or REST API
            parent_protein_id: UniProt ID of the parent protein
            method: Source method ("mcp" or "rest")
            
        Returns:
            Validated ProteinModel instance
            
        Raises:
            DataError: For invalid or incomplete data
        """
        try:
            # Extract common fields with fallbacks for different data structures
            isoform_id = (
                isoform_data.get("isoform_id") or
                isoform_data.get("name", {}).get("value") or
                f"{parent_protein_id}-1"
            )
            
            sequence = (
                isoform_data.get("sequence") or
                isoform_data.get("sequence", {}).get("value", "")
            )
            
            sequence_length = (
                isoform_data.get("length") or
                isoform_data.get("sequence", {}).get("length") or
                len(sequence)
            )
            
            # Extract exon annotations (structure varies by source)
            exon_annotations = {}
            if "exon_annotations" in isoform_data:
                exon_annotations = isoform_data["exon_annotations"]
            elif "features" in isoform_data:
                # Extract exon information from features
                exons = [f for f in isoform_data["features"] if f.get("type") == "exon"]
                if exons:
                    exon_annotations = {"exons": exons}
            
            # Calculate exon count
            exon_count = 0
            if isinstance(exon_annotations, dict):
                if "exons" in exon_annotations:
                    exon_count = len(exon_annotations["exons"])
                elif "exon_count" in exon_annotations:
                    exon_count = exon_annotations["exon_count"]
            
            # Extract TIM barrel location information
            tim_barrel_location = {}
            if "tim_barrel_location" in isoform_data:
                tim_barrel_location = isoform_data["tim_barrel_location"]
            elif "domains" in isoform_data:
                # Look for TIM barrel domains
                for domain in isoform_data["domains"]:
                    if "tim" in domain.get("name", "").lower() or "barrel" in domain.get("name", "").lower():
                        tim_barrel_location = domain
                        break
            
            # Extract organism and names
            organism = (
                isoform_data.get("organism") or
                isoform_data.get("organism_name") or
                "Homo sapiens"
            )
            
            name = (
                isoform_data.get("name") or
                isoform_data.get("protein_name") or
                isoform_data.get("description", "")
            )
            
            description = (
                isoform_data.get("description") or
                isoform_data.get("note", {}).get("texts", [{}])[0].get("value", "") or
                f"Isoform of {parent_protein_id}"
            )
            
            return ProteinModel(
                isoform_id=isoform_id,
                parent_protein_id=parent_protein_id,
                sequence=sequence,
                sequence_length=sequence_length,
                exon_annotations=exon_annotations,
                exon_count=exon_count,
                tim_barrel_location=tim_barrel_location,
                organism=organism,
                name=name,
                description=description
            )
            
        except Exception as e:
            raise DataError(f"Failed to parse protein isoform data from {method}: {str(e)}")
    
    async def get_access_method_status(self) -> Dict[str, Any]:
        """
        Get status information about available access methods.
        
        Returns:
            Dictionary with status information
        """
        mcp_available = await self.mcp_client.is_server_available()
        
        return {
            "mcp_enabled": self.mcp_config.enabled,
            "mcp_available": mcp_available,
            "fallback_enabled": self.mcp_config.fallback_to_rest,
            "current_method": "mcp" if mcp_available and self.mcp_config.enabled else "rest",
            "server_path": self.mcp_config.server_path
        }
    
    async def get_rate_limiting_stats(self) -> Dict[str, Any]:
        """
        Get comprehensive rate limiting statistics for all APIs.
        
        Returns:
            Dictionary with rate limiting statistics and reports
        """
        rate_limit_manager = get_rate_limit_manager()
        
        return {
            "all_stats": rate_limit_manager.get_all_stats(),
            "recent_violations": rate_limit_manager.get_recent_violations(60),
            "comprehensive_report": rate_limit_manager.generate_report()
        }
    
    async def get_performance_metrics(self) -> Dict[str, Any]:
        """
        Get comprehensive performance metrics for all UniProt clients.
        
        Returns:
            Dictionary with performance metrics
        """
        cache_metrics = get_global_cache().get_metrics()
        rate_limiting_stats = await self.get_rate_limiting_stats()
        access_method_status = await self.get_access_method_status()
        
        return {
            "cache_metrics": cache_metrics,
            "rate_limiting": rate_limiting_stats,
            "access_method": access_method_status,
            "mcp_client_stats": {
                "enabled": self.mcp_config.enabled,
                "server_path": self.mcp_config.server_path,
                "fallback_enabled": self.mcp_config.fallback_to_rest
            },
            "rest_client_stats": {
                "base_url": self.api_config.uniprot_base_url,
                "connection_timeout": self.api_config.connection_timeout,
                "request_timeout": self.api_config.request_timeout
            }
        }
    
    async def close(self) -> None:
        """Close all clients."""
        await self.rest_client.close()


# Convenience functions for common operations
async def get_protein_with_isoforms(
    uniprot_id: str,
    client: Optional[UnifiedUniProtClient] = None
) -> List[ProteinModel]:
    """
    Get all isoforms for a protein using the unified client.
    
    Args:
        uniprot_id: UniProt protein identifier
        client: Unified client, creates new one if not provided
        
    Returns:
        List of validated ProteinModel instances
    """
    if client is None:
        async with UnifiedUniProtClient() as client:
            return await get_protein_with_isoforms(uniprot_id, client)
    
    try:
        isoforms_data = await client.get_protein_isoforms(uniprot_id)
        proteins = []
        
        method = await client._determine_access_method()
        
        for isoform_data in isoforms_data:
            try:
                protein = client.parse_protein_isoform_data(isoform_data, uniprot_id, method)
                proteins.append(protein)
            except DataError as e:
                client.logger.warning(
                    "Skipping invalid isoform data for protein %s: %s",
                    uniprot_id,
                    str(e),
                    extra={"uniprot_id": uniprot_id, "isoform_data": isoform_data}
                )
        
        return proteins
        
    except Exception as e:
        client.logger.error(
            "Failed to get protein isoforms for %s: %s",
            uniprot_id,
            str(e),
            extra={"uniprot_id": uniprot_id}
        )
        return []


async def get_proteins_batch(
    uniprot_ids: List[str],
    client: Optional[UnifiedUniProtClient] = None,
    batch_size: int = 10
) -> List[ProteinModel]:
    """
    Get isoforms for multiple proteins in batches.
    
    Args:
        uniprot_ids: List of UniProt protein identifiers
        client: Unified client, creates new one if not provided
        batch_size: Number of proteins to process concurrently
        
    Returns:
        List of all ProteinModel instances
    """
    if client is None:
        async with UnifiedUniProtClient() as client:
            return await get_proteins_batch(uniprot_ids, client, batch_size)
    
    all_proteins = []
    
    # Process proteins in batches to avoid overwhelming the APIs
    for i in range(0, len(uniprot_ids), batch_size):
        batch = uniprot_ids[i:i + batch_size]
        
        # Create tasks for concurrent processing
        tasks = [get_protein_with_isoforms(uniprot_id, client) for uniprot_id in batch]
        
        # Execute batch concurrently
        batch_results = await asyncio.gather(*tasks, return_exceptions=True)
        
        # Collect results and handle exceptions
        for uniprot_id, result in zip(batch, batch_results):
            if isinstance(result, Exception):
                client.logger.error(
                    "Failed to process protein %s: %s",
                    uniprot_id,
                    str(result),
                    extra={"uniprot_id": uniprot_id}
                )
            else:
                all_proteins.extend(result)
    
    return all_proteins