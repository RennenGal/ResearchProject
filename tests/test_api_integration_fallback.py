"""
Property-based tests for API integration fallback functionality.

Tests Property 9: API Integration Fallback
Validates Requirements 8.1, 8.2, 8.4
"""

import pytest
import asyncio
from unittest.mock import AsyncMock, MagicMock, patch
from hypothesis import given, strategies as st, settings, assume
from typing import Dict, Any, List

from protein_data_collector.api.uniprot_client import UnifiedUniProtClient, MCPServerClient, UniProtRESTClient
from protein_data_collector.config import APIConfig, MCPConfig, RetryConfig, RateLimitingConfig
from protein_data_collector.errors import APIError, NetworkError, ConfigurationError


# Test data generators
@st.composite
def uniprot_id_strategy(draw):
    """Generate valid UniProt IDs for testing."""
    # UniProt IDs are typically 6-10 characters, alphanumeric
    length = draw(st.integers(min_value=6, max_value=10))
    return draw(st.text(alphabet=st.characters(whitelist_categories=('Lu', 'Nd')), min_size=length, max_size=length))


@st.composite
def protein_data_strategy(draw):
    """Generate protein data dictionaries for testing."""
    uniprot_id = draw(uniprot_id_strategy())
    return {
        "accession": uniprot_id,
        "id": f"{uniprot_id}_HUMAN",
        "protein_name": draw(st.text(min_size=5, max_size=100)),
        "organism_name": "Homo sapiens",
        "sequence": {
            "value": draw(st.text(alphabet="ACDEFGHIKLMNPQRSTVWY", min_size=50, max_size=500)),
            "length": draw(st.integers(min_value=50, max_value=500))
        },
        "gene_names": [{"value": draw(st.text(min_size=3, max_size=10))}],
        "comments": []
    }


@st.composite
def isoform_data_strategy(draw):
    """Generate isoform data for testing."""
    uniprot_id = draw(uniprot_id_strategy())
    isoform_num = draw(st.integers(min_value=1, max_value=5))
    sequence = draw(st.text(alphabet="ACDEFGHIKLMNPQRSTVWY", min_size=50, max_size=500))
    
    return {
        "isoform_id": f"{uniprot_id}-{isoform_num}",
        "sequence": sequence,
        "length": len(sequence),
        "is_canonical": isoform_num == 1,
        "description": draw(st.text(min_size=10, max_size=200))
    }


class TestAPIIntegrationFallback:
    """Test suite for API integration fallback functionality."""
    
    def _create_api_config(self):
        """Create API configuration for testing."""
        return APIConfig(
            uniprot_base_url="https://rest.uniprot.org/",
            request_timeout=30,
            connection_timeout=10
        )
    
    def _create_rate_limiting_config(self):
        """Create rate limiting configuration for testing."""
        return RateLimitingConfig(
            uniprot_requests_per_second=5.0
        )
    
    def _create_mcp_config_enabled(self):
        """Create MCP configuration with fallback enabled."""
        return MCPConfig(
            enabled=True,
            server_path="/mock/mcp/server",
            fallback_to_rest=True
        )
    
    def _create_mcp_config_disabled(self):
        """Create MCP configuration disabled."""
        return MCPConfig(
            enabled=False,
            server_path="/mock/mcp/server",
            fallback_to_rest=True
        )
    
    def _create_retry_config(self):
        """Create retry configuration for testing."""
        return RetryConfig(
            max_retries=2,
            initial_delay=0.1,
            backoff_multiplier=1.5,
            max_delay=5.0
        )
    
    def _create_retry_controller(self):
        """Create retry controller for testing."""
        from protein_data_collector.retry import RetryController
        return RetryController(self._create_retry_config())
    
    @given(uniprot_id=uniprot_id_strategy(), protein_data=protein_data_strategy())
    @settings(max_examples=100, deadline=5000)
    @pytest.mark.asyncio
    async def test_mcp_available_uses_mcp_server(self, uniprot_id, protein_data):
        """
        Property 9: API Integration Fallback - MCP server used when available.
        
        For any UniProt data request, when MCP servers are available and functional,
        they should be used.
        
        **Validates: Requirements 8.1, 8.2, 8.4**
        """
        assume(len(uniprot_id) >= 6)
        assume(protein_data["sequence"]["value"])
        
        api_config = self._create_api_config()
        mcp_config_enabled = self._create_mcp_config_enabled()
        retry_controller = self._create_retry_controller()
        
        with patch('protein_data_collector.api.uniprot_client.MCPServerClient') as mock_mcp_class, \
             patch('protein_data_collector.api.uniprot_client.UniProtRESTClient') as mock_rest_class:
            
            # Setup MCP server mock to be available and return data
            mock_mcp = AsyncMock()
            mock_mcp.is_server_available.return_value = True
            mock_mcp.get_protein_data.return_value = protein_data
            mock_mcp_class.return_value = mock_mcp
            
            # Setup REST client mock (should not be called)
            mock_rest = AsyncMock()
            mock_rest_class.return_value = mock_rest
            
            # Create unified client
            async with UnifiedUniProtClient(api_config, mcp_config_enabled, retry_controller) as client:
                # Get protein data
                result = await client.get_protein_data(uniprot_id)
                
                # Verify MCP server was used
                mock_mcp.is_server_available.assert_called()
                mock_mcp.get_protein_data.assert_called_once_with(uniprot_id)
                
                # Verify REST client was not used for the main call
                mock_rest.get_protein_data.assert_not_called()
                
                # Verify result matches expected data
                assert result == protein_data
    
    @given(uniprot_id=uniprot_id_strategy(), protein_data=protein_data_strategy())
    @settings(max_examples=100, deadline=5000)
    @pytest.mark.asyncio
    async def test_mcp_unavailable_falls_back_to_rest(self, uniprot_id, protein_data):
        """
        Property 9: API Integration Fallback - REST API used when MCP unavailable.
        
        For any UniProt data request, when MCP servers are unavailable,
        the system should seamlessly fall back to direct REST API calls.
        
        **Validates: Requirements 8.1, 8.2, 8.4**
        """
        assume(len(uniprot_id) >= 6)
        assume(protein_data["sequence"]["value"])
        
        api_config = self._create_api_config()
        mcp_config_enabled = self._create_mcp_config_enabled()
        retry_controller = self._create_retry_controller()
        
        with patch('protein_data_collector.api.uniprot_client.MCPServerClient') as mock_mcp_class, \
             patch('protein_data_collector.api.uniprot_client.UniProtRESTClient') as mock_rest_class:
            
            # Setup MCP server mock to be unavailable
            mock_mcp = AsyncMock()
            mock_mcp.is_server_available.return_value = False
            mock_mcp_class.return_value = mock_mcp
            
            # Setup REST client mock to return data
            mock_rest = AsyncMock()
            mock_rest.get_protein_data.return_value = protein_data
            mock_rest_class.return_value = mock_rest
            
            # Create unified client
            async with UnifiedUniProtClient(api_config, mcp_config_enabled, retry_controller) as client:
                # Get protein data
                result = await client.get_protein_data(uniprot_id)
                
                # Verify MCP server availability was checked
                mock_mcp.is_server_available.assert_called()
                
                # Verify REST client was used as fallback
                mock_rest.get_protein_data.assert_called_once_with(uniprot_id)
                
                # Verify result matches expected data
                assert result == protein_data
    
    @given(uniprot_id=uniprot_id_strategy(), isoforms=st.lists(isoform_data_strategy(), min_size=1, max_size=5))
    @settings(max_examples=100, deadline=5000)
    @pytest.mark.asyncio
    async def test_mcp_failure_falls_back_to_rest_for_isoforms(self, uniprot_id, isoforms):
        """
        Property 9: API Integration Fallback - REST fallback on MCP failure for isoforms.
        
        For any protein isoform request, when MCP servers fail during operation,
        the system should fall back to REST API and return identical results.
        
        **Validates: Requirements 8.1, 8.2, 8.4**
        """
        assume(len(uniprot_id) >= 6)
        assume(all(isoform["sequence"] for isoform in isoforms))
        
        api_config = self._create_api_config()
        mcp_config_enabled = self._create_mcp_config_enabled()
        retry_controller = self._create_retry_controller()
        
        with patch('protein_data_collector.api.uniprot_client.MCPServerClient') as mock_mcp_class, \
             patch('protein_data_collector.api.uniprot_client.UniProtRESTClient') as mock_rest_class:
            
            # Setup MCP server mock to be available but fail during operation
            mock_mcp = AsyncMock()
            mock_mcp.is_server_available.return_value = True
            mock_mcp.get_protein_isoforms.side_effect = APIError("MCP server error")
            mock_mcp_class.return_value = mock_mcp
            
            # Setup REST client mock to return isoforms
            mock_rest = AsyncMock()
            mock_rest.get_protein_isoforms.return_value = isoforms
            mock_rest_class.return_value = mock_rest
            
            # Create unified client
            async with UnifiedUniProtClient(api_config, mcp_config_enabled, retry_controller) as client:
                # Get protein isoforms
                result = await client.get_protein_isoforms(uniprot_id)
                
                # Verify MCP server was tried first
                mock_mcp.is_server_available.assert_called()
                mock_mcp.get_protein_isoforms.assert_called_once_with(uniprot_id)
                
                # Verify REST client was used as fallback
                mock_rest.get_protein_isoforms.assert_called_once_with(uniprot_id)
                
                # Verify result matches expected data
                assert result == isoforms
    
    @given(uniprot_id=uniprot_id_strategy())
    @settings(max_examples=100, deadline=5000)
    @pytest.mark.asyncio
    async def test_mcp_disabled_uses_rest_directly(self, uniprot_id):
        """
        Property 9: API Integration Fallback - Direct REST when MCP disabled.
        
        For any UniProt data request, when MCP servers are disabled in configuration,
        the system should use REST API directly without checking MCP availability.
        
        **Validates: Requirements 8.1, 8.2, 8.4**
        """
        assume(len(uniprot_id) >= 6)
        
        api_config = self._create_api_config()
        mcp_config_disabled = self._create_mcp_config_disabled()
        retry_controller = self._create_retry_controller()
        
        protein_data = {
            "accession": uniprot_id,
            "sequence": {"value": "ACDEFGHIKLMNPQRSTVWY", "length": 20}
        }
        
        with patch('protein_data_collector.api.uniprot_client.MCPServerClient') as mock_mcp_class, \
             patch('protein_data_collector.api.uniprot_client.UniProtRESTClient') as mock_rest_class:
            
            # Setup MCP server mock (should not be used)
            mock_mcp = AsyncMock()
            mock_mcp_class.return_value = mock_mcp
            
            # Setup REST client mock to return data
            mock_rest = AsyncMock()
            mock_rest.get_protein_data.return_value = protein_data
            mock_rest_class.return_value = mock_rest
            
            # Create unified client with MCP disabled
            async with UnifiedUniProtClient(api_config, mcp_config_disabled, retry_controller) as client:
                # Get protein data
                result = await client.get_protein_data(uniprot_id)
                
                # Verify MCP server availability was not checked (since disabled)
                mock_mcp.is_server_available.assert_not_called()
                
                # Verify REST client was used directly
                mock_rest.get_protein_data.assert_called_once_with(uniprot_id)
                
                # Verify result matches expected data
                assert result == protein_data
    
    @given(uniprot_id=uniprot_id_strategy(), protein_data=protein_data_strategy())
    @settings(max_examples=100, deadline=5000)
    @pytest.mark.asyncio
    async def test_identical_results_from_mcp_and_rest(self, uniprot_id, protein_data):
        """
        Property 9: API Integration Fallback - Identical results from both methods.
        
        For any UniProt data request, the results from MCP server and REST API
        should be functionally equivalent when processed by the unified client.
        
        **Validates: Requirements 8.1, 8.2, 8.4**
        """
        assume(len(uniprot_id) >= 6)
        assume(protein_data["sequence"]["value"])
        
        api_config = self._create_api_config()
        mcp_config_enabled = self._create_mcp_config_enabled()
        retry_controller = self._create_retry_controller()
        
        with patch('protein_data_collector.api.uniprot_client.MCPServerClient') as mock_mcp_class, \
             patch('protein_data_collector.api.uniprot_client.UniProtRESTClient') as mock_rest_class:
            
            # Test MCP path
            mock_mcp = AsyncMock()
            mock_mcp.is_server_available.return_value = True
            mock_mcp.get_protein_data.return_value = protein_data
            mock_mcp_class.return_value = mock_mcp
            
            mock_rest = AsyncMock()
            mock_rest_class.return_value = mock_rest
            
            async with UnifiedUniProtClient(api_config, mcp_config_enabled, retry_controller) as client:
                mcp_result = await client.get_protein_data(uniprot_id)
            
            # Reset mocks for REST path
            mock_mcp.reset_mock()
            mock_rest.reset_mock()
            
            # Test REST path
            mock_mcp.is_server_available.return_value = False
            mock_rest.get_protein_data.return_value = protein_data
            
            async with UnifiedUniProtClient(api_config, mcp_config_enabled, retry_controller) as client:
                rest_result = await client.get_protein_data(uniprot_id)
            
            # Verify both methods return identical results
            assert mcp_result == rest_result == protein_data
    
    @given(uniprot_id=uniprot_id_strategy())
    @settings(max_examples=100, deadline=5000)
    @pytest.mark.asyncio
    async def test_fallback_disabled_raises_error_when_mcp_unavailable(self, uniprot_id):
        """
        Property 9: API Integration Fallback - Error when fallback disabled.
        
        For any UniProt data request, when MCP servers are unavailable and
        fallback is disabled, the system should raise a configuration error.
        
        **Validates: Requirements 8.1, 8.2, 8.4**
        """
        assume(len(uniprot_id) >= 6)
        
        api_config = self._create_api_config()
        retry_controller = self._create_retry_controller()
        
        # Create MCP config with fallback disabled
        mcp_config_no_fallback = MCPConfig(
            enabled=True,
            server_path="/mock/mcp/server",
            fallback_to_rest=False
        )
        
        with patch('protein_data_collector.api.uniprot_client.MCPServerClient') as mock_mcp_class:
            # Setup MCP server mock to be unavailable
            mock_mcp = AsyncMock()
            mock_mcp.is_server_available.return_value = False
            mock_mcp_class.return_value = mock_mcp
            
            # Create unified client with fallback disabled
            async with UnifiedUniProtClient(api_config, mcp_config_no_fallback, retry_controller) as client:
                # Attempt to get protein data should raise ConfigurationError
                with pytest.raises(ConfigurationError, match="MCP server unavailable and fallback disabled"):
                    await client.get_protein_data(uniprot_id)
    
    @given(uniprot_id=uniprot_id_strategy())
    @settings(max_examples=100, deadline=5000)
    @pytest.mark.asyncio
    async def test_access_method_status_reflects_current_state(self, uniprot_id):
        """
        Property 9: API Integration Fallback - Status reflects actual access method.
        
        For any system state, the access method status should accurately reflect
        which method (MCP or REST) is currently being used.
        
        **Validates: Requirements 8.1, 8.2, 8.4**
        """
        assume(len(uniprot_id) >= 6)
        
        api_config = self._create_api_config()
        mcp_config_enabled = self._create_mcp_config_enabled()
        retry_controller = self._create_retry_controller()
        
        with patch('protein_data_collector.api.uniprot_client.MCPServerClient') as mock_mcp_class:
            # Test with MCP available
            mock_mcp = AsyncMock()
            mock_mcp.is_server_available.return_value = True
            mock_mcp_class.return_value = mock_mcp
            
            async with UnifiedUniProtClient(api_config, mcp_config_enabled, retry_controller) as client:
                status = await client.get_access_method_status()
                
                assert status["mcp_enabled"] == True
                assert status["mcp_available"] == True
                assert status["fallback_enabled"] == True
                assert status["current_method"] == "mcp"
            
            # Test with MCP unavailable
            mock_mcp.is_server_available.return_value = False
            
            async with UnifiedUniProtClient(api_config, mcp_config_enabled, retry_controller) as client:
                status = await client.get_access_method_status()
                
                assert status["mcp_enabled"] == True
                assert status["mcp_available"] == False
                assert status["fallback_enabled"] == True
                assert status["current_method"] == "rest"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])