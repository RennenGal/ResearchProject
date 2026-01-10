"""
Property-based tests for complete isoform data collection functionality.

Tests Property 4: Complete Isoform Data Collection
Validates Requirements 3.2, 3.3, 3.4, 3.5, 3.6
"""

import pytest
import asyncio
from unittest.mock import AsyncMock, MagicMock, patch
from hypothesis import given, strategies as st, settings, assume
from typing import Dict, Any, List

from protein_data_collector.collector.uniprot_collector import UniProtIsoformCollector, collect_all_isoforms
from protein_data_collector.models.entities import InterProProteinModel, ProteinModel
from protein_data_collector.api.uniprot_client import UnifiedUniProtClient
from protein_data_collector.config import APIConfig, RetryConfig
from protein_data_collector.retry import RetryController


# Test data generators
@st.composite
def uniprot_id_strategy(draw):
    """Generate valid UniProt IDs for testing."""
    # Use integers to ensure uniqueness, then convert to string format
    unique_num = draw(st.integers(min_value=100000, max_value=999999))
    return f"P{unique_num}"


@st.composite
def interpro_protein_strategy(draw):
    """Generate InterProProteinModel instances for testing."""
    uniprot_id = draw(uniprot_id_strategy())
    tim_barrel_accession = draw(st.text(min_size=5, max_size=15))
    name = draw(st.text(min_size=5, max_size=100))
    
    return InterProProteinModel(
        uniprot_id=uniprot_id,
        tim_barrel_accession=tim_barrel_accession,
        name=name,
        organism="Homo sapiens",
        basic_metadata={}
    )


@st.composite
def complete_isoform_data_strategy(draw):
    """Generate complete isoform data with all required fields."""
    uniprot_id = draw(uniprot_id_strategy())
    isoform_num = draw(st.integers(min_value=1, max_value=5))
    sequence = draw(st.text(alphabet="ACDEFGHIKLMNPQRSTVWY", min_size=50, max_size=500))
    
    # Generate exon annotations
    exon_count = draw(st.integers(min_value=1, max_value=10))
    exons = []
    current_pos = 1
    for i in range(exon_count):
        exon_length = draw(st.integers(min_value=10, max_value=50))
        exons.append({
            "start": current_pos,
            "end": current_pos + exon_length - 1,
            "type": "exon"
        })
        current_pos += exon_length
    
    # Generate TIM barrel location within sequence bounds
    tim_start = draw(st.integers(min_value=1, max_value=max(1, len(sequence) - 50)))
    tim_end = draw(st.integers(min_value=tim_start + 10, max_value=len(sequence)))
    
    return {
        "isoform_id": f"{uniprot_id}-{isoform_num}",
        "sequence": sequence,
        "length": len(sequence),
        "exon_annotations": {"exons": exons},
        "exon_count": exon_count,
        "tim_barrel_location": {
            "start": tim_start,
            "end": tim_end,
            "confidence": draw(st.floats(min_value=0.5, max_value=1.0))
        },
        "organism": "Homo sapiens",
        "name": draw(st.text(min_size=5, max_size=100)),
        "description": draw(st.text(min_size=10, max_size=200))
    }


class TestCompleteIsoformDataCollection:
    """Test suite for complete isoform data collection functionality."""
    
    def _create_api_config(self):
        """Create API configuration for testing."""
        return APIConfig(
            uniprot_base_url="https://rest.uniprot.org/",
            request_timeout=30,
            connection_timeout=10
        )
    
    def _create_retry_controller(self):
        """Create retry controller for testing."""
        retry_config = RetryConfig(
            max_retries=2,
            initial_delay=0.1,
            backoff_multiplier=1.5,
            max_delay=5.0
        )
        return RetryController(retry_config)
    
    @given(
        interpro_protein=interpro_protein_strategy(),
        isoforms_data=st.lists(complete_isoform_data_strategy(), min_size=1, max_size=5)
    )
    @settings(max_examples=100, deadline=5000)
    @pytest.mark.asyncio
    async def test_complete_isoform_data_collection_all_fields(self, interpro_protein, isoforms_data):
        """
        Property 4: Complete Isoform Data Collection - All required fields collected.
        
        For any protein isoform retrieved from UniProt, the system should collect
        all required fields: sequence information, exon annotations, TIM barrel location,
        protein name, description, and organism information.
        
        **Validates: Requirements 3.2, 3.3, 3.4, 3.5, 3.6**
        """
        assume(len(interpro_protein.uniprot_id) >= 6)
        assume(all(isoform["sequence"] for isoform in isoforms_data))
        assume(all(isoform["exon_annotations"] for isoform in isoforms_data))
        
        api_config = self._create_api_config()
        retry_controller = self._create_retry_controller()
        
        with patch('protein_data_collector.api.uniprot_client.UnifiedUniProtClient') as mock_client_class:
            # Setup mock client
            mock_client = AsyncMock()
            mock_client.get_protein_isoforms.return_value = isoforms_data
            mock_client._determine_access_method.return_value = "rest"
            
            # Mock the parse method to return valid ProteinModel instances
            def mock_parse_isoform_data(isoform_data, parent_id, method):
                return ProteinModel(
                    isoform_id=isoform_data["isoform_id"],
                    parent_protein_id=parent_id,
                    parent_tim_barrel_accession="PF00001",  # Default test value
                    sequence=isoform_data["sequence"],
                    sequence_length=isoform_data["length"],
                    exon_annotations=isoform_data["exon_annotations"],
                    exon_count=isoform_data["exon_count"],
                    tim_barrel_location=isoform_data["tim_barrel_location"],
                    organism_name=isoform_data["organism"],  # Use organism_name instead of organism
                    protein_name=isoform_data["name"],  # Use protein_name instead of name
                    # Note: No description field in ProteinModel
                )
            
            # Set the mock method directly instead of using side_effect
            mock_client.parse_protein_isoform_data = mock_parse_isoform_data
            mock_client_class.return_value = mock_client
            
            # Test the collector
            async with UniProtIsoformCollector(mock_client) as collector:
                collector.reset_collection_state()
                result_proteins = await collector.collect_protein_isoforms(interpro_protein)
                
                # Verify all isoforms were collected
                assert len(result_proteins) == len(isoforms_data)
                
                # Verify all required fields are present for each isoform
                for protein in result_proteins:
                    # Requirement 3.2: Sequence information
                    assert protein.sequence is not None
                    assert len(protein.sequence) > 0
                    assert protein.sequence_length > 0
                    assert protein.sequence_length == len(protein.sequence)
                    
                    # Requirement 3.3: Exon annotation data
                    assert protein.exon_annotations is not None
                    assert isinstance(protein.exon_annotations, dict)
                    assert protein.exon_count is not None
                    assert protein.exon_count >= 0
                    
                    # Requirement 3.4: TIM barrel location information
                    assert protein.tim_barrel_location is not None
                    assert isinstance(protein.tim_barrel_location, dict)
                    
                    # Requirement 3.5: Protein name and description
                    assert protein.protein_name is not None
                    assert len(protein.protein_name) > 0
                    # Note: ProteinModel doesn't have a description field, it uses various comment fields
                    
                    # Requirement 3.6: Organism information
                    assert protein.organism_name is not None
                    assert protein.organism_name == "Homo sapiens"
                    
                    # Verify parent relationship
                    assert protein.parent_protein_id == interpro_protein.uniprot_id
    
    @given(
        interpro_proteins=st.lists(interpro_protein_strategy(), min_size=1, max_size=10, unique_by=lambda x: x.uniprot_id),
        isoforms_per_protein=st.integers(min_value=1, max_value=3)
    )
    @settings(max_examples=100, deadline=5000)
    @pytest.mark.asyncio
    async def test_batch_collection_completeness(self, interpro_proteins, isoforms_per_protein):
        """
        Property 4: Complete Isoform Data Collection - Batch processing completeness.
        
        For any batch of proteins, the system should collect complete isoform data
        for each protein, maintaining data integrity across the entire batch.
        
        **Validates: Requirements 3.2, 3.3, 3.4, 3.5, 3.6**
        """
        assume(len(interpro_proteins) >= 1)
        assume(all(len(p.uniprot_id) >= 6 for p in interpro_proteins))
        
        # Generate isoform data for each protein
        all_isoforms_data = {}
        for protein in interpro_proteins:
            protein_isoforms = []
            for i in range(isoforms_per_protein):
                isoform_data = {
                    "isoform_id": f"{protein.uniprot_id}-{i+1}",
                    "sequence": "ACDEFGHIKLMNPQRSTVWY" * 10,  # 200 AA sequence
                    "length": 200,
                    "exon_annotations": {"exons": [{"start": 1, "end": 200, "type": "exon"}]},
                    "exon_count": 1,
                    "tim_barrel_location": {"start": 10, "end": 190, "confidence": 0.9},
                    "organism": "Homo sapiens",
                    "name": f"Test protein {protein.uniprot_id}",
                    "description": f"Test description for {protein.uniprot_id}"
                }
                protein_isoforms.append(isoform_data)
            all_isoforms_data[protein.uniprot_id] = protein_isoforms
        
        with patch('protein_data_collector.api.uniprot_client.UnifiedUniProtClient') as mock_client_class:
            # Setup mock client
            mock_client = AsyncMock()
            
            def mock_get_isoforms(uniprot_id):
                return all_isoforms_data.get(uniprot_id, [])
            
            mock_client.get_protein_isoforms.side_effect = mock_get_isoforms
            mock_client._determine_access_method.return_value = "rest"
            
            # Mock the parse method
            def mock_parse_isoform_data(isoform_data, parent_id, method):
                return ProteinModel(
                    isoform_id=isoform_data["isoform_id"],
                    parent_protein_id=parent_id,
            parent_tim_barrel_accession="PF00001",  # Default test value
                    sequence=isoform_data["sequence"],
                    sequence_length=isoform_data["length"],
                    exon_annotations=isoform_data["exon_annotations"],
                    exon_count=isoform_data["exon_count"],
                    tim_barrel_location=isoform_data["tim_barrel_location"],
                    organism_name=isoform_data["organism"],
                    protein_name=isoform_data["name"],
                )
            
            # Set the mock method directly
            mock_client.parse_protein_isoform_data = mock_parse_isoform_data
            mock_client_class.return_value = mock_client
            
            # Test batch collection
            async with UniProtIsoformCollector(mock_client) as collector:
                collector.reset_collection_state()
                result_proteins = await collector.collect_isoforms_batch(interpro_proteins)
                
                # Verify total count
                expected_total = len(interpro_proteins) * isoforms_per_protein
                assert len(result_proteins) == expected_total
                
                # Group results by parent protein
                proteins_by_parent = {}
                for protein in result_proteins:
                    parent_id = protein.parent_protein_id
                    if parent_id not in proteins_by_parent:
                        proteins_by_parent[parent_id] = []
                    proteins_by_parent[parent_id].append(protein)
                
                # Verify each protein has the expected number of isoforms
                for interpro_protein in interpro_proteins:
                    parent_id = interpro_protein.uniprot_id
                    assert parent_id in proteins_by_parent
                    assert len(proteins_by_parent[parent_id]) == isoforms_per_protein
                    
                    # Verify completeness for each isoform
                    for protein in proteins_by_parent[parent_id]:
                        assert protein.sequence is not None and len(protein.sequence) > 0
                        assert protein.exon_annotations is not None
                        assert protein.tim_barrel_location is not None
                        assert protein.protein_name is not None and len(protein.protein_name) > 0
                        # Note: ProteinModel doesn't have description field
                        assert protein.organism_name == "Homo sapiens"
    
    @given(
        interpro_protein=interpro_protein_strategy(),
        isoform_data=complete_isoform_data_strategy()
    )
    @settings(max_examples=100, deadline=5000)
    @pytest.mark.asyncio
    async def test_exon_count_calculation_accuracy(self, interpro_protein, isoform_data):
        """
        Property 4: Complete Isoform Data Collection - Exon count calculation.
        
        For any isoform with exon annotations, the system should accurately
        calculate the exon count from the exon annotation data.
        
        **Validates: Requirements 3.3**
        """
        assume(len(interpro_protein.uniprot_id) >= 6)
        assume(isoform_data["exon_annotations"])
        assume(isoform_data["sequence"])
        
        with patch('protein_data_collector.api.uniprot_client.UnifiedUniProtClient') as mock_client_class:
            # Setup mock client
            mock_client = AsyncMock()
            mock_client.get_protein_isoforms.return_value = [isoform_data]
            mock_client._determine_access_method.return_value = "rest"
            
            # Mock the parse method to return data as-is initially
            def mock_parse_isoform_data(isoform_data, parent_id, method):
                return ProteinModel(
                    isoform_id=isoform_data["isoform_id"],
                    parent_protein_id=parent_id,
            parent_tim_barrel_accession="PF00001",  # Default test value
                    sequence=isoform_data["sequence"],
                    sequence_length=isoform_data["length"],
                    exon_annotations=isoform_data["exon_annotations"],
                    exon_count=None,  # Let the collector calculate this
                    tim_barrel_location=isoform_data["tim_barrel_location"],
                    organism_name=isoform_data["organism"],
                    protein_name=isoform_data["name"],
                )
            
            # Set the mock method directly
            mock_client.parse_protein_isoform_data = mock_parse_isoform_data
            mock_client_class.return_value = mock_client
            
            # Test the collector
            async with UniProtIsoformCollector(mock_client) as collector:
                collector.reset_collection_state()
                result_proteins = await collector.collect_protein_isoforms(interpro_protein)
                
                assert len(result_proteins) == 1
                protein = result_proteins[0]
                
                # Verify exon count was calculated correctly
                expected_exon_count = len(isoform_data["exon_annotations"]["exons"])
                assert protein.exon_count == expected_exon_count
                
                # Verify exon annotations are preserved
                assert protein.exon_annotations == isoform_data["exon_annotations"]
    
    @given(
        interpro_protein=interpro_protein_strategy(),
        sequence_length=st.integers(min_value=100, max_value=1000)
    )
    @settings(max_examples=100, deadline=5000)
    @pytest.mark.asyncio
    async def test_tim_barrel_location_validation(self, interpro_protein, sequence_length):
        """
        Property 4: Complete Isoform Data Collection - TIM barrel location validation.
        
        For any isoform with TIM barrel location data, the coordinates should be
        within the bounds of the protein sequence length.
        
        **Validates: Requirements 3.4**
        """
        assume(len(interpro_protein.uniprot_id) >= 6)
        assume(sequence_length >= 100)
        
        # Generate valid TIM barrel location within sequence bounds
        tim_start = sequence_length // 4
        tim_end = (sequence_length * 3) // 4
        
        isoform_data = {
            "isoform_id": f"{interpro_protein.uniprot_id}-1",
            "sequence": "A" * sequence_length,
            "length": sequence_length,
            "exon_annotations": {"exons": [{"start": 1, "end": sequence_length}]},
            "exon_count": 1,
            "tim_barrel_location": {
                "start": tim_start,
                "end": tim_end,
                "confidence": 0.95
            },
            "organism": "Homo sapiens",
            "name": "Test protein",
            "description": "Test description"
        }
        
        with patch('protein_data_collector.api.uniprot_client.UnifiedUniProtClient') as mock_client_class:
            # Setup mock client
            mock_client = AsyncMock()
            mock_client.get_protein_isoforms.return_value = [isoform_data]
            mock_client._determine_access_method.return_value = "rest"
            
            # Mock the parse method
            def mock_parse_isoform_data(isoform_data, parent_id, method):
                return ProteinModel(
                    isoform_id=isoform_data["isoform_id"],
                    parent_protein_id=parent_id,
            parent_tim_barrel_accession="PF00001",  # Default test value
                    sequence=isoform_data["sequence"],
                    sequence_length=isoform_data["length"],
                    exon_annotations=isoform_data["exon_annotations"],
                    exon_count=isoform_data["exon_count"],
                    tim_barrel_location=isoform_data["tim_barrel_location"],
                    organism_name=isoform_data["organism"],
                    protein_name=isoform_data["name"],
                )
            
            # Set the mock method directly
            mock_client.parse_protein_isoform_data = mock_parse_isoform_data
            mock_client_class.return_value = mock_client
            
            # Test the collector
            async with UniProtIsoformCollector(mock_client) as collector:
                collector.reset_collection_state()
                result_proteins = await collector.collect_protein_isoforms(interpro_protein)
                
                assert len(result_proteins) == 1
                protein = result_proteins[0]
                
                # Verify TIM barrel location is within sequence bounds
                tim_location = protein.tim_barrel_location
                assert tim_location is not None
                assert isinstance(tim_location, dict)
                
                start = tim_location.get("start", 0)
                end = tim_location.get("end", 0)
                
                assert start > 0
                assert end > start
                assert start <= protein.sequence_length
                assert end <= protein.sequence_length
    
    @given(interpro_proteins=st.lists(interpro_protein_strategy(), min_size=1, max_size=5, unique_by=lambda x: x.uniprot_id))
    @settings(max_examples=100, deadline=5000)
    @pytest.mark.asyncio
    async def test_collection_report_accuracy(self, interpro_proteins):
        """
        Property 4: Complete Isoform Data Collection - Collection reporting accuracy.
        
        For any collection operation, the system should provide accurate reporting
        of collection statistics including success rates and isoform counts.
        
        **Validates: Requirements 3.2, 3.3, 3.4, 3.5, 3.6**
        """
        assume(len(interpro_proteins) >= 1)
        assume(all(len(p.uniprot_id) >= 6 for p in interpro_proteins))
        
        # Generate isoform data - some proteins succeed, some fail
        successful_proteins = interpro_proteins[:len(interpro_proteins)//2 + 1]
        failed_proteins = interpro_proteins[len(interpro_proteins)//2 + 1:]
        
        isoforms_per_successful = 2
        total_expected_isoforms = len(successful_proteins) * isoforms_per_successful
        
        with patch('protein_data_collector.api.uniprot_client.UnifiedUniProtClient') as mock_client_class:
            # Setup mock client
            mock_client = AsyncMock()
            
            def mock_get_isoforms(uniprot_id):
                if any(p.uniprot_id == uniprot_id for p in successful_proteins):
                    return [
                        {
                            "isoform_id": f"{uniprot_id}-{i+1}",
                            "sequence": "ACDEFGHIKLMNPQRSTVWY" * 10,
                            "length": 200,
                            "exon_annotations": {"exons": [{"start": 1, "end": 200}]},
                            "exon_count": 1,
                            "tim_barrel_location": {"start": 10, "end": 190},
                            "organism": "Homo sapiens",
                            "name": f"Protein {uniprot_id}",
                            "description": f"Description {uniprot_id}"
                        }
                        for i in range(isoforms_per_successful)
                    ]
                else:
                    # Simulate failure for some proteins
                    raise Exception(f"Failed to get isoforms for {uniprot_id}")
            
            mock_client.get_protein_isoforms.side_effect = mock_get_isoforms
            mock_client._determine_access_method.return_value = "rest"
            
            # Mock the parse method
            def mock_parse_isoform_data(isoform_data, parent_id, method):
                return ProteinModel(
                    isoform_id=isoform_data["isoform_id"],
                    parent_protein_id=parent_id,
            parent_tim_barrel_accession="PF00001",  # Default test value
                    sequence=isoform_data["sequence"],
                    sequence_length=isoform_data["length"],
                    exon_annotations=isoform_data["exon_annotations"],
                    exon_count=isoform_data["exon_count"],
                    tim_barrel_location=isoform_data["tim_barrel_location"],
                    organism_name=isoform_data["organism"],
                    protein_name=isoform_data["name"],
                )
            
            # Set the mock method directly
            mock_client.parse_protein_isoform_data = mock_parse_isoform_data
            mock_client_class.return_value = mock_client
            
            # Test collection with reporting
            async with UniProtIsoformCollector(mock_client) as collector:
                collector.reset_collection_state()
                proteins, report = await collector.collect_isoforms_with_report(interpro_proteins)
                
                # Verify report accuracy
                assert report.total_proteins_processed == len(interpro_proteins)
                assert report.successful_proteins == len(successful_proteins)
                assert report.failed_proteins == len(failed_proteins)
                assert report.total_isoforms_collected == total_expected_isoforms
                
                # Verify success rate calculation
                expected_success_rate = (len(successful_proteins) / len(interpro_proteins)) * 100
                assert abs(report.success_rate - expected_success_rate) < 0.01
                
                # Verify average isoforms per protein
                if report.successful_proteins > 0:
                    expected_avg = report.total_isoforms_collected / report.successful_proteins
                    assert abs(report.average_isoforms_per_protein - expected_avg) < 0.01
                
                # Verify timing information
                assert report.start_time is not None
                assert report.end_time is not None
                assert report.collection_duration is not None
                assert report.collection_duration > 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])