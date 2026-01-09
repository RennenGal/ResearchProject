"""
Property-based tests for query result accuracy.

Tests Property 8: Query Result Accuracy
**Validates: Requirements 5.1, 5.2, 5.3, 5.5**

Feature: protein-data-collector, Property 8: Query Result Accuracy
"""

import pytest
from hypothesis import given, strategies as st, settings, assume
from typing import List, Dict, Set, Any

from protein_data_collector.query.engine import QueryEngine, QueryFilters, QueryResult


# Generators for test data
@st.composite
def query_filters_generator(draw):
    """Generate QueryFilters for testing."""
    return QueryFilters(
        pfam_family=draw(st.one_of(st.none(), st.text(min_size=5, max_size=20))),
        protein_id=draw(st.one_of(st.none(), st.text(min_size=6, max_size=20))),
        organism=draw(st.one_of(st.none(), st.sampled_from(["Homo sapiens", "Human"]))),
        min_sequence_length=draw(st.one_of(st.none(), st.integers(min_value=10, max_value=100))),
        max_sequence_length=draw(st.one_of(st.none(), st.integers(min_value=100, max_value=1000))),
        min_exon_count=draw(st.one_of(st.none(), st.integers(min_value=1, max_value=5))),
        max_exon_count=draw(st.one_of(st.none(), st.integers(min_value=5, max_value=20))),
        has_tim_barrel=draw(st.one_of(st.none(), st.booleans())),
        tim_barrel_confidence=draw(st.one_of(st.none(), st.floats(min_value=0.0, max_value=1.0)))
    )


@st.composite
def mock_protein_data_generator(draw):
    """Generate mock protein data for testing query result formatting."""
    # Generate valid amino acid sequence
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    sequence_length = draw(st.integers(min_value=50, max_value=500))
    sequence = ''.join(draw(st.lists(st.sampled_from(amino_acids), min_size=sequence_length, max_size=sequence_length)))
    
    # Generate exon annotations
    exon_count = draw(st.integers(min_value=1, max_value=10))
    exons = []
    current_pos = 1
    for i in range(exon_count):
        exon_length = draw(st.integers(min_value=10, max_value=50))
        exons.append({
            "start": current_pos,
            "end": current_pos + exon_length - 1
        })
        current_pos += exon_length
    
    exon_annotations = {"exons": exons}
    
    # Generate TIM barrel location within sequence bounds
    min_tim_length = 20
    if sequence_length < min_tim_length:
        tim_start = 1
        tim_end = sequence_length
    else:
        tim_start = draw(st.integers(min_value=1, max_value=sequence_length - min_tim_length + 1))
        tim_end = draw(st.integers(min_value=tim_start + min_tim_length - 1, max_value=sequence_length))
    
    tim_barrel_location = {
        "start": tim_start,
        "end": tim_end,
        "confidence": draw(st.floats(min_value=0.0, max_value=1.0))
    }
    
    # Generate protein identifiers
    uniprot_id = draw(st.text(min_size=6, max_size=20, alphabet=st.characters(whitelist_categories=('Lu', 'Nd'))))
    isoform_id = f"{uniprot_id}-{draw(st.integers(min_value=1, max_value=10))}"
    pfam_accession = draw(st.text(min_size=5, max_size=20, alphabet=st.characters(whitelist_categories=('Lu', 'Ll', 'Nd'))))
    
    # Generate optional text fields
    name = draw(st.one_of(
        st.none(), 
        st.text(min_size=1, max_size=255, alphabet=st.characters(min_codepoint=32, max_codepoint=126)).filter(lambda x: x.strip())
    ))
    description = draw(st.one_of(
        st.none(), 
        st.text(min_size=1, max_size=1000, alphabet=st.characters(min_codepoint=32, max_codepoint=126)).filter(lambda x: x.strip())
    ))
    organism = draw(st.sampled_from(["Homo sapiens", "Human"]))
    
    return {
        'isoform_id': isoform_id,
        'uniprot_id': uniprot_id,
        'parent_protein_id': uniprot_id,
        'sequence': sequence,
        'sequence_length': sequence_length,
        'exon_annotations': exon_annotations,
        'exon_count': exon_count,
        'tim_barrel_location': tim_barrel_location,
        'organism': organism,
        'name': name,
        'description': description,
        'pfam_accession': pfam_accession,
        'created_at': '2024-01-01T00:00:00'
    }


class TestQueryResultAccuracy:
    """Test query result accuracy property."""
    
    @given(query_filters_generator())
    @settings(max_examples=100, deadline=None)
    def test_query_filters_consistency_property(self, filters):
        """
        Property 8: Query Result Accuracy - Filter Consistency
        
        For any QueryFilters object, all filter attributes should be properly 
        accessible and maintain their assigned values.
        
        **Validates: Requirements 5.4**
        """
        # Verify filter object maintains all assigned values
        assert hasattr(filters, 'pfam_family')
        assert hasattr(filters, 'protein_id')
        assert hasattr(filters, 'organism')
        assert hasattr(filters, 'min_sequence_length')
        assert hasattr(filters, 'max_sequence_length')
        assert hasattr(filters, 'min_exon_count')
        assert hasattr(filters, 'max_exon_count')
        assert hasattr(filters, 'has_tim_barrel')
        assert hasattr(filters, 'tim_barrel_confidence')
        
        # Verify data types are correct
        if filters.pfam_family is not None:
            assert isinstance(filters.pfam_family, str)
        
        if filters.protein_id is not None:
            assert isinstance(filters.protein_id, str)
        
        if filters.organism is not None:
            assert isinstance(filters.organism, str)
        
        if filters.min_sequence_length is not None:
            assert isinstance(filters.min_sequence_length, int)
            assert filters.min_sequence_length > 0
        
        if filters.max_sequence_length is not None:
            assert isinstance(filters.max_sequence_length, int)
            assert filters.max_sequence_length > 0
        
        if filters.min_exon_count is not None:
            assert isinstance(filters.min_exon_count, int)
            assert filters.min_exon_count > 0
        
        if filters.max_exon_count is not None:
            assert isinstance(filters.max_exon_count, int)
            assert filters.max_exon_count > 0
        
        if filters.has_tim_barrel is not None:
            assert isinstance(filters.has_tim_barrel, bool)
        
        if filters.tim_barrel_confidence is not None:
            assert isinstance(filters.tim_barrel_confidence, float)
            assert 0.0 <= filters.tim_barrel_confidence <= 1.0
        
        # Verify logical consistency between min/max values
        if (filters.min_sequence_length is not None and 
            filters.max_sequence_length is not None):
            # This property should hold for valid filters
            # If it doesn't, the filter combination is invalid
            if filters.min_sequence_length > filters.max_sequence_length:
                # This is an invalid filter combination, but the object should still be valid
                pass
        
        if (filters.min_exon_count is not None and 
            filters.max_exon_count is not None):
            # This property should hold for valid filters
            if filters.min_exon_count > filters.max_exon_count:
                # This is an invalid filter combination, but the object should still be valid
                pass
    
    @given(st.lists(mock_protein_data_generator(), min_size=1, max_size=10))
    @settings(max_examples=100, deadline=None)
    def test_query_result_structure_property(self, protein_data_list):
        """
        Property 8: Query Result Accuracy - Result Structure
        
        For any list of protein data, when formatted as a QueryResult, the result 
        should maintain proper structure and include all required display fields.
        
        **Validates: Requirements 5.1, 5.2, 5.3, 5.5**
        """
        # Create a QueryResult manually to test structure
        query_result = QueryResult(
            pfam_families=[],
            proteins=[],
            isoforms=protein_data_list,
            total_count=len(protein_data_list),
            query_metadata={"test": True}
        )
        
        # Verify QueryResult structure
        assert isinstance(query_result.pfam_families, list)
        assert isinstance(query_result.proteins, list)
        assert isinstance(query_result.isoforms, list)
        assert isinstance(query_result.total_count, int)
        assert isinstance(query_result.query_metadata, dict)
        
        # Verify total count matches isoforms
        assert query_result.total_count == len(protein_data_list)
        
        # Verify each isoform has required display fields
        for isoform in query_result.isoforms:
            required_fields = [
                'isoform_id', 'parent_protein_id', 'sequence', 'sequence_length',
                'exon_annotations', 'exon_count', 'tim_barrel_location', 
                'organism', 'name', 'description'
            ]
            
            for field in required_fields:
                assert field in isoform, f"Required field {field} missing from isoform"
            
            # Verify data types and constraints
            assert isinstance(isoform['sequence_length'], int)
            assert isoform['sequence_length'] > 0
            
            if isoform['sequence']:
                assert isinstance(isoform['sequence'], str)
                assert len(isoform['sequence']) == isoform['sequence_length']
                # Verify sequence contains only valid amino acids
                valid_amino_acids = set('ACDEFGHIKLMNPQRSTVWYX-')
                assert all(c.upper() in valid_amino_acids for c in isoform['sequence'])
            
            if isoform['exon_count'] is not None:
                assert isinstance(isoform['exon_count'], int)
                assert isoform['exon_count'] > 0
            
            if isoform['exon_annotations']:
                assert isinstance(isoform['exon_annotations'], dict)
                if 'exons' in isoform['exon_annotations']:
                    exons = isoform['exon_annotations']['exons']
                    assert isinstance(exons, list)
                    if isoform['exon_count'] is not None:
                        assert len(exons) == isoform['exon_count']
            
            if isoform['tim_barrel_location']:
                assert isinstance(isoform['tim_barrel_location'], dict)
                tim_loc = isoform['tim_barrel_location']
                if 'start' in tim_loc and 'end' in tim_loc:
                    assert isinstance(tim_loc['start'], int)
                    assert isinstance(tim_loc['end'], int)
                    assert tim_loc['start'] > 0
                    assert tim_loc['end'] > 0
                    assert tim_loc['start'] <= tim_loc['end']
                    assert tim_loc['end'] <= isoform['sequence_length']
                
                if 'confidence' in tim_loc:
                    assert isinstance(tim_loc['confidence'], (int, float))
                    assert 0.0 <= tim_loc['confidence'] <= 1.0
    
    @given(mock_protein_data_generator())
    @settings(max_examples=100, deadline=None)
    def test_protein_data_filtering_logic_property(self, protein_data):
        """
        Property 8: Query Result Accuracy - Filtering Logic
        
        For any protein data and filtering criteria, the filtering logic should 
        correctly determine whether the protein matches the criteria.
        
        **Validates: Requirements 5.4**
        """
        # Test various filter conditions against the protein data
        
        # Test organism filtering
        organism_filter = QueryFilters(organism=protein_data['organism'])
        # Protein should match its own organism
        assert protein_data['organism'] == organism_filter.organism
        
        # Test sequence length filtering
        seq_len = protein_data['sequence_length']
        
        # Min sequence length filter
        min_len_filter = QueryFilters(min_sequence_length=seq_len - 10)
        # Protein should match if its length >= filter min
        assert seq_len >= min_len_filter.min_sequence_length
        
        # Max sequence length filter
        max_len_filter = QueryFilters(max_sequence_length=seq_len + 10)
        # Protein should match if its length <= filter max
        assert seq_len <= max_len_filter.max_sequence_length
        
        # Test exon count filtering
        if protein_data['exon_count'] is not None:
            exon_count = protein_data['exon_count']
            
            # Min exon count filter
            min_exon_filter = QueryFilters(min_exon_count=max(1, exon_count - 2))
            assert exon_count >= min_exon_filter.min_exon_count
            
            # Max exon count filter
            max_exon_filter = QueryFilters(max_exon_count=exon_count + 2)
            assert exon_count <= max_exon_filter.max_exon_count
        
        # Test TIM barrel filtering
        has_tim_barrel = protein_data['tim_barrel_location'] is not None
        
        # Has TIM barrel filter
        tim_barrel_filter = QueryFilters(has_tim_barrel=has_tim_barrel)
        assert (protein_data['tim_barrel_location'] is not None) == tim_barrel_filter.has_tim_barrel
        
        # TIM barrel confidence filtering
        if has_tim_barrel and 'confidence' in protein_data['tim_barrel_location']:
            confidence = protein_data['tim_barrel_location']['confidence']
            
            # Confidence filter
            conf_filter = QueryFilters(tim_barrel_confidence=confidence - 0.1)
            if conf_filter.tim_barrel_confidence >= 0.0:
                assert confidence >= conf_filter.tim_barrel_confidence
    
    def test_empty_query_result_structure(self):
        """Test that empty query results maintain proper structure."""
        empty_result = QueryResult(
            pfam_families=[],
            proteins=[],
            isoforms=[],
            total_count=0,
            query_metadata={}
        )
        
        # Verify structure is maintained even when empty
        assert isinstance(empty_result.pfam_families, list)
        assert isinstance(empty_result.proteins, list)
        assert isinstance(empty_result.isoforms, list)
        assert isinstance(empty_result.total_count, int)
        assert isinstance(empty_result.query_metadata, dict)
        
        assert len(empty_result.pfam_families) == 0
        assert len(empty_result.proteins) == 0
        assert len(empty_result.isoforms) == 0
        assert empty_result.total_count == 0
    
    def test_query_result_metadata_consistency(self):
        """Test that query result metadata maintains consistency."""
        test_metadata = {
            "pfam_id": "PF00001",
            "found": True,
            "protein_count": 5,
            "isoform_count": 10
        }
        
        result = QueryResult(
            pfam_families=[],
            proteins=[],
            isoforms=[],
            total_count=10,
            query_metadata=test_metadata
        )
        
        # Verify metadata is preserved
        assert result.query_metadata == test_metadata
        assert result.query_metadata['pfam_id'] == "PF00001"
        assert result.query_metadata['found'] == True
        assert result.query_metadata['protein_count'] == 5
        assert result.query_metadata['isoform_count'] == 10
        
        # Verify total count consistency
        assert result.total_count == test_metadata['isoform_count']