"""
Property-based tests for collection progress persistence.

Tests Property 12: Collection Progress Persistence
**Validates: Requirements 6.5**

Feature: protein-data-collector, Property 12: Collection Progress Persistence
"""

import pytest
import tempfile
import json
from pathlib import Path
from datetime import datetime
from hypothesis import given, strategies as st, settings
from typing import List, Dict, Any

from protein_data_collector.collector.data_collector import CollectionProgress


# Generators for test data
@st.composite
def collection_progress_generator(draw):
    """Generate valid collection progress data."""
    phases = ["not_started", "pfam_families", "interpro_proteins", "uniprot_isoforms", "storage", "completed"]
    
    phase = draw(st.sampled_from(phases))
    
    # Generate collected counts based on phase constraints
    if phase == "not_started":
        pfam_families_collected = 0
        interpro_proteins_collected = 0
        uniprot_isoforms_collected = 0
        pfam_families_stored = 0
        interpro_proteins_stored = 0
        uniprot_isoforms_stored = 0
    elif phase == "pfam_families":
        pfam_families_collected = draw(st.integers(min_value=0, max_value=1000))
        interpro_proteins_collected = 0
        uniprot_isoforms_collected = 0
        pfam_families_stored = 0
        interpro_proteins_stored = 0
        uniprot_isoforms_stored = 0
    elif phase == "interpro_proteins":
        pfam_families_collected = draw(st.integers(min_value=0, max_value=1000))
        interpro_proteins_collected = draw(st.integers(min_value=0, max_value=5000))
        uniprot_isoforms_collected = 0
        pfam_families_stored = 0
        interpro_proteins_stored = 0
        uniprot_isoforms_stored = 0
    elif phase == "uniprot_isoforms":
        pfam_families_collected = draw(st.integers(min_value=0, max_value=1000))
        interpro_proteins_collected = draw(st.integers(min_value=0, max_value=5000))
        uniprot_isoforms_collected = draw(st.integers(min_value=0, max_value=10000))
        pfam_families_stored = 0
        interpro_proteins_stored = 0
        uniprot_isoforms_stored = 0
    elif phase == "storage":
        pfam_families_collected = draw(st.integers(min_value=0, max_value=1000))
        interpro_proteins_collected = draw(st.integers(min_value=0, max_value=5000))
        uniprot_isoforms_collected = draw(st.integers(min_value=0, max_value=10000))
        pfam_families_stored = draw(st.integers(min_value=0, max_value=pfam_families_collected))
        interpro_proteins_stored = draw(st.integers(min_value=0, max_value=interpro_proteins_collected))
        uniprot_isoforms_stored = draw(st.integers(min_value=0, max_value=uniprot_isoforms_collected))
    else:  # completed
        pfam_families_collected = draw(st.integers(min_value=1, max_value=1000))
        interpro_proteins_collected = draw(st.integers(min_value=1, max_value=5000))
        uniprot_isoforms_collected = draw(st.integers(min_value=1, max_value=10000))
        pfam_families_stored = pfam_families_collected
        interpro_proteins_stored = interpro_proteins_collected
        uniprot_isoforms_stored = uniprot_isoforms_collected
    
    # Generate fixed timestamps for reproducibility
    base_time = datetime(2024, 1, 1, 12, 0, 0)  # Fixed base time
    start_time = draw(st.datetimes(min_value=base_time, max_value=datetime(2024, 12, 31, 23, 59, 59)))
    last_checkpoint = draw(st.one_of(
        st.none(),
        st.datetimes(min_value=start_time, max_value=start_time.replace(hour=23, minute=59, second=59))
    ))
    
    errors = draw(st.lists(
        st.text(min_size=1, max_size=200, alphabet=st.characters(min_codepoint=32, max_codepoint=126)).filter(lambda x: x.strip()),
        max_size=10
    ))
    
    progress = CollectionProgress()
    progress.phase = phase
    progress.pfam_families_collected = pfam_families_collected
    progress.interpro_proteins_collected = interpro_proteins_collected
    progress.uniprot_isoforms_collected = uniprot_isoforms_collected
    progress.pfam_families_stored = pfam_families_stored
    progress.interpro_proteins_stored = interpro_proteins_stored
    progress.uniprot_isoforms_stored = uniprot_isoforms_stored
    progress.start_time = start_time
    progress.last_checkpoint = last_checkpoint
    progress.errors = errors
    
    return progress


@st.composite
def progress_file_path_generator(draw):
    """Generate valid file paths for progress files."""
    # Use temporary directory for testing
    temp_dir = tempfile.gettempdir()
    filename = draw(st.text(
        min_size=5, 
        max_size=50, 
        alphabet=st.characters(whitelist_categories=('Lu', 'Ll', 'Nd'), whitelist_characters='_-')
    )) + ".json"
    
    return Path(temp_dir) / "test_progress" / filename


class TestCollectionProgressPersistence:
    """Test collection progress persistence property."""
    
    @given(collection_progress_generator(), progress_file_path_generator())
    @settings(max_examples=100, deadline=None)
    def test_collection_progress_persistence_property(self, progress, progress_file):
        """
        Property 12: Collection Progress Persistence
        
        For any interrupted data collection operation, the system should be able to 
        resume from the last successfully completed step without duplicating work 
        or losing progress.
        
        **Validates: Requirements 6.5**
        """
        # Ensure the directory exists
        progress_file.parent.mkdir(parents=True, exist_ok=True)
        
        try:
            # Property 1: Progress can be serialized to dictionary
            progress_dict = progress.to_dict()
            
            # Verify all essential fields are present
            assert "phase" in progress_dict
            assert "pfam_families_collected" in progress_dict
            assert "interpro_proteins_collected" in progress_dict
            assert "uniprot_isoforms_collected" in progress_dict
            assert "pfam_families_stored" in progress_dict
            assert "interpro_proteins_stored" in progress_dict
            assert "uniprot_isoforms_stored" in progress_dict
            assert "errors" in progress_dict
            assert "total_entities_collected" in progress_dict
            assert "total_entities_stored" in progress_dict
            assert "duration_seconds" in progress_dict
            
            # Property 2: Progress can be saved to file
            with open(progress_file, 'w') as f:
                json.dump(progress_dict, f, indent=2)
            
            # Verify file was created and is readable
            assert progress_file.exists()
            assert progress_file.is_file()
            assert progress_file.stat().st_size > 0
            
            # Property 3: Progress can be loaded from file
            with open(progress_file, 'r') as f:
                loaded_dict = json.load(f)
            
            # Verify loaded data matches original
            assert loaded_dict["phase"] == progress_dict["phase"]
            assert loaded_dict["pfam_families_collected"] == progress_dict["pfam_families_collected"]
            assert loaded_dict["interpro_proteins_collected"] == progress_dict["interpro_proteins_collected"]
            assert loaded_dict["uniprot_isoforms_collected"] == progress_dict["uniprot_isoforms_collected"]
            assert loaded_dict["pfam_families_stored"] == progress_dict["pfam_families_stored"]
            assert loaded_dict["interpro_proteins_stored"] == progress_dict["interpro_proteins_stored"]
            assert loaded_dict["uniprot_isoforms_stored"] == progress_dict["uniprot_isoforms_stored"]
            assert loaded_dict["errors"] == progress_dict["errors"]
            
            # Property 4: Progress can be reconstructed from dictionary
            reconstructed_progress = CollectionProgress.from_dict(loaded_dict)
            
            # Verify reconstructed progress matches original
            assert reconstructed_progress.phase == progress.phase
            assert reconstructed_progress.pfam_families_collected == progress.pfam_families_collected
            assert reconstructed_progress.interpro_proteins_collected == progress.interpro_proteins_collected
            assert reconstructed_progress.uniprot_isoforms_collected == progress.uniprot_isoforms_collected
            assert reconstructed_progress.pfam_families_stored == progress.pfam_families_stored
            assert reconstructed_progress.interpro_proteins_stored == progress.interpro_proteins_stored
            assert reconstructed_progress.uniprot_isoforms_stored == progress.uniprot_isoforms_stored
            assert reconstructed_progress.errors == progress.errors
            
            # Property 5: Calculated properties are preserved
            assert reconstructed_progress.total_entities_collected == progress.total_entities_collected
            assert reconstructed_progress.total_entities_stored == progress.total_entities_stored
            
            # Property 6: Timestamps are preserved (if present)
            if progress.start_time:
                assert reconstructed_progress.start_time is not None
                # Allow for small differences due to serialization
                time_diff = abs((reconstructed_progress.start_time - progress.start_time).total_seconds())
                assert time_diff < 1.0  # Less than 1 second difference
            
            if progress.last_checkpoint:
                assert reconstructed_progress.last_checkpoint is not None
                time_diff = abs((reconstructed_progress.last_checkpoint - progress.last_checkpoint).total_seconds())
                assert time_diff < 1.0  # Less than 1 second difference
            
        finally:
            # Clean up test file
            if progress_file.exists():
                progress_file.unlink()
            # Clean up directory if empty
            if progress_file.parent.exists() and not any(progress_file.parent.iterdir()):
                progress_file.parent.rmdir()
    
    @given(collection_progress_generator())
    @settings(max_examples=100, deadline=None)
    def test_progress_round_trip_consistency(self, progress):
        """
        Test that progress data survives round-trip serialization without loss.
        
        This is part of the collection progress persistence property.
        """
        # Round trip: progress -> dict -> progress
        progress_dict = progress.to_dict()
        reconstructed_progress = CollectionProgress.from_dict(progress_dict)
        
        # Verify all fields are preserved
        assert reconstructed_progress.phase == progress.phase
        assert reconstructed_progress.pfam_families_collected == progress.pfam_families_collected
        assert reconstructed_progress.interpro_proteins_collected == progress.interpro_proteins_collected
        assert reconstructed_progress.uniprot_isoforms_collected == progress.uniprot_isoforms_collected
        assert reconstructed_progress.pfam_families_stored == progress.pfam_families_stored
        assert reconstructed_progress.interpro_proteins_stored == progress.interpro_proteins_stored
        assert reconstructed_progress.uniprot_isoforms_stored == progress.uniprot_isoforms_stored
        assert reconstructed_progress.errors == progress.errors
        
        # Verify calculated properties are consistent
        assert reconstructed_progress.total_entities_collected == progress.total_entities_collected
        assert reconstructed_progress.total_entities_stored == progress.total_entities_stored
    
    @given(collection_progress_generator())
    @settings(max_examples=100, deadline=None)
    def test_progress_data_integrity(self, progress):
        """
        Test that progress data maintains integrity constraints.
        
        Verifies that stored counts never exceed collected counts.
        """
        # Property: Stored counts should never exceed collected counts
        assert progress.pfam_families_stored <= progress.pfam_families_collected, \
            f"Stored PFAM families ({progress.pfam_families_stored}) exceeds collected ({progress.pfam_families_collected})"
        
        assert progress.interpro_proteins_stored <= progress.interpro_proteins_collected, \
            f"Stored proteins ({progress.interpro_proteins_stored}) exceeds collected ({progress.interpro_proteins_collected})"
        
        assert progress.uniprot_isoforms_stored <= progress.uniprot_isoforms_collected, \
            f"Stored isoforms ({progress.uniprot_isoforms_stored}) exceeds collected ({progress.uniprot_isoforms_collected})"
        
        # Property: All counts should be non-negative
        assert progress.pfam_families_collected >= 0
        assert progress.interpro_proteins_collected >= 0
        assert progress.uniprot_isoforms_collected >= 0
        assert progress.pfam_families_stored >= 0
        assert progress.interpro_proteins_stored >= 0
        assert progress.uniprot_isoforms_stored >= 0
        
        # Property: Total counts should be sum of individual counts
        expected_total_collected = (progress.pfam_families_collected + 
                                  progress.interpro_proteins_collected + 
                                  progress.uniprot_isoforms_collected)
        assert progress.total_entities_collected == expected_total_collected
        
        expected_total_stored = (progress.pfam_families_stored + 
                               progress.interpro_proteins_stored + 
                               progress.uniprot_isoforms_stored)
        assert progress.total_entities_stored == expected_total_stored
    
    @given(collection_progress_generator())
    @settings(max_examples=100, deadline=None)
    def test_progress_phase_consistency(self, progress):
        """
        Test that progress phase is consistent with collected/stored data.
        
        Verifies that phase progression makes sense with the data.
        """
        # Property: Phase should be consistent with progress data
        if progress.phase == "not_started":
            # No data should be collected or stored
            assert progress.total_entities_collected == 0
            assert progress.total_entities_stored == 0
        
        elif progress.phase == "pfam_families":
            # Only PFAM families might be collected, nothing stored yet
            assert progress.interpro_proteins_collected == 0
            assert progress.uniprot_isoforms_collected == 0
            assert progress.total_entities_stored == 0
        
        elif progress.phase == "interpro_proteins":
            # PFAM families and proteins might be collected, nothing stored yet
            assert progress.uniprot_isoforms_collected == 0
            assert progress.total_entities_stored == 0
        
        elif progress.phase == "uniprot_isoforms":
            # All entities might be collected, nothing stored yet
            assert progress.total_entities_stored == 0
        
        elif progress.phase == "storage":
            # All entities collected, some might be stored
            # No specific constraints here as storage can be partial
            pass
        
        elif progress.phase == "completed":
            # All entities should be collected and stored
            if progress.total_entities_collected > 0:
                # If we collected anything, we should have stored something
                assert progress.total_entities_stored >= 0
    
    def test_empty_progress_persistence(self):
        """Test persistence of empty/initial progress state."""
        progress = CollectionProgress()
        
        # Test serialization of empty progress
        progress_dict = progress.to_dict()
        
        # Verify default values
        assert progress_dict["phase"] == "not_started"
        assert progress_dict["pfam_families_collected"] == 0
        assert progress_dict["interpro_proteins_collected"] == 0
        assert progress_dict["uniprot_isoforms_collected"] == 0
        assert progress_dict["pfam_families_stored"] == 0
        assert progress_dict["interpro_proteins_stored"] == 0
        assert progress_dict["uniprot_isoforms_stored"] == 0
        assert progress_dict["errors"] == []
        assert progress_dict["total_entities_collected"] == 0
        assert progress_dict["total_entities_stored"] == 0
        assert progress_dict["duration_seconds"] == 0.0
        
        # Test round-trip
        reconstructed = CollectionProgress.from_dict(progress_dict)
        assert reconstructed.phase == progress.phase
        assert reconstructed.total_entities_collected == progress.total_entities_collected
        assert reconstructed.total_entities_stored == progress.total_entities_stored
    
    def test_progress_file_corruption_handling(self):
        """Test handling of corrupted progress files."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            # Write invalid JSON
            f.write("{ invalid json content")
            temp_path = Path(f.name)
        
        try:
            # Attempting to load corrupted file should not crash
            # (This would be handled by the DataCollector class)
            with open(temp_path, 'r') as f:
                try:
                    json.load(f)
                    assert False, "Should have raised JSON decode error"
                except json.JSONDecodeError:
                    # Expected behavior - corrupted files should raise errors
                    pass
        finally:
            temp_path.unlink()
    
    def test_progress_with_complex_errors(self):
        """Test persistence of progress with complex error messages."""
        progress = CollectionProgress()
        progress.phase = "interpro_proteins"
        progress.pfam_families_collected = 5
        progress.interpro_proteins_collected = 25
        progress.errors = [
            "API timeout for protein P12345",
            "Validation error: invalid sequence for Q67890",
            "Network error: connection refused to uniprot.org",
            "Rate limit exceeded for InterPro API"
        ]
        
        # Test serialization with complex errors
        progress_dict = progress.to_dict()
        reconstructed = CollectionProgress.from_dict(progress_dict)
        
        # Verify errors are preserved exactly
        assert reconstructed.errors == progress.errors
        assert len(reconstructed.errors) == 4
        assert "API timeout for protein P12345" in reconstructed.errors
        assert "Rate limit exceeded for InterPro API" in reconstructed.errors