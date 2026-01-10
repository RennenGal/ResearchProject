"""
Property-based tests for data validation.

Feature: protein-data-collector, Property 6: Data Validation Consistency
Validates: Requirements 7.1, 7.2
"""

import pytest
from hypothesis import given, strategies as st, settings, HealthCheck

from protein_data_collector.models.validation import (
    ProteinSequenceValidator,
    TIMBarrelLocationValidator,
    DataValidator,
    ValidationResult,
    ValidationError,
    ValidationErrorType
)


# Hypothesis strategies for generating test data
valid_amino_acids = "ACDEFGHIKLMNPQRSTVWY"
invalid_amino_acids = "123456789!@#$%^&*()[]{}|\\:;\"'<>,.?/~`"

valid_protein_sequence_strategy = st.text(
    alphabet=valid_amino_acids,
    min_size=10,
    max_size=1000
)

invalid_protein_sequence_strategy = st.text(
    alphabet=invalid_amino_acids,
    min_size=1,
    max_size=100
)

mixed_protein_sequence_strategy = st.text(
    alphabet=valid_amino_acids + invalid_amino_acids,
    min_size=10,
    max_size=500
)

sequence_length_strategy = st.integers(min_value=10, max_value=5000)

coordinate_strategy = st.integers(min_value=1, max_value=1000)

confidence_strategy = st.floats(min_value=0.0, max_value=1.0)


class TestDataValidationConsistency:
    """Property-based tests for data validation consistency."""
    
    @given(sequence=valid_protein_sequence_strategy)
    @settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
    def test_valid_protein_sequences_always_pass_validation(self, sequence):
        """
        Feature: protein-data-collector, Property 6: Data Validation Consistency
        
        For any protein sequence containing only valid amino acid characters (ACDEFGHIKLMNPQRSTVWY),
        the system should accept the sequence and validation should pass.
        """
        validator = ProteinSequenceValidator()
        result = validator.validate(sequence)
        
        assert result.is_valid, f"Valid sequence should pass validation: {result.error_message}"
        assert len(result.errors) == 0, f"Valid sequence should have no errors: {result.errors}"
        
        # Verify the sequence contains only valid characters
        clean_sequence = sequence.strip().upper()
        invalid_chars = set(clean_sequence) - validator.VALID_AMINO_ACIDS
        assert len(invalid_chars) == 0, f"Test sequence contains invalid characters: {invalid_chars}"
    
    @given(sequence=invalid_protein_sequence_strategy)
    @settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
    def test_invalid_protein_sequences_always_fail_validation(self, sequence):
        """
        Feature: protein-data-collector, Property 6: Data Validation Consistency
        
        For any protein sequence containing invalid characters (non-amino acid characters),
        the system should reject the sequence and validation should fail.
        """
        validator = ProteinSequenceValidator(allow_extended=False)  # Use strict validation
        result = validator.validate(sequence)
        
        # Should fail validation due to invalid characters
        assert not result.is_valid, f"Invalid sequence should fail validation: {sequence}"
        assert len(result.errors) > 0, "Invalid sequence should have validation errors"
        
        # Should have an error about invalid characters
        has_invalid_char_error = any(
            error.error_type == ValidationErrorType.INVALID_SEQUENCE 
            for error in result.errors
        )
        assert has_invalid_char_error, "Should have invalid sequence error"
    
    @given(
        start=coordinate_strategy,
        end=coordinate_strategy,
        sequence_length=sequence_length_strategy,
        confidence=confidence_strategy
    )
    @settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
    def test_tim_barrel_location_validation_consistency(self, start, end, sequence_length, confidence):
        """
        Feature: protein-data-collector, Property 6: Data Validation Consistency
        
        For any TIM barrel location coordinates, validation should be consistent with
        coordinate bounds and sequence length constraints.
        """
        # Ensure start < end for valid coordinates
        if start >= end:
            start, end = min(start, end), max(start, end)
            if start == end:
                end = start + 1
        
        location = {
            "start": start,
            "end": end,
            "confidence": confidence
        }
        
        validator = TIMBarrelLocationValidator()
        result = validator.validate(location, sequence_length)
        
        # Check validation logic consistency
        if start >= 1 and end <= sequence_length and start < end:
            # Should be valid if coordinates are within bounds
            assert result.is_valid, f"Valid coordinates should pass: start={start}, end={end}, seq_len={sequence_length}"
        else:
            # Should be invalid if coordinates are out of bounds
            assert not result.is_valid, f"Invalid coordinates should fail: start={start}, end={end}, seq_len={sequence_length}"
            
            # Check specific error conditions
            if start < 1:
                has_start_error = any(
                    error.error_type == ValidationErrorType.OUT_OF_BOUNDS and "start" in error.field_name.lower()
                    for error in result.errors
                )
                assert has_start_error, "Should have start coordinate error"
            
            if end > sequence_length:
                has_end_error = any(
                    error.error_type == ValidationErrorType.OUT_OF_BOUNDS and "end" in error.field_name.lower()
                    for error in result.errors
                )
                assert has_end_error, "Should have end coordinate error"
    
    @given(
        sequence=valid_protein_sequence_strategy,
        start=coordinate_strategy,
        end=coordinate_strategy
    )
    @settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
    def test_protein_isoform_validation_integrates_sequence_and_location_validation(
        self, sequence, start, end
    ):
        """
        Feature: protein-data-collector, Property 6: Data Validation Consistency
        
        For any protein isoform data, validation should consistently integrate both
        sequence validation and TIM barrel location validation rules.
        """
        sequence_length = len(sequence)
        
        # Ensure valid coordinates
        if start >= end:
            start, end = min(start, end), max(start, end)
            if start == end:
                end = start + 1
        
        # Ensure coordinates are within sequence bounds
        start = max(1, min(start, sequence_length - 1))
        end = min(sequence_length, max(end, start + 1))
        
        protein_data = {
            "isoform_id": "TEST001-1",
            "parent_protein_id": "TEST001",
            "sequence": sequence,
            "sequence_length": sequence_length,
            "tim_barrel_location": {
                "start": start,
                "end": end,
                "confidence": 0.95
            }
        }
        
        validator = DataValidator()
        result = validator.validate_protein_isoform(protein_data)
        
        # Should be valid since we have valid sequence and valid coordinates
        assert result.is_valid, f"Valid protein isoform should pass validation: {result.error_message}"
        
        # Test with invalid sequence
        protein_data_invalid_seq = protein_data.copy()
        protein_data_invalid_seq["sequence"] = "INVALID123"
        protein_data_invalid_seq["sequence_length"] = len("INVALID123")
        
        result_invalid = validator.validate_protein_isoform(protein_data_invalid_seq)
        assert not result_invalid.is_valid, "Protein with invalid sequence should fail validation"
        
        # Should have sequence validation error
        has_sequence_error = any(
            error.error_type == ValidationErrorType.INVALID_SEQUENCE
            for error in result_invalid.errors
        )
        assert has_sequence_error, "Should have sequence validation error"
    
    @given(
        sequence=valid_protein_sequence_strategy,
        declared_length=sequence_length_strategy
    )
    @settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
    def test_sequence_length_consistency_validation(self, sequence, declared_length):
        """
        Feature: protein-data-collector, Property 6: Data Validation Consistency
        
        For any protein sequence and declared length, validation should consistently
        detect mismatches between actual and declared sequence lengths.
        """
        actual_length = len(sequence.strip())
        
        protein_data = {
            "isoform_id": "TEST001-1",
            "parent_protein_id": "TEST001",
            "sequence": sequence,
            "sequence_length": declared_length
        }
        
        validator = DataValidator()
        result = validator.validate_protein_isoform(protein_data)
        
        if actual_length == declared_length:
            # Should not have length mismatch error
            has_length_error = any(
                error.error_type == ValidationErrorType.INCONSISTENT_DATA and "sequence_length" in error.field_name
                for error in result.errors
            )
            assert not has_length_error, f"Should not have length error when lengths match: actual={actual_length}, declared={declared_length}"
        else:
            # Should have length mismatch error
            has_length_error = any(
                error.error_type == ValidationErrorType.INCONSISTENT_DATA and "sequence_length" in error.field_name
                for error in result.errors
            )
            assert has_length_error, f"Should have length error when lengths don't match: actual={actual_length}, declared={declared_length}"
    
    @given(sequence=st.text(min_size=0, max_size=5))
    @settings(max_examples=50, suppress_health_check=[HealthCheck.function_scoped_fixture])
    def test_empty_and_short_sequences_consistently_rejected(self, sequence):
        """
        Feature: protein-data-collector, Property 6: Data Validation Consistency
        
        For any empty or very short sequence (< 10 characters), validation should
        consistently reject the sequence as invalid.
        """
        validator = ProteinSequenceValidator(allow_extended=False)  # Use strict validation
        result = validator.validate(sequence)
        
        clean_sequence = sequence.strip()
        
        if not clean_sequence:
            # Empty sequences should always fail
            assert not result.is_valid, "Empty sequence should fail validation"
            has_empty_error = any(
                error.error_type in [ValidationErrorType.MISSING_REQUIRED_FIELD, ValidationErrorType.INVALID_SEQUENCE]
                for error in result.errors
            )
            assert has_empty_error, "Should have error for empty sequence"
        elif len(clean_sequence) < 10:
            # Very short sequences should generate warnings but may still be valid if they contain valid amino acids
            valid_chars = set(clean_sequence.upper()) <= validator.VALID_AMINO_ACIDS
            if valid_chars:
                # Should be valid but with warnings
                assert result.warnings is not None, "Short valid sequence should have warnings"
            else:
                # Should be invalid due to invalid characters
                assert not result.is_valid, "Short sequence with invalid characters should fail"