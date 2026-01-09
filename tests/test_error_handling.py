"""
Tests for comprehensive error handling system.

Tests error classification, logging, and recovery strategies
for different types of errors encountered during data collection.
"""

import pytest
import json
from unittest.mock import Mock, patch
from json import JSONDecodeError

from protein_data_collector.errors import (
    ErrorHandler, ErrorContext, ErrorCategory, ErrorSeverity, ErrorAction,
    NetworkError, APIError, DataError, ValidationError, DatabaseError, ConfigurationError,
    handle_error, create_error_context
)


class TestErrorClassification:
    """Test error classification and handling strategies."""
    
    def test_network_error_classification(self):
        """Test that network errors are properly classified."""
        handler = ErrorHandler()
        context = ErrorContext(operation="test_operation", database="TestDB")
        
        # Test ConnectionError
        error = ConnectionError("Connection failed")
        error_info = handler.classify_error(error, context)
        
        assert error_info.category == ErrorCategory.NETWORK
        assert error_info.severity == ErrorSeverity.MEDIUM
        assert error_info.action == ErrorAction.RETRY
        assert "network connectivity" in " ".join(error_info.recovery_suggestions).lower()
    
    def test_api_error_classification(self):
        """Test that API errors are properly classified."""
        handler = ErrorHandler()
        context = ErrorContext(operation="test_operation", database="TestDB")
        
        # Test generic HTTP error using a mock response
        class MockHTTPError(Exception):
            def __init__(self, message):
                super().__init__(message)
                self.response = Mock()
                self.response.status_code = 500
        
        MockHTTPError.__name__ = "HTTPError"
        
        error = MockHTTPError("HTTP error occurred")
        error_info = handler.classify_error(error, context)
        
        assert error_info.category == ErrorCategory.API
        assert error_info.severity == ErrorSeverity.MEDIUM
        assert error_info.action == ErrorAction.RETRY
        assert "api" in " ".join(error_info.recovery_suggestions).lower()
    
    def test_data_error_classification(self):
        """Test that data errors are properly classified."""
        handler = ErrorHandler()
        context = ErrorContext(operation="test_operation", database="TestDB")
        
        # Test JSON decode error
        error = JSONDecodeError("Invalid JSON", "test", 0)
        error_info = handler.classify_error(error, context)
        
        assert error_info.category == ErrorCategory.DATA
        assert error_info.severity == ErrorSeverity.HIGH
        assert error_info.action == ErrorAction.SKIP
        assert "data format" in " ".join(error_info.recovery_suggestions).lower()
    
    def test_custom_error_classification(self):
        """Test that custom errors maintain their classification."""
        handler = ErrorHandler()
        context = ErrorContext(operation="test_operation", database="TestDB")
        
        # Test custom NetworkError
        error = NetworkError("Custom network error", context=context)
        error_info = handler.classify_error(error, context)
        
        assert error_info.category == ErrorCategory.NETWORK
        assert error_info.severity == ErrorSeverity.MEDIUM
        assert error_info.action == ErrorAction.RETRY
        assert error_info.message == "Custom network error"
    
    def test_unknown_error_classification(self):
        """Test that unknown errors get default classification."""
        handler = ErrorHandler()
        context = ErrorContext(operation="test_operation", database="TestDB")
        
        # Test unknown error type
        error = RuntimeError("Unknown error")
        error_info = handler.classify_error(error, context)
        
        assert error_info.category == ErrorCategory.UNKNOWN
        assert error_info.severity == ErrorSeverity.MEDIUM
        assert error_info.action == ErrorAction.LOG_AND_CONTINUE


class TestErrorLogging:
    """Test error logging functionality."""
    
    @patch('protein_data_collector.errors.logging.getLogger')
    def test_error_logging_levels(self, mock_get_logger):
        """Test that errors are logged at appropriate levels."""
        mock_logger = Mock()
        mock_get_logger.return_value = mock_logger
        
        handler = ErrorHandler()
        context = ErrorContext(operation="test_operation", database="TestDB")
        
        # Test critical error logging
        critical_error = ConfigurationError("Critical config error", context=context)
        handler.handle_error(critical_error, context)
        mock_logger.critical.assert_called_once()
        
        # Reset mock
        mock_logger.reset_mock()
        
        # Test high severity error logging
        high_error = DataError("Data processing error", context=context)
        handler.handle_error(high_error, context)
        mock_logger.error.assert_called_once()
        
        # Reset mock
        mock_logger.reset_mock()
        
        # Test medium severity error logging
        medium_error = NetworkError("Network connection error", context=context)
        handler.handle_error(medium_error, context)
        mock_logger.warning.assert_called_once()
    
    @patch('protein_data_collector.errors.logging.getLogger')
    def test_contextual_logging(self, mock_get_logger):
        """Test that contextual information is included in logs."""
        mock_logger = Mock()
        mock_get_logger.return_value = mock_logger
        
        handler = ErrorHandler()
        context = ErrorContext(
            operation="test_operation",
            database="TestDB",
            entity_id="test_entity_123",
            entity_type="protein",
            request_url="https://api.example.com/test",
            response_status=500,
            additional_data={"custom_field": "custom_value"}
        )
        
        error = NetworkError("Test error", context=context)
        handler.handle_error(error, context)
        
        # Verify logger was called with contextual information
        mock_logger.warning.assert_called_once()
        call_args = mock_logger.warning.call_args
        
        # Check that extra data contains expected fields
        extra_data = call_args[1]['extra']
        assert extra_data['operation'] == "test_operation"
        assert extra_data['database'] == "TestDB"
        assert extra_data['entity_id'] == "test_entity_123"
        assert extra_data['entity_type'] == "protein"
        assert extra_data['request_url'] == "https://api.example.com/test"
        assert extra_data['response_status'] == 500
        assert extra_data['custom_field'] == "custom_value"


class TestErrorStatistics:
    """Test error statistics tracking."""
    
    def test_error_count_tracking(self):
        """Test that error counts are tracked correctly."""
        handler = ErrorHandler()
        context1 = ErrorContext(operation="operation1", database="DB1")
        context2 = ErrorContext(operation="operation2", database="DB2")
        
        # Generate some errors
        error1 = NetworkError("Error 1")
        error2 = APIError("Error 2")
        error3 = NetworkError("Error 3")
        
        handler.handle_error(error1, context1)
        handler.handle_error(error2, context2)
        handler.handle_error(error3, context1)  # Same operation as first error
        
        stats = handler.get_error_statistics()
        
        # Should have two different operation keys
        assert len(stats) == 2
        assert stats["network:operation1"] == 2  # Two network errors for operation1
        assert stats["api:operation2"] == 1     # One API error for operation2
    
    def test_error_statistics_reset(self):
        """Test that error statistics can be reset."""
        handler = ErrorHandler()
        context = ErrorContext(operation="test_operation", database="TestDB")
        
        # Generate an error
        error = NetworkError("Test error")
        handler.handle_error(error, context)
        
        # Verify statistics exist
        stats = handler.get_error_statistics()
        assert len(stats) > 0
        
        # Reset and verify empty
        handler.reset_error_statistics()
        stats = handler.get_error_statistics()
        assert len(stats) == 0


class TestErrorContext:
    """Test error context creation and usage."""
    
    def test_error_context_creation(self):
        """Test error context creation with various parameters."""
        context = create_error_context(
            operation="test_operation",
            database="TestDB",
            entity_id="entity_123",
            entity_type="protein",
            request_url="https://api.example.com/test",
            request_params={"param1": "value1"},
            response_status=404,
            custom_field="custom_value"
        )
        
        assert context.operation == "test_operation"
        assert context.database == "TestDB"
        assert context.entity_id == "entity_123"
        assert context.entity_type == "protein"
        assert context.request_url == "https://api.example.com/test"
        assert context.request_params == {"param1": "value1"}
        assert context.response_status == 404
        assert context.additional_data["custom_field"] == "custom_value"
        assert context.timestamp is not None
    
    def test_minimal_error_context(self):
        """Test error context creation with minimal parameters."""
        context = create_error_context(operation="minimal_operation")
        
        assert context.operation == "minimal_operation"
        assert context.database is None
        assert context.entity_id is None
        assert context.timestamp is not None
        assert context.additional_data == {}


class TestGlobalErrorHandler:
    """Test global error handler functions."""
    
    def test_global_error_handler_singleton(self):
        """Test that global error handler returns the same instance."""
        from protein_data_collector.errors import get_error_handler, set_error_handler
        
        handler1 = get_error_handler()
        handler2 = get_error_handler()
        
        assert handler1 is handler2
        
        # Test setting custom handler
        custom_handler = ErrorHandler()
        set_error_handler(custom_handler)
        
        handler3 = get_error_handler()
        assert handler3 is custom_handler
    
    def test_global_handle_error_function(self):
        """Test global handle_error convenience function."""
        context = ErrorContext(operation="test_operation", database="TestDB")
        error = NetworkError("Test error")
        
        error_info = handle_error(error, context)
        
        assert error_info.category == ErrorCategory.NETWORK
        assert error_info.severity == ErrorSeverity.MEDIUM
        assert error_info.action == ErrorAction.RETRY


class TestCustomExceptions:
    """Test custom exception classes."""
    
    def test_network_error_properties(self):
        """Test NetworkError properties."""
        context = ErrorContext(operation="test_op")
        original_error = OSError("Original error")  # Use OSError instead of ConnectionError
        
        error = NetworkError("Network failed", context=context, original_exception=original_error)
        
        assert error.category == ErrorCategory.NETWORK
        assert error.severity == ErrorSeverity.MEDIUM
        assert error.context is context
        assert error.original_exception is original_error
        assert str(error) == "Network failed"
    
    def test_api_error_with_status_code(self):
        """Test APIError with status code."""
        context = ErrorContext(operation="test_op")
        
        error = APIError("API failed", status_code=404, context=context)
        
        assert error.category == ErrorCategory.API
        assert error.severity == ErrorSeverity.MEDIUM
        assert error.status_code == 404
        assert error.context is context
    
    def test_validation_error_properties(self):
        """Test ValidationError properties."""
        context = ErrorContext(operation="validation_op")
        
        error = ValidationError("Validation failed", context=context)
        
        assert error.category == ErrorCategory.VALIDATION
        assert error.severity == ErrorSeverity.HIGH
        assert error.context is context
    
    def test_configuration_error_critical_severity(self):
        """Test that ConfigurationError has critical severity."""
        context = ErrorContext(operation="config_op")
        
        error = ConfigurationError("Config failed", context=context)
        
        assert error.category == ErrorCategory.CONFIGURATION
        assert error.severity == ErrorSeverity.CRITICAL
        assert error.context is context