"""
Unit tests for the logging configuration system.
"""

import pytest
import logging
import json
import tempfile
import os
from io import StringIO
from protein_data_collector.logging_config import (
    JSONFormatter, PerformanceFilter, setup_logging, get_logger,
    log_api_call, log_collection_progress, ContextualFormatter, log_database_operation
)
from protein_data_collector.config import LoggingConfig


class TestJSONFormatter:
    """Test JSON log formatter functionality."""
    
    def test_basic_formatting(self):
        """Test basic JSON log formatting."""
        formatter = JSONFormatter()
        record = logging.LogRecord(
            name="test_logger",
            level=logging.INFO,
            pathname="/test/path.py",
            lineno=42,
            msg="Test message",
            args=(),
            exc_info=None
        )
        record.module = "test_module"
        record.funcName = "test_function"
        
        formatted = formatter.format(record)
        log_data = json.loads(formatted)
        
        assert log_data["level"] == "INFO"
        assert log_data["logger"] == "test_logger"
        assert log_data["message"] == "Test message"
        assert log_data["module"] == "test_module"
        assert log_data["function"] == "test_function"
        assert log_data["line"] == 42
        assert "timestamp" in log_data
    
    def test_formatting_with_extra_fields(self):
        """Test JSON formatting with extra fields."""
        formatter = JSONFormatter()
        record = logging.LogRecord(
            name="test_logger",
            level=logging.INFO,
            pathname="/test/path.py",
            lineno=42,
            msg="Test message",
            args=(),
            exc_info=None
        )
        record.module = "test_module"
        record.funcName = "test_function"
        record.api_method = "GET"
        record.api_url = "https://example.com"
        
        formatted = formatter.format(record)
        log_data = json.loads(formatted)
        
        assert "extra" in log_data
        assert log_data["extra"]["api_method"] == "GET"
        assert log_data["extra"]["api_url"] == "https://example.com"
    
    def test_formatting_with_exception(self):
        """Test JSON formatting with exception information."""
        formatter = JSONFormatter()
        
        try:
            raise ValueError("Test exception")
        except ValueError:
            import sys
            exc_info = sys.exc_info()
        
        record = logging.LogRecord(
            name="test_logger",
            level=logging.ERROR,
            pathname="/test/path.py",
            lineno=42,
            msg="Error occurred",
            args=(),
            exc_info=exc_info
        )
        record.module = "test_module"
        record.funcName = "test_function"
        
        formatted = formatter.format(record)
        log_data = json.loads(formatted)
        
        assert "exception" in log_data
        assert "ValueError: Test exception" in log_data["exception"]


class TestPerformanceFilter:
    """Test performance filter functionality."""
    
    def test_performance_metrics_addition(self):
        """Test adding performance metrics to log records."""
        filter_obj = PerformanceFilter()
        
        record = logging.LogRecord(
            name="test_logger",
            level=logging.INFO,
            pathname="/test/path.py",
            lineno=42,
            msg="Test message",
            args=(),
            exc_info=None
        )
        
        result = filter_obj.filter(record)
        
        assert result is True
        assert hasattr(record, 'cpu_percent')
        assert hasattr(record, 'memory_mb')
        assert hasattr(record, 'uptime_seconds')
        assert hasattr(record, 'process_id')
        assert hasattr(record, 'iso_timestamp')
    
    def test_filter_resilience(self):
        """Test filter handles errors gracefully."""
        filter_obj = PerformanceFilter()
        
        record = logging.LogRecord(
            name="test_logger",
            level=logging.INFO,
            pathname="/test/path.py",
            lineno=42,
            msg="Test message",
            args=(),
            exc_info=None
        )
        
        # Should not raise an exception even if performance metrics fail
        result = filter_obj.filter(record)
        assert result is True


class TestLoggingSetup:
    """Test logging system setup."""
    
    def test_setup_logging_basic(self):
        """Test basic logging setup."""
        config = LoggingConfig(level="DEBUG", format="json")
        setup_logging(config)
        
        root_logger = logging.getLogger()
        assert root_logger.level == logging.DEBUG
        assert len(root_logger.handlers) >= 1
    
    def test_setup_logging_with_file(self):
        """Test logging setup with file handler."""
        with tempfile.NamedTemporaryFile(delete=False) as f:
            log_file = f.name
        
        try:
            config = LoggingConfig(
                level="INFO",
                format="json",
                log_file=log_file,
                max_file_size_mb=1,
                backup_count=3
            )
            setup_logging(config)
            
            root_logger = logging.getLogger()
            assert len(root_logger.handlers) >= 2  # Console + file
            
            # Test logging to file
            test_logger = logging.getLogger("test")
            test_logger.info("Test message")
            
            # Verify file was written
            with open(log_file, 'r') as f:
                content = f.read()
                assert "Test message" in content
        finally:
            if os.path.exists(log_file):
                os.unlink(log_file)
    
    def test_get_logger_basic(self):
        """Test getting logger."""
        logger = get_logger("test_logger")
        
        assert logger.name == "test_logger"
        assert isinstance(logger, logging.Logger)


class TestLoggingUtilities:
    """Test logging utility functions."""
    
    def test_log_api_call_success(self, caplog):
        """Test logging successful API call."""
        logger = logging.getLogger("test_api")
        
        with caplog.at_level(logging.INFO):
            log_api_call(
                logger=logger,
                api_name="test_api",
                endpoint="/api/test",
                method="GET",
                status_code=200,
                duration=0.5
            )
        
        assert len(caplog.records) == 1
        record = caplog.records[0]
        assert "API call: test_api GET /api/test" in record.message
    
    def test_log_api_call_error(self, caplog):
        """Test logging failed API call."""
        logger = logging.getLogger("test_api")
        
        with caplog.at_level(logging.WARNING):
            log_api_call(
                logger=logger,
                api_name="test_api",
                endpoint="/api/test",
                method="POST",
                status_code=500,
                duration=1.0
            )
        
        assert len(caplog.records) == 1
        record = caplog.records[0]
        assert "API call: test_api POST /api/test" in record.message
    
    def test_log_database_operation_success(self, caplog):
        """Test logging successful database operation."""
        logger = logging.getLogger("test_db")
        
        with caplog.at_level(logging.INFO):
            log_database_operation(
                logger=logger,
                operation="INSERT",
                table="proteins",
                count=10,
                duration=0.2
            )
        
        assert len(caplog.records) == 1
        record = caplog.records[0]
        assert "Database INSERT on proteins" in record.message
    
    def test_log_database_operation_error(self, caplog):
        """Test logging failed database operation."""
        logger = logging.getLogger("test_db")
        
        with caplog.at_level(logging.ERROR):
            log_database_operation(
                logger=logger,
                operation="UPDATE",
                table="pfam_families",
                error="Constraint violation"
            )
        
        assert len(caplog.records) == 1
        record = caplog.records[0]
        assert "Database UPDATE on pfam_families failed" in record.message
    
    def test_log_collection_progress(self, caplog):
        """Test logging collection progress."""
        logger = logging.getLogger("test_collection")
        
        with caplog.at_level(logging.INFO):
            log_collection_progress(
                logger=logger,
                phase="proteins",
                progress={"current": 25, "total": 100, "entity_type": "InterProProtein"}
            )
        
        assert len(caplog.records) == 1
        record = caplog.records[0]
        assert "Collection progress: proteins" in record.message