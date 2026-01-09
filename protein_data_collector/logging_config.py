"""
Logging configuration for the Protein Data Collector system.

This module provides structured logging with performance metrics, contextual information,
and configurable output formats for comprehensive monitoring and debugging.
"""

import logging
import logging.handlers
import json
import sys
import time
from datetime import datetime
from typing import Dict, Any, Optional
from pathlib import Path
import psutil
import os

from .config import LoggingConfig


class PerformanceFilter(logging.Filter):
    """Filter to add performance metrics to log records."""
    
    def __init__(self):
        super().__init__()
        self.process = psutil.Process()
        self.start_time = time.time()
    
    def filter(self, record):
        """Add performance metrics to the log record."""
        try:
            # Add performance metrics
            record.cpu_percent = self.process.cpu_percent()
            record.memory_mb = self.process.memory_info().rss / 1024 / 1024
            record.uptime_seconds = time.time() - self.start_time
            
            # Add process information
            record.process_id = os.getpid()
            record.thread_id = record.thread
            
            # Add timestamp in ISO format
            record.iso_timestamp = datetime.fromtimestamp(record.created).isoformat()
            
        except Exception:
            # Don't fail logging if performance metrics can't be collected
            pass
        
        return True


class JSONFormatter(logging.Formatter):
    """JSON formatter for structured logging."""
    
    def __init__(self, include_performance=True):
        super().__init__()
        self.include_performance = include_performance
    
    def format(self, record):
        """Format log record as JSON."""
        log_entry = {
            "timestamp": datetime.fromtimestamp(record.created).isoformat(),
            "level": record.levelname,
            "logger": record.name,
            "message": record.getMessage(),
            "module": record.module,
            "function": record.funcName,
            "line": record.lineno,
        }
        
        # Add exception information if present
        if record.exc_info:
            log_entry["exception"] = self.formatException(record.exc_info)
        
        # Add extra fields from the record
        extra_fields = {}
        for key, value in record.__dict__.items():
            if key not in [
                'name', 'msg', 'args', 'levelname', 'levelno', 'pathname',
                'filename', 'module', 'exc_info', 'exc_text', 'stack_info',
                'lineno', 'funcName', 'created', 'msecs', 'relativeCreated',
                'thread', 'threadName', 'processName', 'process', 'message'
            ]:
                extra_fields[key] = value
        
        if extra_fields:
            log_entry["extra"] = extra_fields
        
        # Add performance metrics if enabled
        if self.include_performance and hasattr(record, 'cpu_percent'):
            log_entry["performance"] = {
                "cpu_percent": getattr(record, 'cpu_percent', 0),
                "memory_mb": getattr(record, 'memory_mb', 0),
                "uptime_seconds": getattr(record, 'uptime_seconds', 0),
                "process_id": getattr(record, 'process_id', 0),
                "thread_id": getattr(record, 'thread_id', 0)
            }
        
        return json.dumps(log_entry, default=str)


class ContextualFormatter(logging.Formatter):
    """Human-readable formatter with contextual information."""
    
    def __init__(self, include_performance=True):
        super().__init__()
        self.include_performance = include_performance
        
        # Define format string
        format_str = (
            "%(iso_timestamp)s [%(levelname)s] %(name)s:%(funcName)s:%(lineno)d - %(message)s"
        )
        
        if include_performance:
            format_str += " [CPU: %(cpu_percent).1f%% MEM: %(memory_mb).1fMB]"
        
        self._formatter = logging.Formatter(format_str)
    
    def format(self, record):
        """Format log record with contextual information."""
        # Ensure required attributes exist
        if not hasattr(record, 'iso_timestamp'):
            record.iso_timestamp = datetime.fromtimestamp(record.created).isoformat()
        if not hasattr(record, 'cpu_percent'):
            record.cpu_percent = 0.0
        if not hasattr(record, 'memory_mb'):
            record.memory_mb = 0.0
        
        return self._formatter.format(record)


class MetricsHandler(logging.Handler):
    """Handler to collect logging metrics."""
    
    def __init__(self):
        super().__init__()
        self.metrics = {
            'total_logs': 0,
            'error_count': 0,
            'warning_count': 0,
            'info_count': 0,
            'debug_count': 0,
            'last_error': None,
            'last_warning': None,
            'start_time': time.time()
        }
    
    def emit(self, record):
        """Collect metrics from log records."""
        self.metrics['total_logs'] += 1
        
        if record.levelno >= logging.ERROR:
            self.metrics['error_count'] += 1
            self.metrics['last_error'] = {
                'timestamp': datetime.fromtimestamp(record.created).isoformat(),
                'message': record.getMessage(),
                'logger': record.name,
                'module': record.module,
                'function': record.funcName,
                'line': record.lineno
            }
        elif record.levelno >= logging.WARNING:
            self.metrics['warning_count'] += 1
            self.metrics['last_warning'] = {
                'timestamp': datetime.fromtimestamp(record.created).isoformat(),
                'message': record.getMessage(),
                'logger': record.name,
                'module': record.module,
                'function': record.funcName,
                'line': record.lineno
            }
        elif record.levelno >= logging.INFO:
            self.metrics['info_count'] += 1
        else:
            self.metrics['debug_count'] += 1
    
    def get_metrics(self) -> Dict[str, Any]:
        """Get collected metrics."""
        uptime = time.time() - self.metrics['start_time']
        return {
            **self.metrics,
            'uptime_seconds': uptime,
            'logs_per_second': self.metrics['total_logs'] / uptime if uptime > 0 else 0,
            'error_rate': self.metrics['error_count'] / self.metrics['total_logs'] if self.metrics['total_logs'] > 0 else 0
        }


# Global metrics handler instance
_metrics_handler = None


def get_logging_metrics() -> Dict[str, Any]:
    """Get logging metrics."""
    global _metrics_handler
    if _metrics_handler:
        return _metrics_handler.get_metrics()
    return {}


def setup_logging(config: LoggingConfig) -> None:
    """
    Set up comprehensive logging configuration.
    
    Args:
        config: Logging configuration object
    """
    global _metrics_handler
    
    # Create logs directory if it doesn't exist
    if config.log_file:
        log_path = Path(config.log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(getattr(logging, config.level.upper()))
    
    # Clear existing handlers
    root_logger.handlers.clear()
    
    # Create performance filter
    perf_filter = PerformanceFilter()
    
    # Create metrics handler
    _metrics_handler = MetricsHandler()
    root_logger.addHandler(_metrics_handler)
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(getattr(logging, config.level.upper()))
    console_handler.addFilter(perf_filter)
    
    if config.format.lower() == 'json' and config.structured:
        console_formatter = JSONFormatter(include_performance=True)
    else:
        console_formatter = ContextualFormatter(include_performance=True)
    
    console_handler.setFormatter(console_formatter)
    root_logger.addHandler(console_handler)
    
    # File handler (if configured)
    if config.log_file:
        # Parse file size
        max_bytes = _parse_file_size(f"{config.max_file_size_mb}MB")
        
        file_handler = logging.handlers.RotatingFileHandler(
            filename=config.log_file,
            maxBytes=max_bytes,
            backupCount=config.backup_count,
            encoding='utf-8'
        )
        file_handler.setLevel(getattr(logging, config.level.upper()))
        file_handler.addFilter(perf_filter)
        
        if config.format.lower() == 'json' and config.structured:
            file_formatter = JSONFormatter(include_performance=True)
        else:
            file_formatter = ContextualFormatter(include_performance=True)
        
        file_handler.setFormatter(file_formatter)
        root_logger.addHandler(file_handler)
    
    # Configure specific loggers
    _configure_library_loggers()
    
    # Log startup message
    logger = logging.getLogger(__name__)
    logger.info(
        "Logging system initialized",
        extra={
            "log_level": config.level,
            "log_format": config.format,
            "structured": config.structured,
            "file_logging": bool(config.log_file),
            "performance_metrics": True
        }
    )


def _parse_file_size(size_str: str) -> int:
    """Parse file size string to bytes."""
    size_str = size_str.upper().strip()
    
    if size_str.endswith('KB'):
        return int(size_str[:-2]) * 1024
    elif size_str.endswith('MB'):
        return int(size_str[:-2]) * 1024 * 1024
    elif size_str.endswith('GB'):
        return int(size_str[:-2]) * 1024 * 1024 * 1024
    else:
        # Assume bytes
        return int(size_str)


def _configure_library_loggers():
    """Configure logging levels for third-party libraries."""
    # Reduce noise from third-party libraries
    logging.getLogger('httpx').setLevel(logging.WARNING)
    logging.getLogger('httpcore').setLevel(logging.WARNING)
    logging.getLogger('urllib3').setLevel(logging.WARNING)
    logging.getLogger('mysql.connector').setLevel(logging.WARNING)
    logging.getLogger('sqlalchemy').setLevel(logging.WARNING)
    logging.getLogger('uvicorn').setLevel(logging.INFO)
    logging.getLogger('fastapi').setLevel(logging.INFO)


def get_logger(name: str) -> logging.Logger:
    """
    Get a logger with the specified name.
    
    Args:
        name: Logger name (typically __name__)
        
    Returns:
        Configured logger instance
    """
    return logging.getLogger(name)


def log_performance_metrics(logger: logging.Logger, operation: str, duration: float, **kwargs):
    """
    Log performance metrics for an operation.
    
    Args:
        logger: Logger instance
        operation: Operation name
        duration: Operation duration in seconds
        **kwargs: Additional metrics to log
    """
    logger.info(
        f"Performance: {operation} completed",
        extra={
            "operation": operation,
            "duration_seconds": duration,
            "duration_ms": duration * 1000,
            **kwargs
        }
    )


def log_api_call(logger: logging.Logger, api_name: str, endpoint: str, 
                method: str, status_code: int, duration: float, **kwargs):
    """
    Log API call metrics.
    
    Args:
        logger: Logger instance
        api_name: API name (e.g., 'interpro', 'uniprot')
        endpoint: API endpoint
        method: HTTP method
        status_code: HTTP status code
        duration: Request duration in seconds
        **kwargs: Additional context
    """
    level = logging.INFO if 200 <= status_code < 400 else logging.WARNING
    
    logger.log(
        level,
        f"API call: {api_name} {method} {endpoint}",
        extra={
            "api_name": api_name,
            "endpoint": endpoint,
            "method": method,
            "status_code": status_code,
            "duration_seconds": duration,
            "duration_ms": duration * 1000,
            "success": 200 <= status_code < 400,
            **kwargs
        }
    )


def log_collection_progress(logger: logging.Logger, phase: str, progress: Dict[str, Any]):
    """
    Log data collection progress.
    
    Args:
        logger: Logger instance
        phase: Collection phase name
        progress: Progress information
    """
    logger.info(
        f"Collection progress: {phase}",
        extra={
            "collection_phase": phase,
            "progress": progress
        }
    )


def log_database_operation(logger: logging.Logger, operation: str, table: str, 
                          count: Optional[int] = None, duration: Optional[float] = None, 
                          error: Optional[str] = None, **kwargs):
    """
    Log database operation metrics.
    
    Args:
        logger: Logger instance
        operation: Database operation (INSERT, UPDATE, DELETE, SELECT)
        table: Database table name
        count: Number of records affected
        duration: Operation duration in seconds
        error: Error message if operation failed
        **kwargs: Additional context
    """
    level = logging.ERROR if error else logging.INFO
    
    extra_data = {
        "operation": operation,
        "table": table,
        "database_operation": True,
        **kwargs
    }
    
    if count is not None:
        extra_data["record_count"] = count
    if duration is not None:
        extra_data["duration_seconds"] = duration
        extra_data["duration_ms"] = duration * 1000
    if error:
        extra_data["error"] = error
    
    message = f"Database {operation} on {table}"
    if error:
        message += f" failed: {error}"
    elif count is not None:
        message += f" affected {count} records"
    
    logger.log(level, message, extra=extra_data)


def log_error_with_context(logger: logging.Logger, error: Exception, 
                          operation: str, **context):
    """
    Log error with full context information.
    
    Args:
        logger: Logger instance
        error: Exception that occurred
        operation: Operation that failed
        **context: Additional context information
    """
    logger.error(
        f"Error in {operation}: {str(error)}",
        exc_info=True,
        extra={
            "operation": operation,
            "error_type": type(error).__name__,
            "error_message": str(error),
            **context
        }
    )