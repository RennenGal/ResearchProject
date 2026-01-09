"""
Comprehensive error handling system for the Protein Data Collector.

This module provides error classification, logging, and recovery strategies
for different types of errors encountered during data collection operations.
"""

import logging
from enum import Enum
from typing import Any, Dict, Optional, Union
from dataclasses import dataclass, field
from datetime import datetime


class ErrorCategory(Enum):
    """Categories of errors that can occur in the system."""
    NETWORK = "network"
    API = "api"
    DATA = "data"
    VALIDATION = "validation"
    DATABASE = "database"
    CONFIGURATION = "configuration"
    UNKNOWN = "unknown"


class ErrorSeverity(Enum):
    """Severity levels for errors."""
    LOW = "low"
    MEDIUM = "medium"
    HIGH = "high"
    CRITICAL = "critical"


class ErrorAction(Enum):
    """Actions to take when an error occurs."""
    RETRY = "retry"
    SKIP = "skip"
    FAIL = "fail"
    LOG_AND_CONTINUE = "log_and_continue"
    FALLBACK = "fallback"


@dataclass
class ErrorContext:
    """Contextual information about an error."""
    operation: str
    database: Optional[str] = None
    entity_id: Optional[str] = None
    entity_type: Optional[str] = None
    request_url: Optional[str] = None
    request_params: Optional[Dict[str, Any]] = None
    response_status: Optional[int] = None
    timestamp: datetime = field(default_factory=datetime.now)
    additional_data: Dict[str, Any] = field(default_factory=dict)


@dataclass
class ErrorInfo:
    """Complete information about an error occurrence."""
    category: ErrorCategory
    severity: ErrorSeverity
    action: ErrorAction
    message: str
    original_exception: Exception
    context: ErrorContext
    recovery_suggestions: list[str] = field(default_factory=list)


class ProteinDataCollectorError(Exception):
    """Base exception class for Protein Data Collector errors."""
    
    def __init__(
        self,
        message: str,
        category: ErrorCategory = ErrorCategory.UNKNOWN,
        severity: ErrorSeverity = ErrorSeverity.MEDIUM,
        context: Optional[ErrorContext] = None,
        original_exception: Optional[Exception] = None
    ):
        super().__init__(message)
        self.message = message
        self.category = category
        self.severity = severity
        self.context = context or ErrorContext(operation="unknown")
        self.original_exception = original_exception


class NetworkError(ProteinDataCollectorError):
    """Errors related to network connectivity and timeouts."""
    
    def __init__(self, message: str, context: Optional[ErrorContext] = None, original_exception: Optional[Exception] = None):
        super().__init__(
            message=message,
            category=ErrorCategory.NETWORK,
            severity=ErrorSeverity.MEDIUM,
            context=context,
            original_exception=original_exception
        )


class APIError(ProteinDataCollectorError):
    """Errors related to external API interactions."""
    
    def __init__(
        self,
        message: str,
        status_code: Optional[int] = None,
        context: Optional[ErrorContext] = None,
        original_exception: Optional[Exception] = None
    ):
        super().__init__(
            message=message,
            category=ErrorCategory.API,
            severity=ErrorSeverity.MEDIUM,
            context=context,
            original_exception=original_exception
        )
        self.status_code = status_code


class DataError(ProteinDataCollectorError):
    """Errors related to data processing and format issues."""
    
    def __init__(self, message: str, context: Optional[ErrorContext] = None, original_exception: Optional[Exception] = None):
        super().__init__(
            message=message,
            category=ErrorCategory.DATA,
            severity=ErrorSeverity.HIGH,
            context=context,
            original_exception=original_exception
        )


class ValidationError(ProteinDataCollectorError):
    """Errors related to data validation failures."""
    
    def __init__(self, message: str, context: Optional[ErrorContext] = None, original_exception: Optional[Exception] = None):
        super().__init__(
            message=message,
            category=ErrorCategory.VALIDATION,
            severity=ErrorSeverity.HIGH,
            context=context,
            original_exception=original_exception
        )


class DatabaseError(ProteinDataCollectorError):
    """Errors related to database operations."""
    
    def __init__(self, message: str, context: Optional[ErrorContext] = None, original_exception: Optional[Exception] = None):
        super().__init__(
            message=message,
            category=ErrorCategory.DATABASE,
            severity=ErrorSeverity.HIGH,
            context=context,
            original_exception=original_exception
        )


class ConfigurationError(ProteinDataCollectorError):
    """Errors related to system configuration."""
    
    def __init__(self, message: str, context: Optional[ErrorContext] = None, original_exception: Optional[Exception] = None):
        super().__init__(
            message=message,
            category=ErrorCategory.CONFIGURATION,
            severity=ErrorSeverity.CRITICAL,
            context=context,
            original_exception=original_exception
        )


class ErrorHandler:
    """
    Comprehensive error handler with classification and recovery strategies.
    
    Provides centralized error handling with contextual logging and
    appropriate recovery actions for different error types.
    """
    
    def __init__(self):
        """Initialize error handler with logger."""
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        self._error_counts = {}
    
    def classify_error(self, exception: Exception, context: ErrorContext) -> ErrorInfo:
        """
        Classify an exception and determine appropriate handling strategy.
        
        Args:
            exception: The exception to classify
            context: Contextual information about the error
            
        Returns:
            ErrorInfo with classification and recommended action
        """
        # Handle known custom exceptions
        if isinstance(exception, ProteinDataCollectorError):
            return ErrorInfo(
                category=exception.category,
                severity=exception.severity,
                action=self._determine_action(exception.category, exception.severity),
                message=exception.message,
                original_exception=exception,
                context=context,
                recovery_suggestions=self._get_recovery_suggestions(exception.category)
            )
        
        # Classify standard exceptions
        category, severity = self._classify_standard_exception(exception, context)
        action = self._determine_action(category, severity)
        
        return ErrorInfo(
            category=category,
            severity=severity,
            action=action,
            message=str(exception),
            original_exception=exception,
            context=context,
            recovery_suggestions=self._get_recovery_suggestions(category)
        )
    
    def _classify_standard_exception(self, exception: Exception, context: ErrorContext) -> tuple[ErrorCategory, ErrorSeverity]:
        """Classify standard Python exceptions."""
        exception_type = type(exception).__name__
        
        # Network-related errors
        if exception_type in ['ConnectionError', 'TimeoutError', 'ConnectTimeout', 'ReadTimeout']:
            return ErrorCategory.NETWORK, ErrorSeverity.MEDIUM
        
        # HTTP/API errors
        if exception_type in ['HTTPError', 'RequestException']:
            # Check status code if available
            if hasattr(exception, 'response') and exception.response:
                status_code = exception.response.status_code
                if 400 <= status_code < 500:
                    return ErrorCategory.API, ErrorSeverity.HIGH
                elif status_code >= 500:
                    return ErrorCategory.API, ErrorSeverity.MEDIUM
            return ErrorCategory.API, ErrorSeverity.MEDIUM
        
        # Data processing errors
        if exception_type in ['JSONDecodeError', 'ValueError', 'KeyError', 'IndexError']:
            return ErrorCategory.DATA, ErrorSeverity.HIGH
        
        # Database errors
        if exception_type in ['DatabaseError', 'IntegrityError', 'OperationalError']:
            return ErrorCategory.DATABASE, ErrorSeverity.HIGH
        
        # Configuration errors
        if exception_type in ['FileNotFoundError', 'PermissionError'] and 'config' in str(exception).lower():
            return ErrorCategory.CONFIGURATION, ErrorSeverity.CRITICAL
        
        # Default classification
        return ErrorCategory.UNKNOWN, ErrorSeverity.MEDIUM
    
    def _determine_action(self, category: ErrorCategory, severity: ErrorSeverity) -> ErrorAction:
        """Determine appropriate action based on error category and severity."""
        if category == ErrorCategory.NETWORK:
            return ErrorAction.RETRY
        elif category == ErrorCategory.API:
            if severity == ErrorSeverity.HIGH:
                return ErrorAction.SKIP  # Client errors (4xx) - skip problematic requests
            else:
                return ErrorAction.RETRY  # Server errors (5xx) - retry
        elif category == ErrorCategory.DATA:
            return ErrorAction.SKIP  # Skip malformed data
        elif category == ErrorCategory.VALIDATION:
            return ErrorAction.SKIP  # Skip invalid data
        elif category == ErrorCategory.DATABASE:
            if severity == ErrorSeverity.CRITICAL:
                return ErrorAction.FAIL  # Critical DB issues should fail
            else:
                return ErrorAction.RETRY  # Transient DB issues can be retried
        elif category == ErrorCategory.CONFIGURATION:
            return ErrorAction.FAIL  # Configuration errors should fail immediately
        else:
            return ErrorAction.LOG_AND_CONTINUE  # Unknown errors - log and continue
    
    def _get_recovery_suggestions(self, category: ErrorCategory) -> list[str]:
        """Get recovery suggestions for different error categories."""
        suggestions = {
            ErrorCategory.NETWORK: [
                "Check network connectivity",
                "Verify API endpoints are accessible",
                "Consider increasing timeout values",
                "Check for firewall or proxy issues"
            ],
            ErrorCategory.API: [
                "Verify API credentials and permissions",
                "Check API rate limits and quotas",
                "Validate request parameters and format",
                "Check API documentation for changes"
            ],
            ErrorCategory.DATA: [
                "Validate data format and structure",
                "Check for missing required fields",
                "Verify data encoding and character sets",
                "Review data transformation logic"
            ],
            ErrorCategory.VALIDATION: [
                "Review validation rules and constraints",
                "Check data quality at source",
                "Verify field formats and ranges",
                "Consider data cleaning procedures"
            ],
            ErrorCategory.DATABASE: [
                "Check database connectivity and credentials",
                "Verify database schema and constraints",
                "Monitor database performance and resources",
                "Review transaction isolation levels"
            ],
            ErrorCategory.CONFIGURATION: [
                "Verify configuration file format and syntax",
                "Check file permissions and accessibility",
                "Validate configuration values and ranges",
                "Review environment variable settings"
            ]
        }
        return suggestions.get(category, ["Review error details and system logs"])
    
    def handle_error(self, exception: Exception, context: ErrorContext) -> ErrorInfo:
        """
        Handle an error with appropriate logging and classification.
        
        Args:
            exception: The exception to handle
            context: Contextual information about the error
            
        Returns:
            ErrorInfo with handling details
        """
        error_info = self.classify_error(exception, context)
        
        # Track error counts for monitoring
        error_key = f"{error_info.category.value}:{error_info.context.operation}"
        self._error_counts[error_key] = self._error_counts.get(error_key, 0) + 1
        
        # Log error with appropriate level
        self._log_error(error_info)
        
        return error_info
    
    def _log_error(self, error_info: ErrorInfo) -> None:
        """Log error with contextual information."""
        log_data = {
            "error_category": error_info.category.value,
            "error_severity": error_info.severity.value,
            "recommended_action": error_info.action.value,
            "operation": error_info.context.operation,
            "database": error_info.context.database,
            "entity_id": error_info.context.entity_id,
            "entity_type": error_info.context.entity_type,
            "request_url": error_info.context.request_url,
            "response_status": error_info.context.response_status,
            "timestamp": error_info.context.timestamp.isoformat(),
            "exception_type": type(error_info.original_exception).__name__,
            "exception_message": str(error_info.original_exception),
            "recovery_suggestions": error_info.recovery_suggestions
        }
        
        # Add additional context data
        if error_info.context.additional_data:
            log_data.update(error_info.context.additional_data)
        
        # Log at appropriate level based on severity
        if error_info.severity == ErrorSeverity.CRITICAL:
            self.logger.critical("Critical error occurred: %s", error_info.message, extra=log_data)
        elif error_info.severity == ErrorSeverity.HIGH:
            self.logger.error("High severity error: %s", error_info.message, extra=log_data)
        elif error_info.severity == ErrorSeverity.MEDIUM:
            self.logger.warning("Medium severity error: %s", error_info.message, extra=log_data)
        else:
            self.logger.info("Low severity error: %s", error_info.message, extra=log_data)
    
    def get_error_statistics(self) -> Dict[str, int]:
        """Get error count statistics."""
        return self._error_counts.copy()
    
    def reset_error_statistics(self) -> None:
        """Reset error count statistics."""
        self._error_counts.clear()


# Global error handler instance
_error_handler: Optional[ErrorHandler] = None


def get_error_handler() -> ErrorHandler:
    """Get the global error handler instance."""
    global _error_handler
    if _error_handler is None:
        _error_handler = ErrorHandler()
    return _error_handler


def set_error_handler(handler: ErrorHandler) -> None:
    """Set the global error handler instance."""
    global _error_handler
    _error_handler = handler


def handle_error(exception: Exception, context: ErrorContext) -> ErrorInfo:
    """
    Convenience function to handle errors using the global error handler.
    
    Args:
        exception: The exception to handle
        context: Contextual information about the error
        
    Returns:
        ErrorInfo with handling details
    """
    return get_error_handler().handle_error(exception, context)


def create_error_context(
    operation: str,
    database: Optional[str] = None,
    entity_id: Optional[str] = None,
    entity_type: Optional[str] = None,
    request_url: Optional[str] = None,
    request_params: Optional[Dict[str, Any]] = None,
    response_status: Optional[int] = None,
    **additional_data
) -> ErrorContext:
    """
    Convenience function to create error context.
    
    Args:
        operation: Name of the operation being performed
        database: Name of the database being accessed
        entity_id: ID of the entity being processed
        entity_type: Type of the entity being processed
        request_url: URL of the request that failed
        request_params: Parameters of the request
        response_status: HTTP response status code
        **additional_data: Additional contextual data
        
    Returns:
        ErrorContext instance
    """
    return ErrorContext(
        operation=operation,
        database=database,
        entity_id=entity_id,
        entity_type=entity_type,
        request_url=request_url,
        request_params=request_params,
        response_status=response_status,
        additional_data=additional_data
    )