"""
Retry controller and mechanisms for handling external API failures.

This module provides configurable retry logic with exponential backoff
for all external database API calls (InterPro, UniProt).
"""

import asyncio
import logging
import time
from functools import wraps
from typing import Any, Callable, Optional, TypeVar, Union
from dataclasses import dataclass

from .config import RetryConfig
from .errors import ErrorHandler, ErrorContext, ErrorAction, handle_error

# Type variable for generic function return types
T = TypeVar('T')

logger = logging.getLogger(__name__)


@dataclass
class RetryAttempt:
    """Information about a retry attempt."""
    attempt_number: int
    delay: float
    error: Exception
    database: str
    operation: str


class RetryController:
    """
    Configurable retry controller with exponential backoff logic.
    
    Handles retry logic for all external database API calls with
    configurable parameters and comprehensive logging.
    """
    
    def __init__(self, config: RetryConfig):
        """
        Initialize retry controller with configuration.
        
        Args:
            config: RetryConfig instance with retry parameters
        """
        self.config = config
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
    
    def calculate_delay(self, attempt: int) -> float:
        """
        Calculate delay for retry attempt using exponential backoff.
        
        Args:
            attempt: Current attempt number (0-based)
            
        Returns:
            Delay in seconds, capped at max_delay
        """
        delay = self.config.initial_delay * (self.config.backoff_multiplier ** attempt)
        return min(delay, self.config.max_delay)
    
    def log_retry_attempt(self, database: str, operation: str, attempt: int, error: Exception, delay: float) -> None:
        """
        Log retry attempt with contextual information.
        
        Args:
            database: Name of the database being accessed
            operation: Description of the operation being retried
            attempt: Current attempt number (1-based for logging)
            error: Exception that caused the retry
            delay: Delay before next attempt
        """
        self.logger.warning(
            "Retry attempt for %s operation on %s",
            operation,
            database,
            extra={
                "database": database,
                "operation": operation,
                "attempt_number": attempt,
                "max_retries": self.config.max_retries,
                "delay_seconds": delay,
                "error_type": type(error).__name__,
                "error_message": str(error)
            }
        )
    
    def log_final_failure(self, database: str, operation: str, error: Exception) -> None:
        """
        Log final failure after all retries exhausted.
        
        Args:
            database: Name of the database being accessed
            operation: Description of the failed operation
            error: Final exception that caused the failure
        """
        self.logger.error(
            "Final failure for %s operation on %s after %d retries",
            operation,
            database,
            self.config.max_retries,
            extra={
                "database": database,
                "operation": operation,
                "max_retries": self.config.max_retries,
                "error_type": type(error).__name__,
                "error_message": str(error),
                "final_failure": True
            }
        )
    
    def execute_with_retry(
        self,
        operation: Callable[[], T],
        database: str,
        operation_name: str,
        error_handler: Optional[ErrorHandler] = None
    ) -> T:
        """
        Execute operation with retry logic (synchronous version).
        
        Args:
            operation: Function to execute with retry logic
            database: Name of the database being accessed
            operation_name: Description of the operation for logging
            error_handler: Optional error handler for classification
            
        Returns:
            Result of the operation
            
        Raises:
            Exception: The last exception if all retries are exhausted
        """
        last_exception = None
        handler = error_handler or ErrorHandler()
        
        for attempt in range(self.config.max_retries + 1):
            try:
                return operation()
            except Exception as e:
                last_exception = e
                
                # Create error context for classification
                context = ErrorContext(
                    operation=operation_name,
                    database=database,
                    additional_data={"attempt": attempt + 1, "max_retries": self.config.max_retries}
                )
                
                # Classify error to determine if we should retry
                error_info = handler.classify_error(e, context)
                
                if attempt == self.config.max_retries:
                    # Final attempt failed, log and re-raise
                    self.log_final_failure(database, operation_name, e)
                    raise e
                
                # Check if error should be retried
                if error_info.action not in [ErrorAction.RETRY, ErrorAction.FALLBACK]:
                    # Error should not be retried, re-raise immediately
                    self.logger.info(
                        "Error classified as non-retryable for %s operation on %s",
                        operation_name,
                        database,
                        extra={
                            "database": database,
                            "operation": operation_name,
                            "error_category": error_info.category.value,
                            "recommended_action": error_info.action.value,
                            "attempt": attempt + 1
                        }
                    )
                    raise e
                
                # Calculate delay and log retry attempt
                delay = self.calculate_delay(attempt)
                self.log_retry_attempt(database, operation_name, attempt + 1, e, delay)
                
                # Wait before retry
                time.sleep(delay)
        
        # This should never be reached, but included for completeness
        if last_exception:
            raise last_exception
        raise RuntimeError("Unexpected retry loop exit")
    
    async def execute_with_retry_async(
        self,
        operation: Callable[[], T],
        database: str,
        operation_name: str,
        error_handler: Optional[ErrorHandler] = None
    ) -> T:
        """
        Execute operation with retry logic (asynchronous version).
        
        Args:
            operation: Async function to execute with retry logic
            database: Name of the database being accessed
            operation_name: Description of the operation for logging
            error_handler: Optional error handler for classification
            
        Returns:
            Result of the operation
            
        Raises:
            Exception: The last exception if all retries are exhausted
        """
        last_exception = None
        handler = error_handler or ErrorHandler()
        
        for attempt in range(self.config.max_retries + 1):
            try:
                if asyncio.iscoroutinefunction(operation):
                    return await operation()
                else:
                    return operation()
            except Exception as e:
                last_exception = e
                
                # Create error context for classification
                context = ErrorContext(
                    operation=operation_name,
                    database=database,
                    additional_data={"attempt": attempt + 1, "max_retries": self.config.max_retries}
                )
                
                # Classify error to determine if we should retry
                error_info = handler.classify_error(e, context)
                
                if attempt == self.config.max_retries:
                    # Final attempt failed, log and re-raise
                    self.log_final_failure(database, operation_name, e)
                    raise e
                
                # Check if error should be retried
                if error_info.action not in [ErrorAction.RETRY, ErrorAction.FALLBACK]:
                    # Error should not be retried, re-raise immediately
                    self.logger.info(
                        "Error classified as non-retryable for %s operation on %s",
                        operation_name,
                        database,
                        extra={
                            "database": database,
                            "operation": operation_name,
                            "error_category": error_info.category.value,
                            "recommended_action": error_info.action.value,
                            "attempt": attempt + 1
                        }
                    )
                    raise e
                
                # Calculate delay and log retry attempt
                delay = self.calculate_delay(attempt)
                self.log_retry_attempt(database, operation_name, attempt + 1, e, delay)
                
                # Wait before retry (async)
                await asyncio.sleep(delay)
        
        # This should never be reached, but included for completeness
        if last_exception:
            raise last_exception
        raise RuntimeError("Unexpected retry loop exit")


def with_retry(
    database: str,
    operation_name: str,
    config: Optional[RetryConfig] = None
):
    """
    Decorator for adding retry logic to functions.
    
    Args:
        database: Name of the database being accessed
        operation_name: Description of the operation for logging
        config: Optional RetryConfig, uses default if not provided
        
    Returns:
        Decorated function with retry logic
    """
    def decorator(func: Callable[..., T]) -> Callable[..., T]:
        @wraps(func)
        def wrapper(*args, **kwargs) -> T:
            # Import here to avoid circular imports
            from .config import get_config
            
            retry_config = config or get_config().retry
            controller = RetryController(retry_config)
            
            def operation():
                return func(*args, **kwargs)
            
            return controller.execute_with_retry(operation, database, operation_name)
        
        return wrapper
    return decorator


def with_retry_async(
    database: str,
    operation_name: str,
    config: Optional[RetryConfig] = None
):
    """
    Decorator for adding retry logic to async functions.
    
    Args:
        database: Name of the database being accessed
        operation_name: Description of the operation for logging
        config: Optional RetryConfig, uses default if not provided
        
    Returns:
        Decorated async function with retry logic
    """
    def decorator(func: Callable[..., T]) -> Callable[..., T]:
        @wraps(func)
        async def wrapper(*args, **kwargs) -> T:
            # Import here to avoid circular imports
            from .config import get_config
            
            retry_config = config or get_config().retry
            controller = RetryController(retry_config)
            
            async def operation():
                return await func(*args, **kwargs)
            
            return await controller.execute_with_retry_async(operation, database, operation_name)
        
        return wrapper
    return decorator


# Global retry controller instance
_retry_controller: Optional[RetryController] = None


def get_retry_controller() -> RetryController:
    """Get the global retry controller instance."""
    global _retry_controller
    if _retry_controller is None:
        from .config import get_config
        _retry_controller = RetryController(get_config().retry)
    return _retry_controller


def set_retry_controller(controller: RetryController) -> None:
    """Set the global retry controller instance."""
    global _retry_controller
    _retry_controller = controller