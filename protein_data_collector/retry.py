"""Retry decorator using tenacity with exponential backoff."""

import logging
from tenacity import (
    retry,
    stop_after_attempt,
    wait_exponential,
    retry_if_exception_type,
    before_sleep_log,
)
from .errors import NetworkError, APIError

logger = logging.getLogger(__name__)


def with_retry(max_attempts: int = 3, initial_delay: float = 1.0, max_delay: float = 60.0):
    """Retry on transient network/API errors with exponential backoff."""
    return retry(
        stop=stop_after_attempt(max_attempts),
        wait=wait_exponential(multiplier=initial_delay, max=max_delay),
        retry=retry_if_exception_type((NetworkError, APIError)),
        before_sleep=before_sleep_log(logger, logging.WARNING),
        reraise=True,
    )
