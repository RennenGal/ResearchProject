"""
Comprehensive rate limiting system for external API access.

This module provides configurable rate limiting with exponential backoff
for rate limit violations, monitoring, and reporting capabilities.
"""

import asyncio
import logging
import time
from dataclasses import dataclass, field
from typing import Dict, Optional, Any, List
from enum import Enum
from collections import defaultdict, deque


class RateLimitViolationType(Enum):
    """Types of rate limit violations."""
    SOFT_LIMIT = "soft_limit"  # Approaching rate limit
    HARD_LIMIT = "hard_limit"  # Rate limit exceeded
    BURST_LIMIT = "burst_limit"  # Too many requests in short time


@dataclass
class RateLimitConfig:
    """Configuration for rate limiting behavior."""
    requests_per_second: float = 10.0
    burst_limit: int = 50  # Maximum requests in burst window
    burst_window_seconds: int = 60  # Burst window duration
    
    # Exponential backoff for rate limit violations
    violation_initial_delay: float = 1.0
    violation_backoff_multiplier: float = 2.0
    violation_max_delay: float = 300.0  # 5 minutes max
    
    # Monitoring thresholds
    soft_limit_threshold: float = 0.8  # Warn at 80% of rate limit
    enable_monitoring: bool = True
    enable_reporting: bool = True


@dataclass
class RateLimitViolation:
    """Information about a rate limit violation."""
    violation_type: RateLimitViolationType
    timestamp: float
    api_name: str
    current_rate: float
    limit_rate: float
    delay_applied: float
    requests_in_window: int = 0


@dataclass
class RateLimitStats:
    """Statistics for rate limiting monitoring."""
    api_name: str
    total_requests: int = 0
    total_violations: int = 0
    current_rate: float = 0.0
    average_delay: float = 0.0
    last_violation: Optional[RateLimitViolation] = None
    violations_by_type: Dict[RateLimitViolationType, int] = field(default_factory=lambda: defaultdict(int))


class TokenBucket:
    """Token bucket algorithm implementation for rate limiting."""
    
    def __init__(self, rate: float, capacity: Optional[int] = None):
        """
        Initialize token bucket.
        
        Args:
            rate: Tokens per second
            capacity: Maximum tokens in bucket (defaults to rate * 2)
        """
        self.rate = rate
        self.capacity = capacity or max(int(rate * 2), 1)
        self.tokens = float(self.capacity)
        self.last_update = time.time()
        self._lock = asyncio.Lock()
    
    async def acquire(self, tokens: int = 1) -> float:
        """
        Acquire tokens from bucket, waiting if necessary.
        
        Args:
            tokens: Number of tokens to acquire
            
        Returns:
            Delay time in seconds (0 if no delay needed)
        """
        async with self._lock:
            now = time.time()
            
            # Add tokens based on elapsed time
            elapsed = now - self.last_update
            self.tokens = min(self.capacity, self.tokens + elapsed * self.rate)
            self.last_update = now
            
            if self.tokens >= tokens:
                # Sufficient tokens available
                self.tokens -= tokens
                return 0.0
            else:
                # Need to wait for tokens
                needed_tokens = tokens - self.tokens
                delay = needed_tokens / self.rate
                self.tokens = 0.0  # Use all available tokens
                return delay
    
    def get_current_tokens(self) -> float:
        """Get current number of tokens in bucket."""
        now = time.time()
        elapsed = now - self.last_update
        return min(self.capacity, self.tokens + elapsed * self.rate)


class RateLimitMonitor:
    """Monitor and track rate limiting statistics."""
    
    def __init__(self, enable_monitoring: bool = True):
        """
        Initialize rate limit monitor.
        
        Args:
            enable_monitoring: Whether to collect monitoring data
        """
        self.enable_monitoring = enable_monitoring
        self.stats: Dict[str, RateLimitStats] = {}
        self.violations: List[RateLimitViolation] = []
        self.request_history: Dict[str, deque] = defaultdict(lambda: deque(maxlen=1000))
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
    
    def record_request(self, api_name: str, delay: float = 0.0) -> None:
        """
        Record a request for monitoring.
        
        Args:
            api_name: Name of the API
            delay: Delay applied for rate limiting
        """
        if not self.enable_monitoring:
            return
        
        timestamp = time.time()
        
        # Initialize stats if needed
        if api_name not in self.stats:
            self.stats[api_name] = RateLimitStats(api_name=api_name)
        
        # Update stats
        stats = self.stats[api_name]
        stats.total_requests += 1
        
        # Update average delay (exponential moving average)
        alpha = 0.1  # Smoothing factor
        stats.average_delay = (1 - alpha) * stats.average_delay + alpha * delay
        
        # Record request timestamp
        self.request_history[api_name].append(timestamp)
        
        # Calculate current rate (requests per second over last minute)
        cutoff_time = timestamp - 60.0
        recent_requests = [t for t in self.request_history[api_name] if t > cutoff_time]
        stats.current_rate = len(recent_requests) / 60.0
    
    def record_violation(self, violation: RateLimitViolation) -> None:
        """
        Record a rate limit violation.
        
        Args:
            violation: Rate limit violation information
        """
        if not self.enable_monitoring:
            return
        
        self.violations.append(violation)
        
        # Update stats
        if violation.api_name in self.stats:
            stats = self.stats[violation.api_name]
            stats.total_violations += 1
            stats.last_violation = violation
            stats.violations_by_type[violation.violation_type] += 1
        
        # Log violation
        self.logger.warning(
            "Rate limit violation for %s API",
            violation.api_name,
            extra={
                "api_name": violation.api_name,
                "violation_type": violation.violation_type.value,
                "current_rate": violation.current_rate,
                "limit_rate": violation.limit_rate,
                "delay_applied": violation.delay_applied,
                "requests_in_window": violation.requests_in_window
            }
        )
    
    def get_stats(self, api_name: Optional[str] = None) -> Dict[str, RateLimitStats]:
        """
        Get rate limiting statistics.
        
        Args:
            api_name: Specific API name, or None for all APIs
            
        Returns:
            Dictionary of statistics by API name
        """
        if api_name:
            return {api_name: self.stats.get(api_name, RateLimitStats(api_name=api_name))}
        return dict(self.stats)
    
    def get_recent_violations(self, minutes: int = 60) -> List[RateLimitViolation]:
        """
        Get recent rate limit violations.
        
        Args:
            minutes: Number of minutes to look back
            
        Returns:
            List of recent violations
        """
        cutoff_time = time.time() - (minutes * 60)
        return [v for v in self.violations if v.timestamp > cutoff_time]
    
    def generate_report(self) -> Dict[str, Any]:
        """
        Generate comprehensive rate limiting report.
        
        Returns:
            Dictionary with rate limiting report data
        """
        report = {
            "timestamp": time.time(),
            "apis": {},
            "summary": {
                "total_apis": len(self.stats),
                "total_requests": sum(s.total_requests for s in self.stats.values()),
                "total_violations": sum(s.total_violations for s in self.stats.values()),
                "recent_violations": len(self.get_recent_violations(60))
            }
        }
        
        for api_name, stats in self.stats.items():
            report["apis"][api_name] = {
                "total_requests": stats.total_requests,
                "total_violations": stats.total_violations,
                "current_rate": stats.current_rate,
                "average_delay": stats.average_delay,
                "violations_by_type": dict(stats.violations_by_type),
                "last_violation": {
                    "type": stats.last_violation.violation_type.value,
                    "timestamp": stats.last_violation.timestamp,
                    "delay": stats.last_violation.delay_applied
                } if stats.last_violation else None
            }
        
        return report


class APIRateLimiter:
    """Rate limiter for a specific API with comprehensive monitoring."""
    
    def __init__(self, api_name: str, config: RateLimitConfig, monitor: Optional[RateLimitMonitor] = None):
        """
        Initialize API rate limiter.
        
        Args:
            api_name: Name of the API
            config: Rate limiting configuration
            monitor: Optional monitor for statistics collection
        """
        self.api_name = api_name
        self.config = config
        self.monitor = monitor
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        
        # Token bucket for rate limiting
        self.token_bucket = TokenBucket(config.requests_per_second)
        
        # Burst detection
        self.burst_requests: deque = deque(maxlen=config.burst_limit * 2)
        
        # Violation tracking for exponential backoff
        self.consecutive_violations = 0
        self.last_violation_time = 0.0
        
        # Statistics
        self.total_requests = 0
        self.total_delays = 0.0
    
    async def acquire(self, tokens: int = 1) -> float:
        """
        Acquire permission to make API request with rate limiting.
        
        Args:
            tokens: Number of tokens to acquire (usually 1 per request)
            
        Returns:
            Total delay applied in seconds
        """
        now = time.time()
        total_delay = 0.0
        
        # Check for burst limit violation
        burst_delay = await self._check_burst_limit(now)
        if burst_delay > 0:
            total_delay += burst_delay
            await asyncio.sleep(burst_delay)
        
        # Apply token bucket rate limiting
        token_delay = await self.token_bucket.acquire(tokens)
        if token_delay > 0:
            total_delay += token_delay
            await asyncio.sleep(token_delay)
        
        # Check for soft limit warning
        await self._check_soft_limit()
        
        # Record request
        self.total_requests += 1
        self.total_delays += total_delay
        self.burst_requests.append(now)
        
        if self.monitor:
            self.monitor.record_request(self.api_name, total_delay)
        
        # Log if significant delay was applied
        if total_delay > 0.1:  # Log delays over 100ms
            self.logger.debug(
                "Rate limiting delay applied for %s API",
                self.api_name,
                extra={
                    "api_name": self.api_name,
                    "delay_seconds": total_delay,
                    "tokens_requested": tokens,
                    "current_rate": self._calculate_current_rate()
                }
            )
        
        return total_delay
    
    async def _check_burst_limit(self, now: float) -> float:
        """
        Check for burst limit violations and apply exponential backoff.
        
        Args:
            now: Current timestamp
            
        Returns:
            Delay to apply in seconds
        """
        # Count requests in burst window
        window_start = now - self.config.burst_window_seconds
        requests_in_window = sum(1 for t in self.burst_requests if t > window_start)
        
        if requests_in_window >= self.config.burst_limit:
            # Burst limit exceeded
            self.consecutive_violations += 1
            delay = self._calculate_violation_delay()
            
            violation = RateLimitViolation(
                violation_type=RateLimitViolationType.BURST_LIMIT,
                timestamp=now,
                api_name=self.api_name,
                current_rate=self._calculate_current_rate(),
                limit_rate=self.config.requests_per_second,
                delay_applied=delay,
                requests_in_window=requests_in_window
            )
            
            if self.monitor:
                self.monitor.record_violation(violation)
            
            self.last_violation_time = now
            return delay
        
        # Reset consecutive violations if enough time has passed
        if now - self.last_violation_time > self.config.burst_window_seconds:
            self.consecutive_violations = 0
        
        return 0.0
    
    async def _check_soft_limit(self) -> None:
        """Check for soft limit threshold and log warnings."""
        current_rate = self._calculate_current_rate()
        soft_limit = self.config.requests_per_second * self.config.soft_limit_threshold
        
        if current_rate > soft_limit:
            violation = RateLimitViolation(
                violation_type=RateLimitViolationType.SOFT_LIMIT,
                timestamp=time.time(),
                api_name=self.api_name,
                current_rate=current_rate,
                limit_rate=self.config.requests_per_second,
                delay_applied=0.0
            )
            
            if self.monitor:
                self.monitor.record_violation(violation)
    
    def _calculate_violation_delay(self) -> float:
        """Calculate delay for rate limit violation using exponential backoff."""
        delay = self.config.violation_initial_delay * (
            self.config.violation_backoff_multiplier ** (self.consecutive_violations - 1)
        )
        return min(delay, self.config.violation_max_delay)
    
    def _calculate_current_rate(self) -> float:
        """Calculate current request rate (requests per second)."""
        if not self.burst_requests:
            return 0.0
        
        now = time.time()
        window_start = now - 60.0  # Calculate rate over last minute
        recent_requests = [t for t in self.burst_requests if t > window_start]
        return len(recent_requests) / 60.0
    
    def get_stats(self) -> Dict[str, Any]:
        """Get statistics for this rate limiter."""
        return {
            "api_name": self.api_name,
            "total_requests": self.total_requests,
            "average_delay": self.total_delays / max(self.total_requests, 1),
            "current_rate": self._calculate_current_rate(),
            "consecutive_violations": self.consecutive_violations,
            "current_tokens": self.token_bucket.get_current_tokens(),
            "config": {
                "requests_per_second": self.config.requests_per_second,
                "burst_limit": self.config.burst_limit,
                "burst_window_seconds": self.config.burst_window_seconds
            }
        }


class RateLimitManager:
    """Central manager for all API rate limiters."""
    
    def __init__(self, enable_monitoring: bool = True, enable_reporting: bool = True):
        """
        Initialize rate limit manager.
        
        Args:
            enable_monitoring: Whether to enable monitoring
            enable_reporting: Whether to enable reporting
        """
        self.limiters: Dict[str, APIRateLimiter] = {}
        self.monitor = RateLimitMonitor(enable_monitoring) if enable_monitoring else None
        self.enable_reporting = enable_reporting
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
    
    def create_limiter(self, api_name: str, config: RateLimitConfig) -> APIRateLimiter:
        """
        Create or get rate limiter for an API.
        
        Args:
            api_name: Name of the API
            config: Rate limiting configuration
            
        Returns:
            APIRateLimiter instance
        """
        if api_name not in self.limiters:
            self.limiters[api_name] = APIRateLimiter(api_name, config, self.monitor)
            self.logger.info(
                "Created rate limiter for %s API",
                api_name,
                extra={
                    "api_name": api_name,
                    "requests_per_second": config.requests_per_second,
                    "burst_limit": config.burst_limit
                }
            )
        
        return self.limiters[api_name]
    
    def get_limiter(self, api_name: str) -> Optional[APIRateLimiter]:
        """
        Get existing rate limiter for an API.
        
        Args:
            api_name: Name of the API
            
        Returns:
            APIRateLimiter instance or None if not found
        """
        return self.limiters.get(api_name)
    
    async def acquire(self, api_name: str, tokens: int = 1) -> float:
        """
        Acquire permission for API request.
        
        Args:
            api_name: Name of the API
            tokens: Number of tokens to acquire
            
        Returns:
            Delay applied in seconds
            
        Raises:
            ValueError: If no rate limiter exists for the API
        """
        limiter = self.limiters.get(api_name)
        if not limiter:
            raise ValueError(f"No rate limiter configured for API: {api_name}")
        
        return await limiter.acquire(tokens)
    
    def get_all_stats(self) -> Dict[str, Any]:
        """Get statistics for all rate limiters."""
        stats = {}
        for api_name, limiter in self.limiters.items():
            stats[api_name] = limiter.get_stats()
        
        if self.monitor:
            stats["monitor"] = self.monitor.get_stats()
        
        return stats
    
    def generate_report(self) -> Dict[str, Any]:
        """Generate comprehensive rate limiting report."""
        if not self.enable_reporting or not self.monitor:
            return {"error": "Reporting not enabled"}
        
        return self.monitor.generate_report()
    
    def get_recent_violations(self, minutes: int = 60) -> List[RateLimitViolation]:
        """Get recent rate limit violations across all APIs."""
        if not self.monitor:
            return []
        
        return self.monitor.get_recent_violations(minutes)


# Global rate limit manager instance
_rate_limit_manager: Optional[RateLimitManager] = None


def get_rate_limit_manager() -> RateLimitManager:
    """Get the global rate limit manager instance."""
    global _rate_limit_manager
    if _rate_limit_manager is None:
        _rate_limit_manager = RateLimitManager()
    return _rate_limit_manager


def set_rate_limit_manager(manager: RateLimitManager) -> None:
    """Set the global rate limit manager instance."""
    global _rate_limit_manager
    _rate_limit_manager = manager


# Convenience functions
async def acquire_rate_limit(api_name: str, tokens: int = 1) -> float:
    """
    Acquire rate limit permission for an API.
    
    Args:
        api_name: Name of the API
        tokens: Number of tokens to acquire
        
    Returns:
        Delay applied in seconds
    """
    manager = get_rate_limit_manager()
    return await manager.acquire(api_name, tokens)


def create_api_rate_limiter(api_name: str, requests_per_second: float, **kwargs) -> APIRateLimiter:
    """
    Create rate limiter for an API with default configuration.
    
    Args:
        api_name: Name of the API
        requests_per_second: Rate limit in requests per second
        **kwargs: Additional configuration parameters
        
    Returns:
        APIRateLimiter instance
    """
    config = RateLimitConfig(requests_per_second=requests_per_second, **kwargs)
    manager = get_rate_limit_manager()
    return manager.create_limiter(api_name, config)