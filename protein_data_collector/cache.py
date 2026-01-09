"""
Response caching and performance optimization system.

This module provides configurable response caching with TTL, performance monitoring,
and metrics collection for the protein data collector system.
"""

import asyncio
import hashlib
import json
import logging
import time
from dataclasses import dataclass, field
from typing import Any, Dict, Optional, List, Tuple, Union, Callable
from enum import Enum
import pickle
from collections import defaultdict, OrderedDict
from threading import RLock


class CacheStrategy(Enum):
    """Cache eviction strategies."""
    LRU = "lru"  # Least Recently Used
    LFU = "lfu"  # Least Frequently Used
    TTL = "ttl"  # Time To Live only
    FIFO = "fifo"  # First In, First Out


@dataclass
class CacheConfig:
    """Configuration for response caching."""
    enabled: bool = True
    default_ttl_seconds: int = 3600  # 1 hour default
    max_cache_size: int = 10000  # Maximum number of cached items
    max_memory_mb: int = 500  # Maximum memory usage in MB
    strategy: CacheStrategy = CacheStrategy.LRU
    
    # API-specific TTL settings
    interpro_ttl_seconds: int = 7200  # 2 hours for InterPro
    uniprot_ttl_seconds: int = 3600   # 1 hour for UniProt
    mcp_ttl_seconds: int = 1800       # 30 minutes for MCP
    
    # Performance settings
    cleanup_interval_seconds: int = 300  # 5 minutes
    enable_compression: bool = True
    enable_metrics: bool = True


@dataclass
class CacheEntry:
    """Individual cache entry with metadata."""
    key: str
    value: Any
    created_at: float
    last_accessed: float
    access_count: int
    ttl_seconds: int
    size_bytes: int
    compressed: bool = False
    
    def is_expired(self) -> bool:
        """Check if cache entry has expired."""
        return time.time() - self.created_at > self.ttl_seconds
    
    def touch(self) -> None:
        """Update last accessed time and increment access count."""
        self.last_accessed = time.time()
        self.access_count += 1


@dataclass
class CacheMetrics:
    """Cache performance metrics."""
    hits: int = 0
    misses: int = 0
    evictions: int = 0
    expired_removals: int = 0
    total_requests: int = 0
    total_size_bytes: int = 0
    average_response_time_ms: float = 0.0
    
    @property
    def hit_rate(self) -> float:
        """Calculate cache hit rate."""
        if self.total_requests == 0:
            return 0.0
        return self.hits / self.total_requests
    
    @property
    def miss_rate(self) -> float:
        """Calculate cache miss rate."""
        return 1.0 - self.hit_rate


class ResponseCache:
    """High-performance response cache with configurable eviction strategies."""
    
    def __init__(self, config: CacheConfig):
        """
        Initialize response cache.
        
        Args:
            config: Cache configuration
        """
        self.config = config
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        
        # Thread-safe cache storage
        self._cache: Dict[str, CacheEntry] = {}
        self._access_order: OrderedDict = OrderedDict()  # For LRU
        self._frequency_counter: Dict[str, int] = defaultdict(int)  # For LFU
        self._lock = RLock()
        
        # Metrics
        self.metrics = CacheMetrics()
        self._response_times: List[float] = []
        
        # Background cleanup task
        self._cleanup_task: Optional[asyncio.Task] = None
        if config.enabled:
            self._start_cleanup_task()
    
    def _start_cleanup_task(self) -> None:
        """Start background cleanup task."""
        async def cleanup_loop():
            while True:
                try:
                    await asyncio.sleep(self.config.cleanup_interval_seconds)
                    self._cleanup_expired()
                except asyncio.CancelledError:
                    break
                except Exception as e:
                    self.logger.error(f"Cache cleanup error: {e}")
        
        self._cleanup_task = asyncio.create_task(cleanup_loop())
    
    def _generate_key(self, api_name: str, endpoint: str, params: Optional[Dict[str, Any]] = None) -> str:
        """
        Generate cache key from API call parameters.
        
        Args:
            api_name: Name of the API
            endpoint: API endpoint
            params: Request parameters
            
        Returns:
            Cache key string
        """
        # Create deterministic key from parameters
        key_data = {
            "api": api_name,
            "endpoint": endpoint,
            "params": params or {}
        }
        
        # Sort parameters for consistent key generation
        key_json = json.dumps(key_data, sort_keys=True, separators=(',', ':'))
        key_hash = hashlib.sha256(key_json.encode()).hexdigest()[:16]
        
        return f"{api_name}:{endpoint}:{key_hash}"
    
    def _serialize_value(self, value: Any) -> Tuple[bytes, bool]:
        """
        Serialize value for storage with optional compression.
        
        Args:
            value: Value to serialize
            
        Returns:
            Tuple of (serialized_data, is_compressed)
        """
        try:
            # Serialize using pickle
            data = pickle.dumps(value)
            
            # Compress if enabled and data is large enough
            if self.config.enable_compression and len(data) > 1024:  # 1KB threshold
                import gzip
                compressed_data = gzip.compress(data)
                if len(compressed_data) < len(data) * 0.9:  # Only use if 10%+ savings
                    return compressed_data, True
            
            return data, False
            
        except Exception as e:
            self.logger.warning(f"Failed to serialize cache value: {e}")
            return b"", False
    
    def _deserialize_value(self, data: bytes, compressed: bool) -> Any:
        """
        Deserialize value from storage.
        
        Args:
            data: Serialized data
            compressed: Whether data is compressed
            
        Returns:
            Deserialized value
        """
        try:
            if compressed:
                import gzip
                data = gzip.decompress(data)
            
            return pickle.loads(data)
            
        except Exception as e:
            self.logger.warning(f"Failed to deserialize cache value: {e}")
            return None
    
    def _get_ttl_for_api(self, api_name: str) -> int:
        """
        Get TTL for specific API.
        
        Args:
            api_name: Name of the API
            
        Returns:
            TTL in seconds
        """
        api_ttl_map = {
            "InterPro": self.config.interpro_ttl_seconds,
            "UniProt": self.config.uniprot_ttl_seconds,
            "UniProt_REST": self.config.uniprot_ttl_seconds,
            "UniProt_MCP": self.config.mcp_ttl_seconds,
        }
        
        return api_ttl_map.get(api_name, self.config.default_ttl_seconds)
    
    def _evict_if_needed(self) -> None:
        """Evict entries if cache is over limits."""
        # Check size limits
        if len(self._cache) <= self.config.max_cache_size:
            total_size_mb = self.metrics.total_size_bytes / (1024 * 1024)
            if total_size_mb <= self.config.max_memory_mb:
                return
        
        # Determine how many entries to evict
        target_size = int(self.config.max_cache_size * 0.8)  # Evict to 80% capacity
        entries_to_evict = len(self._cache) - target_size
        
        if entries_to_evict <= 0:
            return
        
        # Select entries to evict based on strategy
        if self.config.strategy == CacheStrategy.LRU:
            # Evict least recently used
            keys_to_evict = list(self._access_order.keys())[:entries_to_evict]
        elif self.config.strategy == CacheStrategy.LFU:
            # Evict least frequently used
            sorted_by_frequency = sorted(
                self._cache.keys(),
                key=lambda k: self._cache[k].access_count
            )
            keys_to_evict = sorted_by_frequency[:entries_to_evict]
        elif self.config.strategy == CacheStrategy.FIFO:
            # Evict oldest entries
            sorted_by_age = sorted(
                self._cache.keys(),
                key=lambda k: self._cache[k].created_at
            )
            keys_to_evict = sorted_by_age[:entries_to_evict]
        else:  # TTL strategy - evict expired first, then oldest
            expired_keys = [k for k, v in self._cache.items() if v.is_expired()]
            if len(expired_keys) >= entries_to_evict:
                keys_to_evict = expired_keys[:entries_to_evict]
            else:
                # Add oldest non-expired entries
                non_expired = [k for k in self._cache.keys() if k not in expired_keys]
                sorted_non_expired = sorted(
                    non_expired,
                    key=lambda k: self._cache[k].created_at
                )
                remaining_to_evict = entries_to_evict - len(expired_keys)
                keys_to_evict = expired_keys + sorted_non_expired[:remaining_to_evict]
        
        # Evict selected entries
        for key in keys_to_evict:
            self._remove_entry(key)
            self.metrics.evictions += 1
        
        self.logger.debug(
            f"Evicted {len(keys_to_evict)} cache entries using {self.config.strategy.value} strategy"
        )
    
    def _remove_entry(self, key: str) -> None:
        """Remove entry from all cache structures."""
        if key in self._cache:
            entry = self._cache[key]
            self.metrics.total_size_bytes -= entry.size_bytes
            del self._cache[key]
        
        if key in self._access_order:
            del self._access_order[key]
        
        if key in self._frequency_counter:
            del self._frequency_counter[key]
    
    def _cleanup_expired(self) -> None:
        """Remove expired entries from cache."""
        with self._lock:
            expired_keys = [k for k, v in self._cache.items() if v.is_expired()]
            
            for key in expired_keys:
                self._remove_entry(key)
                self.metrics.expired_removals += 1
            
            if expired_keys:
                self.logger.debug(f"Cleaned up {len(expired_keys)} expired cache entries")
    
    def get(
        self,
        api_name: str,
        endpoint: str,
        params: Optional[Dict[str, Any]] = None
    ) -> Optional[Any]:
        """
        Get cached response.
        
        Args:
            api_name: Name of the API
            endpoint: API endpoint
            params: Request parameters
            
        Returns:
            Cached response or None if not found/expired
        """
        if not self.config.enabled:
            return None
        
        key = self._generate_key(api_name, endpoint, params)
        
        with self._lock:
            self.metrics.total_requests += 1
            
            if key not in self._cache:
                self.metrics.misses += 1
                return None
            
            entry = self._cache[key]
            
            # Check if expired
            if entry.is_expired():
                self._remove_entry(key)
                self.metrics.expired_removals += 1
                self.metrics.misses += 1
                return None
            
            # Update access tracking
            entry.touch()
            
            # Update LRU order
            if self.config.strategy == CacheStrategy.LRU:
                self._access_order.move_to_end(key)
            
            # Update frequency counter
            if self.config.strategy == CacheStrategy.LFU:
                self._frequency_counter[key] += 1
            
            self.metrics.hits += 1
            
            # Deserialize value
            value = self._deserialize_value(
                entry.value if isinstance(entry.value, bytes) else pickle.dumps(entry.value),
                entry.compressed
            )
            
            self.logger.debug(
                f"Cache hit for {api_name}:{endpoint}",
                extra={
                    "api_name": api_name,
                    "endpoint": endpoint,
                    "cache_key": key,
                    "age_seconds": time.time() - entry.created_at
                }
            )
            
            return value
    
    def put(
        self,
        api_name: str,
        endpoint: str,
        params: Optional[Dict[str, Any]],
        value: Any,
        ttl_seconds: Optional[int] = None
    ) -> None:
        """
        Store response in cache.
        
        Args:
            api_name: Name of the API
            endpoint: API endpoint
            params: Request parameters
            value: Response value to cache
            ttl_seconds: Custom TTL, uses API default if None
        """
        if not self.config.enabled:
            return
        
        key = self._generate_key(api_name, endpoint, params)
        ttl = ttl_seconds or self._get_ttl_for_api(api_name)
        
        # Serialize value
        serialized_data, compressed = self._serialize_value(value)
        if not serialized_data:
            return  # Failed to serialize
        
        with self._lock:
            # Remove existing entry if present
            if key in self._cache:
                self._remove_entry(key)
            
            # Create new entry
            entry = CacheEntry(
                key=key,
                value=serialized_data,
                created_at=time.time(),
                last_accessed=time.time(),
                access_count=1,
                ttl_seconds=ttl,
                size_bytes=len(serialized_data),
                compressed=compressed
            )
            
            # Store entry
            self._cache[key] = entry
            self.metrics.total_size_bytes += entry.size_bytes
            
            # Update tracking structures
            if self.config.strategy == CacheStrategy.LRU:
                self._access_order[key] = True
            
            if self.config.strategy == CacheStrategy.LFU:
                self._frequency_counter[key] = 1
            
            # Evict if needed
            self._evict_if_needed()
            
            self.logger.debug(
                f"Cached response for {api_name}:{endpoint}",
                extra={
                    "api_name": api_name,
                    "endpoint": endpoint,
                    "cache_key": key,
                    "ttl_seconds": ttl,
                    "size_bytes": entry.size_bytes,
                    "compressed": compressed
                }
            )
    
    def invalidate(
        self,
        api_name: Optional[str] = None,
        endpoint: Optional[str] = None,
        params: Optional[Dict[str, Any]] = None
    ) -> int:
        """
        Invalidate cache entries.
        
        Args:
            api_name: API name to invalidate (None for all)
            endpoint: Endpoint to invalidate (None for all in API)
            params: Specific parameters to invalidate (None for all in endpoint)
            
        Returns:
            Number of entries invalidated
        """
        if not self.config.enabled:
            return 0
        
        with self._lock:
            if api_name is None:
                # Invalidate all
                count = len(self._cache)
                self._cache.clear()
                self._access_order.clear()
                self._frequency_counter.clear()
                self.metrics.total_size_bytes = 0
                return count
            
            if endpoint is None:
                # Invalidate all for API
                keys_to_remove = [k for k in self._cache.keys() if k.startswith(f"{api_name}:")]
            elif params is None:
                # Invalidate all for API and endpoint
                prefix = f"{api_name}:{endpoint}:"
                keys_to_remove = [k for k in self._cache.keys() if k.startswith(prefix)]
            else:
                # Invalidate specific entry
                key = self._generate_key(api_name, endpoint, params)
                keys_to_remove = [key] if key in self._cache else []
            
            for key in keys_to_remove:
                self._remove_entry(key)
            
            return len(keys_to_remove)
    
    def get_metrics(self) -> Dict[str, Any]:
        """Get cache performance metrics."""
        with self._lock:
            # Update average response time
            if self._response_times:
                self.metrics.average_response_time_ms = sum(self._response_times) / len(self._response_times)
                # Keep only recent response times
                if len(self._response_times) > 1000:
                    self._response_times = self._response_times[-500:]
            
            return {
                "enabled": self.config.enabled,
                "total_entries": len(self._cache),
                "total_size_mb": self.metrics.total_size_bytes / (1024 * 1024),
                "hit_rate": self.metrics.hit_rate,
                "miss_rate": self.metrics.miss_rate,
                "hits": self.metrics.hits,
                "misses": self.metrics.misses,
                "evictions": self.metrics.evictions,
                "expired_removals": self.metrics.expired_removals,
                "total_requests": self.metrics.total_requests,
                "average_response_time_ms": self.metrics.average_response_time_ms,
                "config": {
                    "max_cache_size": self.config.max_cache_size,
                    "max_memory_mb": self.config.max_memory_mb,
                    "strategy": self.config.strategy.value,
                    "default_ttl_seconds": self.config.default_ttl_seconds
                }
            }
    
    def record_response_time(self, response_time_ms: float) -> None:
        """Record API response time for metrics."""
        if self.config.enable_metrics:
            self._response_times.append(response_time_ms)
    
    async def close(self) -> None:
        """Close cache and cleanup resources."""
        if self._cleanup_task:
            self._cleanup_task.cancel()
            try:
                await self._cleanup_task
            except asyncio.CancelledError:
                pass
        
        with self._lock:
            self._cache.clear()
            self._access_order.clear()
            self._frequency_counter.clear()


class CachedAPIClient:
    """Wrapper for API clients that adds caching functionality."""
    
    def __init__(self, cache: ResponseCache):
        """
        Initialize cached API client wrapper.
        
        Args:
            cache: Response cache instance
        """
        self.cache = cache
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
    
    async def cached_request(
        self,
        api_name: str,
        endpoint: str,
        params: Optional[Dict[str, Any]],
        request_func: Callable,
        ttl_seconds: Optional[int] = None
    ) -> Any:
        """
        Make cached API request.
        
        Args:
            api_name: Name of the API
            endpoint: API endpoint
            params: Request parameters
            request_func: Function to call if cache miss
            ttl_seconds: Custom TTL for this request
            
        Returns:
            API response (from cache or fresh request)
        """
        # Try cache first
        cached_response = self.cache.get(api_name, endpoint, params)
        if cached_response is not None:
            return cached_response
        
        # Cache miss - make actual request
        start_time = time.time()
        
        try:
            response = await request_func()
            
            # Record response time
            response_time_ms = (time.time() - start_time) * 1000
            self.cache.record_response_time(response_time_ms)
            
            # Cache the response
            self.cache.put(api_name, endpoint, params, response, ttl_seconds)
            
            self.logger.debug(
                f"Fresh API request for {api_name}:{endpoint}",
                extra={
                    "api_name": api_name,
                    "endpoint": endpoint,
                    "response_time_ms": response_time_ms,
                    "cached": True
                }
            )
            
            return response
            
        except Exception as e:
            # Record failed response time
            response_time_ms = (time.time() - start_time) * 1000
            self.cache.record_response_time(response_time_ms)
            
            self.logger.warning(
                f"API request failed for {api_name}:{endpoint}: {e}",
                extra={
                    "api_name": api_name,
                    "endpoint": endpoint,
                    "response_time_ms": response_time_ms,
                    "error": str(e)
                }
            )
            
            raise


# Global cache instance
_global_cache: Optional[ResponseCache] = None


def get_global_cache() -> ResponseCache:
    """Get the global cache instance."""
    global _global_cache
    if _global_cache is None:
        from .config import get_config
        config = get_config()
        
        # Create cache config from system config
        cache_config = CacheConfig(
            enabled=True,
            default_ttl_seconds=config.collection.cache_ttl_hours * 3600,
            interpro_ttl_seconds=config.collection.cache_ttl_hours * 3600 * 2,  # 2x default for InterPro
            uniprot_ttl_seconds=config.collection.cache_ttl_hours * 3600,
            mcp_ttl_seconds=config.collection.cache_ttl_hours * 3600 // 2,  # Half default for MCP
            enable_metrics=True
        )
        
        _global_cache = ResponseCache(cache_config)
    
    return _global_cache


def set_global_cache(cache: ResponseCache) -> None:
    """Set the global cache instance."""
    global _global_cache
    _global_cache = cache


# Convenience functions
def cached_api_call(api_name: str, endpoint: str, params: Optional[Dict[str, Any]] = None):
    """
    Decorator for caching API calls.
    
    Args:
        api_name: Name of the API
        endpoint: API endpoint
        params: Request parameters
        
    Returns:
        Decorator function
    """
    def decorator(func):
        async def wrapper(*args, **kwargs):
            cache = get_global_cache()
            cached_client = CachedAPIClient(cache)
            
            async def request_func():
                return await func(*args, **kwargs)
            
            return await cached_client.cached_request(
                api_name, endpoint, params, request_func
            )
        
        return wrapper
    return decorator


async def get_cache_stats() -> Dict[str, Any]:
    """Get global cache statistics."""
    cache = get_global_cache()
    return cache.get_metrics()


async def invalidate_cache(
    api_name: Optional[str] = None,
    endpoint: Optional[str] = None,
    params: Optional[Dict[str, Any]] = None
) -> int:
    """Invalidate global cache entries."""
    cache = get_global_cache()
    return cache.invalidate(api_name, endpoint, params)