"""
Performance monitoring and metrics collection system.

This module provides comprehensive performance monitoring for the protein data
collector system, including API response times, database query performance,
and system resource usage.
"""

import asyncio
import logging
import psutil
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any, Callable
from collections import defaultdict, deque
from enum import Enum
import statistics


class MetricType(Enum):
    """Types of performance metrics."""
    RESPONSE_TIME = "response_time"
    THROUGHPUT = "throughput"
    ERROR_RATE = "error_rate"
    RESOURCE_USAGE = "resource_usage"
    CACHE_PERFORMANCE = "cache_performance"
    DATABASE_PERFORMANCE = "database_performance"


@dataclass
class PerformanceMetric:
    """Individual performance metric data point."""
    metric_type: MetricType
    name: str
    value: float
    timestamp: float
    tags: Dict[str, str] = field(default_factory=dict)
    unit: str = ""


@dataclass
class SystemResourceMetrics:
    """System resource usage metrics."""
    cpu_percent: float
    memory_percent: float
    memory_used_mb: float
    memory_available_mb: float
    disk_usage_percent: float
    network_bytes_sent: int
    network_bytes_recv: int
    timestamp: float


class PerformanceMonitor:
    """Comprehensive performance monitoring system."""
    
    def __init__(self, enable_monitoring: bool = True, max_metrics: int = 10000):
        """
        Initialize performance monitor.
        
        Args:
            enable_monitoring: Whether to collect metrics
            max_metrics: Maximum number of metrics to keep in memory
        """
        self.enable_monitoring = enable_monitoring
        self.max_metrics = max_metrics
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        
        # Metric storage
        self.metrics: Dict[str, deque] = defaultdict(lambda: deque(maxlen=max_metrics))
        self.counters: Dict[str, int] = defaultdict(int)
        self.gauges: Dict[str, float] = defaultdict(float)
        
        # System monitoring
        self.system_metrics: deque = deque(maxlen=1000)
        self._monitoring_task: Optional[asyncio.Task] = None
        
        # Performance tracking
        self.active_operations: Dict[str, float] = {}  # operation_id -> start_time
        
        if enable_monitoring:
            self._start_system_monitoring()
    
    def _start_system_monitoring(self) -> None:
        """Start background system monitoring task."""
        try:
            # Only start monitoring if there's an event loop running
            loop = asyncio.get_running_loop()
            
            async def monitor_system():
                while True:
                    try:
                        await asyncio.sleep(30)  # Monitor every 30 seconds
                        self._collect_system_metrics()
                    except asyncio.CancelledError:
                        break
                    except Exception as e:
                        self.logger.error(f"System monitoring error: {e}")
            
            self._monitoring_task = asyncio.create_task(monitor_system())
        except RuntimeError:
            # No event loop running, skip monitoring
            self.logger.debug("No event loop running, skipping system monitoring")
            self._monitoring_task = None
    
    def _collect_system_metrics(self) -> None:
        """Collect current system resource metrics."""
        try:
            # CPU and memory
            cpu_percent = psutil.cpu_percent(interval=1)
            memory = psutil.virtual_memory()
            
            # Disk usage for current directory
            disk = psutil.disk_usage('.')
            
            # Network I/O
            network = psutil.net_io_counters()
            
            metrics = SystemResourceMetrics(
                cpu_percent=cpu_percent,
                memory_percent=memory.percent,
                memory_used_mb=memory.used / (1024 * 1024),
                memory_available_mb=memory.available / (1024 * 1024),
                disk_usage_percent=(disk.used / disk.total) * 100,
                network_bytes_sent=network.bytes_sent,
                network_bytes_recv=network.bytes_recv,
                timestamp=time.time()
            )
            
            self.system_metrics.append(metrics)
            
            # Update gauges
            self.gauges["system.cpu_percent"] = cpu_percent
            self.gauges["system.memory_percent"] = memory.percent
            self.gauges["system.memory_used_mb"] = metrics.memory_used_mb
            self.gauges["system.disk_usage_percent"] = metrics.disk_usage_percent
            
        except Exception as e:
            self.logger.warning(f"Failed to collect system metrics: {e}")
    
    def record_metric(
        self,
        metric_type: MetricType,
        name: str,
        value: float,
        tags: Optional[Dict[str, str]] = None,
        unit: str = ""
    ) -> None:
        """
        Record a performance metric.
        
        Args:
            metric_type: Type of metric
            name: Metric name
            value: Metric value
            tags: Optional tags for categorization
            unit: Unit of measurement
        """
        if not self.enable_monitoring:
            return
        
        metric = PerformanceMetric(
            metric_type=metric_type,
            name=name,
            value=value,
            timestamp=time.time(),
            tags=tags or {},
            unit=unit
        )
        
        self.metrics[name].append(metric)
        
        # Update gauges for latest values
        self.gauges[name] = value
        
        self.logger.debug(
            f"Recorded metric: {name} = {value} {unit}",
            extra={
                "metric_type": metric_type.value,
                "metric_name": name,
                "metric_value": value,
                "metric_unit": unit,
                "tags": tags
            }
        )
    
    def increment_counter(self, name: str, value: int = 1, tags: Optional[Dict[str, str]] = None) -> None:
        """
        Increment a counter metric.
        
        Args:
            name: Counter name
            value: Increment value
            tags: Optional tags
        """
        if not self.enable_monitoring:
            return
        
        self.counters[name] += value
        
        # Also record as a metric for time series
        self.record_metric(
            MetricType.THROUGHPUT,
            f"{name}.count",
            self.counters[name],
            tags,
            "count"
        )
    
    def set_gauge(self, name: str, value: float, tags: Optional[Dict[str, str]] = None) -> None:
        """
        Set a gauge metric value.
        
        Args:
            name: Gauge name
            value: Gauge value
            tags: Optional tags
        """
        if not self.enable_monitoring:
            return
        
        self.gauges[name] = value
        
        # Also record as a metric for time series
        self.record_metric(
            MetricType.RESOURCE_USAGE,
            name,
            value,
            tags
        )
    
    def start_operation(self, operation_id: str) -> None:
        """
        Start timing an operation.
        
        Args:
            operation_id: Unique identifier for the operation
        """
        if not self.enable_monitoring:
            return
        
        self.active_operations[operation_id] = time.time()
    
    def end_operation(self, operation_id: str, tags: Optional[Dict[str, str]] = None) -> float:
        """
        End timing an operation and record the duration.
        
        Args:
            operation_id: Unique identifier for the operation
            tags: Optional tags for the metric
            
        Returns:
            Operation duration in seconds
        """
        if not self.enable_monitoring:
            return 0.0
        
        if operation_id not in self.active_operations:
            self.logger.warning(f"Operation {operation_id} was not started")
            return 0.0
        
        start_time = self.active_operations.pop(operation_id)
        duration = time.time() - start_time
        
        # Record the duration metric
        self.record_metric(
            MetricType.RESPONSE_TIME,
            f"operation.{operation_id}.duration",
            duration,
            tags,
            "seconds"
        )
        
        return duration
    
    def record_api_call(
        self,
        api_name: str,
        endpoint: str,
        duration_ms: float,
        success: bool,
        status_code: Optional[int] = None
    ) -> None:
        """
        Record API call performance metrics.
        
        Args:
            api_name: Name of the API
            endpoint: API endpoint
            duration_ms: Response time in milliseconds
            success: Whether the call was successful
            status_code: HTTP status code
        """
        tags = {
            "api": api_name,
            "endpoint": endpoint,
            "success": str(success).lower()
        }
        
        if status_code:
            tags["status_code"] = str(status_code)
        
        # Record response time
        self.record_metric(
            MetricType.RESPONSE_TIME,
            f"api.{api_name}.response_time",
            duration_ms,
            tags,
            "ms"
        )
        
        # Increment call counter
        self.increment_counter(f"api.{api_name}.calls", 1, tags)
        
        # Track error rate
        if not success:
            self.increment_counter(f"api.{api_name}.errors", 1, tags)
    
    def record_database_query(
        self,
        query_type: str,
        table: str,
        duration_ms: float,
        rows_affected: Optional[int] = None
    ) -> None:
        """
        Record database query performance metrics.
        
        Args:
            query_type: Type of query (SELECT, INSERT, UPDATE, DELETE)
            table: Database table name
            duration_ms: Query duration in milliseconds
            rows_affected: Number of rows affected
        """
        tags = {
            "query_type": query_type.upper(),
            "table": table
        }
        
        # Record query time
        self.record_metric(
            MetricType.DATABASE_PERFORMANCE,
            f"database.{table}.query_time",
            duration_ms,
            tags,
            "ms"
        )
        
        # Record rows affected if provided
        if rows_affected is not None:
            self.record_metric(
                MetricType.DATABASE_PERFORMANCE,
                f"database.{table}.rows_affected",
                rows_affected,
                tags,
                "rows"
            )
        
        # Increment query counter
        self.increment_counter(f"database.{table}.queries", 1, tags)
    
    def get_metric_summary(self, metric_name: str, minutes: int = 60) -> Dict[str, Any]:
        """
        Get summary statistics for a metric over a time period.
        
        Args:
            metric_name: Name of the metric
            minutes: Time period in minutes
            
        Returns:
            Dictionary with summary statistics
        """
        if metric_name not in self.metrics:
            return {"error": f"Metric {metric_name} not found"}
        
        cutoff_time = time.time() - (minutes * 60)
        recent_metrics = [
            m for m in self.metrics[metric_name]
            if m.timestamp > cutoff_time
        ]
        
        if not recent_metrics:
            return {"error": f"No recent data for metric {metric_name}"}
        
        values = [m.value for m in recent_metrics]
        
        return {
            "metric_name": metric_name,
            "time_period_minutes": minutes,
            "count": len(values),
            "min": min(values),
            "max": max(values),
            "mean": statistics.mean(values),
            "median": statistics.median(values),
            "std_dev": statistics.stdev(values) if len(values) > 1 else 0.0,
            "percentile_95": statistics.quantiles(values, n=20)[18] if len(values) >= 20 else max(values),
            "percentile_99": statistics.quantiles(values, n=100)[98] if len(values) >= 100 else max(values),
            "unit": recent_metrics[0].unit if recent_metrics else ""
        }
    
    def get_api_performance_report(self, api_name: str, minutes: int = 60) -> Dict[str, Any]:
        """
        Get comprehensive performance report for an API.
        
        Args:
            api_name: Name of the API
            minutes: Time period in minutes
            
        Returns:
            Dictionary with API performance report
        """
        response_time_metric = f"api.{api_name}.response_time"
        calls_metric = f"api.{api_name}.calls.count"
        errors_metric = f"api.{api_name}.errors.count"
        
        response_time_summary = self.get_metric_summary(response_time_metric, minutes)
        
        # Calculate throughput and error rate
        total_calls = self.counters.get(f"api.{api_name}.calls", 0)
        total_errors = self.counters.get(f"api.{api_name}.errors", 0)
        
        error_rate = (total_errors / total_calls * 100) if total_calls > 0 else 0.0
        throughput = total_calls / (minutes * 60) if minutes > 0 else 0.0
        
        return {
            "api_name": api_name,
            "time_period_minutes": minutes,
            "response_time": response_time_summary,
            "total_calls": total_calls,
            "total_errors": total_errors,
            "error_rate_percent": error_rate,
            "throughput_per_second": throughput,
            "current_gauge_values": {
                k: v for k, v in self.gauges.items()
                if k.startswith(f"api.{api_name}")
            }
        }
    
    def get_system_performance_report(self, minutes: int = 60) -> Dict[str, Any]:
        """
        Get system resource performance report.
        
        Args:
            minutes: Time period in minutes
            
        Returns:
            Dictionary with system performance report
        """
        cutoff_time = time.time() - (minutes * 60)
        recent_metrics = [
            m for m in self.system_metrics
            if m.timestamp > cutoff_time
        ]
        
        if not recent_metrics:
            return {"error": "No recent system metrics available"}
        
        # Calculate averages
        avg_cpu = statistics.mean([m.cpu_percent for m in recent_metrics])
        avg_memory = statistics.mean([m.memory_percent for m in recent_metrics])
        avg_disk = statistics.mean([m.disk_usage_percent for m in recent_metrics])
        
        # Get current values
        current = recent_metrics[-1] if recent_metrics else None
        
        return {
            "time_period_minutes": minutes,
            "sample_count": len(recent_metrics),
            "averages": {
                "cpu_percent": avg_cpu,
                "memory_percent": avg_memory,
                "disk_usage_percent": avg_disk
            },
            "current": {
                "cpu_percent": current.cpu_percent if current else 0,
                "memory_percent": current.memory_percent if current else 0,
                "memory_used_mb": current.memory_used_mb if current else 0,
                "memory_available_mb": current.memory_available_mb if current else 0,
                "disk_usage_percent": current.disk_usage_percent if current else 0,
                "timestamp": current.timestamp if current else 0
            } if current else {}
        }
    
    def get_comprehensive_report(self, minutes: int = 60) -> Dict[str, Any]:
        """
        Get comprehensive performance report for all monitored components.
        
        Args:
            minutes: Time period in minutes
            
        Returns:
            Dictionary with comprehensive performance report
        """
        # Get API reports for known APIs
        api_reports = {}
        known_apis = ["InterPro", "UniProt_REST"]
        
        for api_name in known_apis:
            if any(k.startswith(f"api.{api_name}") for k in self.counters.keys()):
                api_reports[api_name] = self.get_api_performance_report(api_name, minutes)
        
        # Get system report
        system_report = self.get_system_performance_report(minutes)
        
        # Get top metrics by activity
        active_metrics = {}
        for metric_name, metric_deque in self.metrics.items():
            if len(metric_deque) > 0:
                recent_count = len([
                    m for m in metric_deque
                    if m.timestamp > time.time() - (minutes * 60)
                ])
                if recent_count > 0:
                    active_metrics[metric_name] = recent_count
        
        # Sort by activity
        top_metrics = sorted(active_metrics.items(), key=lambda x: x[1], reverse=True)[:10]
        
        return {
            "report_timestamp": time.time(),
            "time_period_minutes": minutes,
            "monitoring_enabled": self.enable_monitoring,
            "api_performance": api_reports,
            "system_performance": system_report,
            "top_active_metrics": dict(top_metrics),
            "total_metrics_collected": sum(len(deque) for deque in self.metrics.values()),
            "active_operations": len(self.active_operations),
            "counters_summary": dict(self.counters),
            "gauges_summary": dict(self.gauges)
        }
    
    async def close(self) -> None:
        """Close performance monitor and cleanup resources."""
        if self._monitoring_task:
            self._monitoring_task.cancel()
            try:
                await self._monitoring_task
            except asyncio.CancelledError:
                pass


# Global performance monitor instance
_global_monitor: Optional[PerformanceMonitor] = None


def get_performance_monitor() -> PerformanceMonitor:
    """Get the global performance monitor instance."""
    global _global_monitor
    if _global_monitor is None:
        from .config import get_config
        config = get_config()
        _global_monitor = PerformanceMonitor(
            enable_monitoring=config.collection.enable_performance_monitoring
        )
    return _global_monitor


def set_performance_monitor(monitor: PerformanceMonitor) -> None:
    """Set the global performance monitor instance."""
    global _global_monitor
    _global_monitor = monitor


# Convenience functions and decorators
def record_api_performance(api_name: str, endpoint: str):
    """
    Decorator for recording API call performance.
    
    Args:
        api_name: Name of the API
        endpoint: API endpoint
        
    Returns:
        Decorator function
    """
    def decorator(func):
        async def wrapper(*args, **kwargs):
            monitor = get_performance_monitor()
            start_time = time.time()
            success = False
            status_code = None
            
            try:
                result = await func(*args, **kwargs)
                success = True
                
                # Try to extract status code from result if it's a response object
                if hasattr(result, 'status_code'):
                    status_code = result.status_code
                
                return result
                
            except Exception as e:
                # Try to extract status code from exception
                if hasattr(e, 'status_code'):
                    status_code = e.status_code
                raise
                
            finally:
                duration_ms = (time.time() - start_time) * 1000
                monitor.record_api_call(api_name, endpoint, duration_ms, success, status_code)
        
        return wrapper
    return decorator


def record_database_performance(query_type: str, table: str):
    """
    Decorator for recording database query performance.
    
    Args:
        query_type: Type of query
        table: Database table name
        
    Returns:
        Decorator function
    """
    def decorator(func):
        async def wrapper(*args, **kwargs):
            monitor = get_performance_monitor()
            start_time = time.time()
            
            try:
                result = await func(*args, **kwargs)
                
                # Try to extract rows affected from result
                rows_affected = None
                if hasattr(result, 'rowcount'):
                    rows_affected = result.rowcount
                elif isinstance(result, (list, tuple)):
                    rows_affected = len(result)
                
                return result
                
            finally:
                duration_ms = (time.time() - start_time) * 1000
                monitor.record_database_query(query_type, table, duration_ms, rows_affected)
        
        return wrapper
    return decorator


class PerformanceTimer:
    """Context manager for timing operations."""
    
    def __init__(self, operation_name: str, tags: Optional[Dict[str, str]] = None):
        """
        Initialize performance timer.
        
        Args:
            operation_name: Name of the operation
            tags: Optional tags for the metric
        """
        self.operation_name = operation_name
        self.tags = tags
        self.monitor = get_performance_monitor()
        self.start_time = None
    
    def __enter__(self):
        """Start timing."""
        self.start_time = time.time()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """End timing and record metric."""
        if self.start_time:
            duration = time.time() - self.start_time
            self.monitor.record_metric(
                MetricType.RESPONSE_TIME,
                f"operation.{self.operation_name}.duration",
                duration,
                self.tags,
                "seconds"
            )


# Convenience functions
async def get_performance_report(minutes: int = 60) -> Dict[str, Any]:
    """Get comprehensive performance report."""
    monitor = get_performance_monitor()
    return monitor.get_comprehensive_report(minutes)


async def get_api_report(api_name: str, minutes: int = 60) -> Dict[str, Any]:
    """Get performance report for specific API."""
    monitor = get_performance_monitor()
    return monitor.get_api_performance_report(api_name, minutes)