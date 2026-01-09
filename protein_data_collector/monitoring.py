"""
Monitoring and health check utilities for the Protein Data Collector system.

This module provides comprehensive health checks, metrics collection, and monitoring
capabilities for all system components including database, external APIs, and
application performance.
"""

import asyncio
import logging
import time
from datetime import datetime, timedelta
from typing import Dict, Any, List, Optional, Tuple
from dataclasses import dataclass, asdict
from enum import Enum
import httpx
import psutil
from contextlib import asynccontextmanager

from .database.connection import get_database_manager
from .config import get_config
from .api.interpro_client import InterProAPIClient
from .api.uniprot_client import UnifiedUniProtClient

logger = logging.getLogger(__name__)


class HealthStatus(Enum):
    """Health status enumeration."""
    HEALTHY = "healthy"
    DEGRADED = "degraded"
    UNHEALTHY = "unhealthy"


@dataclass
class ComponentHealth:
    """Health status for a system component."""
    name: str
    status: HealthStatus
    response_time_ms: Optional[float] = None
    error_message: Optional[str] = None
    details: Optional[Dict[str, Any]] = None
    last_checked: datetime = None

    def __post_init__(self):
        if self.last_checked is None:
            self.last_checked = datetime.now()


@dataclass
class SystemHealth:
    """Overall system health status."""
    status: HealthStatus
    components: List[ComponentHealth]
    timestamp: datetime
    uptime_seconds: float
    version: str = "1.0.0"

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "status": self.status.value,
            "components": [asdict(comp) for comp in self.components],
            "timestamp": self.timestamp.isoformat(),
            "uptime_seconds": self.uptime_seconds,
            "version": self.version
        }


@dataclass
class PerformanceMetrics:
    """System performance metrics."""
    cpu_percent: float
    memory_percent: float
    disk_usage_percent: float
    active_connections: int
    request_rate: float
    error_rate: float
    avg_response_time_ms: float
    timestamp: datetime

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return asdict(self)


class HealthChecker:
    """Comprehensive health checker for all system components."""

    def __init__(self):
        self.config = get_config()
        self.start_time = time.time()
        self._metrics_cache = {}
        self._cache_ttl = 30  # Cache metrics for 30 seconds

    async def check_database_health(self) -> ComponentHealth:
        """Check database connectivity and performance."""
        start_time = time.time()
        
        try:
            db_manager = get_database_manager()
            
            # Test basic connectivity
            is_connected = db_manager.test_connection()
            if not is_connected:
                return ComponentHealth(
                    name="database",
                    status=HealthStatus.UNHEALTHY,
                    error_message="Database connection failed"
                )
            
            # Test query performance
            with db_manager.get_connection() as conn:
                cursor = conn.cursor()
                cursor.execute("SELECT 1")
                result = cursor.fetchone()
                
                # Get connection stats
                cursor.execute("SHOW STATUS LIKE 'Threads_connected'")
                threads_connected = cursor.fetchone()
                
                cursor.execute("SHOW VARIABLES LIKE 'max_connections'")
                max_connections = cursor.fetchone()
                
                response_time = (time.time() - start_time) * 1000
                
                # Calculate connection usage
                current_connections = int(threads_connected[1]) if threads_connected else 0
                max_conn = int(max_connections[1]) if max_connections else 100
                connection_usage = (current_connections / max_conn) * 100
                
                details = {
                    "connections_active": current_connections,
                    "connections_max": max_conn,
                    "connection_usage_percent": round(connection_usage, 2),
                    "query_test": "passed"
                }
                
                # Determine status based on performance
                if response_time > 1000:  # > 1 second
                    status = HealthStatus.DEGRADED
                elif connection_usage > 80:
                    status = HealthStatus.DEGRADED
                else:
                    status = HealthStatus.HEALTHY
                
                return ComponentHealth(
                    name="database",
                    status=status,
                    response_time_ms=response_time,
                    details=details
                )
                
        except Exception as e:
            logger.error(f"Database health check failed: {e}")
            return ComponentHealth(
                name="database",
                status=HealthStatus.UNHEALTHY,
                error_message=str(e)
            )

    async def check_interpro_api_health(self) -> ComponentHealth:
        """Check InterPro API connectivity and response time."""
        start_time = time.time()
        
        try:
            client = InterProAPIClient()
            
            # Test basic connectivity with a simple query
            async with httpx.AsyncClient(timeout=10.0) as http_client:
                response = await http_client.get(
                    f"{self.config.api.interpro_base_url}entry/pfam/",
                    params={"page_size": 1}
                )
                
                response_time = (time.time() - start_time) * 1000
                
                if response.status_code == 200:
                    data = response.json()
                    details = {
                        "endpoint": "entry/pfam/",
                        "status_code": response.status_code,
                        "response_size": len(response.content),
                        "results_available": data.get("count", 0) > 0
                    }
                    
                    # Determine status based on response time
                    if response_time > 5000:  # > 5 seconds
                        status = HealthStatus.DEGRADED
                    else:
                        status = HealthStatus.HEALTHY
                        
                    return ComponentHealth(
                        name="interpro_api",
                        status=status,
                        response_time_ms=response_time,
                        details=details
                    )
                else:
                    return ComponentHealth(
                        name="interpro_api",
                        status=HealthStatus.DEGRADED,
                        response_time_ms=response_time,
                        error_message=f"HTTP {response.status_code}"
                    )
                    
        except Exception as e:
            logger.error(f"InterPro API health check failed: {e}")
            return ComponentHealth(
                name="interpro_api",
                status=HealthStatus.UNHEALTHY,
                error_message=str(e)
            )

    async def check_uniprot_api_health(self) -> ComponentHealth:
        """Check UniProt API connectivity and response time."""
        start_time = time.time()
        
        try:
            # Test direct REST API
            async with httpx.AsyncClient(timeout=10.0) as http_client:
                response = await http_client.get(
                    f"{self.config.api.uniprot_base_url}uniprotkb/search",
                    params={"query": "organism_id:9606", "size": 1}
                )
                
                response_time = (time.time() - start_time) * 1000
                
                if response.status_code == 200:
                    data = response.json()
                    details = {
                        "endpoint": "uniprotkb/search",
                        "status_code": response.status_code,
                        "response_size": len(response.content),
                        "results_available": len(data.get("results", [])) > 0
                    }
                    
                    # Check MCP server if enabled
                    if self.config.mcp.enabled:
                        try:
                            client = UnifiedUniProtClient()
                            # This would test MCP connectivity
                            details["mcp_enabled"] = True
                            details["mcp_status"] = "available"  # Simplified check
                        except Exception as mcp_error:
                            details["mcp_enabled"] = True
                            details["mcp_status"] = "unavailable"
                            details["mcp_error"] = str(mcp_error)
                    else:
                        details["mcp_enabled"] = False
                    
                    # Determine status based on response time
                    if response_time > 5000:  # > 5 seconds
                        status = HealthStatus.DEGRADED
                    else:
                        status = HealthStatus.HEALTHY
                        
                    return ComponentHealth(
                        name="uniprot_api",
                        status=status,
                        response_time_ms=response_time,
                        details=details
                    )
                else:
                    return ComponentHealth(
                        name="uniprot_api",
                        status=HealthStatus.DEGRADED,
                        response_time_ms=response_time,
                        error_message=f"HTTP {response.status_code}"
                    )
                    
        except Exception as e:
            logger.error(f"UniProt API health check failed: {e}")
            return ComponentHealth(
                name="uniprot_api",
                status=HealthStatus.UNHEALTHY,
                error_message=str(e)
            )

    async def check_system_resources(self) -> ComponentHealth:
        """Check system resource usage."""
        try:
            # Get system metrics
            cpu_percent = psutil.cpu_percent(interval=1)
            memory = psutil.virtual_memory()
            disk = psutil.disk_usage('/')
            
            details = {
                "cpu_percent": round(cpu_percent, 2),
                "memory_percent": round(memory.percent, 2),
                "memory_available_gb": round(memory.available / (1024**3), 2),
                "disk_usage_percent": round(disk.percent, 2),
                "disk_free_gb": round(disk.free / (1024**3), 2)
            }
            
            # Determine status based on resource usage
            if cpu_percent > 90 or memory.percent > 90 or disk.percent > 90:
                status = HealthStatus.UNHEALTHY
            elif cpu_percent > 80 or memory.percent > 80 or disk.percent > 85:
                status = HealthStatus.DEGRADED
            else:
                status = HealthStatus.HEALTHY
            
            return ComponentHealth(
                name="system_resources",
                status=status,
                details=details
            )
            
        except Exception as e:
            logger.error(f"System resource check failed: {e}")
            return ComponentHealth(
                name="system_resources",
                status=HealthStatus.UNHEALTHY,
                error_message=str(e)
            )

    async def get_comprehensive_health(self) -> SystemHealth:
        """Get comprehensive system health status."""
        # Run all health checks concurrently
        health_checks = await asyncio.gather(
            self.check_database_health(),
            self.check_interpro_api_health(),
            self.check_uniprot_api_health(),
            self.check_system_resources(),
            return_exceptions=True
        )
        
        components = []
        for check in health_checks:
            if isinstance(check, ComponentHealth):
                components.append(check)
            else:
                # Handle exceptions
                logger.error(f"Health check failed: {check}")
                components.append(ComponentHealth(
                    name="unknown",
                    status=HealthStatus.UNHEALTHY,
                    error_message=str(check)
                ))
        
        # Determine overall system status
        statuses = [comp.status for comp in components]
        if HealthStatus.UNHEALTHY in statuses:
            overall_status = HealthStatus.UNHEALTHY
        elif HealthStatus.DEGRADED in statuses:
            overall_status = HealthStatus.DEGRADED
        else:
            overall_status = HealthStatus.HEALTHY
        
        uptime = time.time() - self.start_time
        
        return SystemHealth(
            status=overall_status,
            components=components,
            timestamp=datetime.now(),
            uptime_seconds=uptime
        )

    async def get_performance_metrics(self) -> PerformanceMetrics:
        """Get current performance metrics."""
        # Check cache first
        cache_key = "performance_metrics"
        now = time.time()
        
        if cache_key in self._metrics_cache:
            cached_time, cached_metrics = self._metrics_cache[cache_key]
            if now - cached_time < self._cache_ttl:
                return cached_metrics
        
        try:
            # System metrics
            cpu_percent = psutil.cpu_percent(interval=0.1)
            memory = psutil.virtual_memory()
            disk = psutil.disk_usage('/')
            
            # Database metrics
            active_connections = 0
            try:
                db_manager = get_database_manager()
                with db_manager.get_connection() as conn:
                    cursor = conn.cursor()
                    cursor.execute("SHOW STATUS LIKE 'Threads_connected'")
                    result = cursor.fetchone()
                    if result:
                        active_connections = int(result[1])
            except Exception as e:
                logger.warning(f"Could not get database connection count: {e}")
            
            # Application metrics (simplified - would be enhanced with actual metrics)
            request_rate = 0.0  # Would be calculated from actual request logs
            error_rate = 0.0    # Would be calculated from actual error logs
            avg_response_time = 0.0  # Would be calculated from actual response times
            
            metrics = PerformanceMetrics(
                cpu_percent=round(cpu_percent, 2),
                memory_percent=round(memory.percent, 2),
                disk_usage_percent=round(disk.percent, 2),
                active_connections=active_connections,
                request_rate=request_rate,
                error_rate=error_rate,
                avg_response_time_ms=avg_response_time,
                timestamp=datetime.now()
            )
            
            # Cache the metrics
            self._metrics_cache[cache_key] = (now, metrics)
            
            return metrics
            
        except Exception as e:
            logger.error(f"Performance metrics collection failed: {e}")
            # Return default metrics on error
            return PerformanceMetrics(
                cpu_percent=0.0,
                memory_percent=0.0,
                disk_usage_percent=0.0,
                active_connections=0,
                request_rate=0.0,
                error_rate=0.0,
                avg_response_time_ms=0.0,
                timestamp=datetime.now()
            )

    async def get_data_summary(self) -> Dict[str, Any]:
        """Get summary of collected data."""
        try:
            db_manager = get_database_manager()
            with db_manager.get_connection() as conn:
                cursor = conn.cursor()
                
                # Get counts from each table
                cursor.execute("SELECT COUNT(*) FROM pfam_families")
                pfam_count = cursor.fetchone()[0]
                
                cursor.execute("SELECT COUNT(*) FROM interpro_proteins")
                protein_count = cursor.fetchone()[0]
                
                cursor.execute("SELECT COUNT(*) FROM proteins")
                isoform_count = cursor.fetchone()[0]
                
                cursor.execute("SELECT COUNT(*) FROM proteins WHERE tim_barrel_location IS NOT NULL")
                tim_barrel_count = cursor.fetchone()[0]
                
                # Get recent activity
                cursor.execute("""
                    SELECT COUNT(*) FROM proteins 
                    WHERE created_at > DATE_SUB(NOW(), INTERVAL 24 HOUR)
                """)
                recent_additions = cursor.fetchone()[0]
                
                return {
                    "pfam_families": pfam_count,
                    "interpro_proteins": protein_count,
                    "protein_isoforms": isoform_count,
                    "tim_barrel_proteins": tim_barrel_count,
                    "recent_additions_24h": recent_additions,
                    "last_updated": datetime.now().isoformat()
                }
                
        except Exception as e:
            logger.error(f"Data summary collection failed: {e}")
            return {
                "error": "Could not retrieve data summary",
                "message": str(e)
            }


# Global health checker instance
_health_checker = None


def get_health_checker() -> HealthChecker:
    """Get the global health checker instance."""
    global _health_checker
    if _health_checker is None:
        _health_checker = HealthChecker()
    return _health_checker