"""
REST API endpoints for the Protein Data Collector system.

This module provides FastAPI endpoints for collection triggering, status monitoring,
querying, and data export with filtering and pagination support.
"""

import asyncio
import logging
from datetime import datetime
from typing import List, Dict, Any, Optional, Union
from pathlib import Path

from fastapi import FastAPI, HTTPException, Query, BackgroundTasks, Depends, status
from fastapi.responses import StreamingResponse, PlainTextResponse
from pydantic import BaseModel, Field
from enum import Enum
import json
import io

from contextlib import asynccontextmanager

from ..collector.data_collector import DataCollector, CollectionReport, run_complete_collection
from ..query.engine import QueryEngine, QueryFilters, QueryResult
from ..query.export import DataExporter, ExportFormat, ExportOptions
from ..config import get_config, SystemConfig
from ..database.connection import get_database_manager
from ..monitoring import get_health_checker, SystemHealth, PerformanceMetrics
from ..alerting import get_alert_manager, Alert

logger = logging.getLogger(__name__)

@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application lifespan context manager."""
    # Startup
    logger.info("Protein Data Collector API starting up")
    
    # Test database connection
    try:
        db_manager = get_database_manager()
        if db_manager.test_connection():
            logger.info("Database connection successful")
        else:
            logger.warning("Database connection failed")
    except Exception as e:
        logger.error(f"Database connection error: {e}")
    
    yield
    
    # Shutdown
    logger.info("Protein Data Collector API shutting down")
    
    # Clean up active collections
    active_collections.clear()

# Pydantic models for API requests/responses

class CollectionRequest(BaseModel):
    """Request model for starting data collection."""
    tim_barrel_query: str = Field(default="TIM barrel", description="Query string for TIM barrel annotations")
    page_size: int = Field(default=200, ge=1, le=1000, description="Number of results per page for API queries")
    batch_size: Optional[int] = Field(default=None, ge=1, le=1000, description="Batch size for processing operations")
    store_data: bool = Field(default=True, description="Whether to store collected data in database")
    progress_file: Optional[str] = Field(default=None, description="Progress file path for resume functionality")


class CollectionStatus(BaseModel):
    """Response model for collection status."""
    current_phase: str
    progress: Dict[str, Any]
    batch_size: int
    progress_file: Optional[str]
    database_connected: bool


class CollectionSummary(BaseModel):
    """Response model for collection summary."""
    collection_phase: str
    duration_seconds: float
    pfam_families_collected: int
    interpro_proteins_collected: int
    uniprot_isoforms_collected: int
    entities_stored: int
    success_rate: float
    validation_errors: int
    api_errors: int
    storage_errors: int
    total_errors: int


class QueryRequest(BaseModel):
    """Request model for protein queries."""
    pfam_family: Optional[str] = Field(default=None, description="PFAM family identifier")
    protein_id: Optional[str] = Field(default=None, description="UniProt protein identifier")
    organism: Optional[str] = Field(default=None, description="Organism filter")
    min_sequence_length: Optional[int] = Field(default=None, ge=1, description="Minimum sequence length")
    max_sequence_length: Optional[int] = Field(default=None, ge=1, description="Maximum sequence length")
    min_exon_count: Optional[int] = Field(default=None, ge=0, description="Minimum exon count")
    max_exon_count: Optional[int] = Field(default=None, ge=0, description="Maximum exon count")
    has_tim_barrel: Optional[bool] = Field(default=None, description="Filter for TIM barrel presence")
    tim_barrel_confidence: Optional[float] = Field(default=None, ge=0.0, le=1.0, description="Minimum TIM barrel confidence")
    include_isoforms: bool = Field(default=True, description="Include isoform details")
    limit: Optional[int] = Field(default=100, ge=1, le=10000, description="Maximum number of results")


class ExportRequest(BaseModel):
    """Request model for data export."""
    format: str = Field(..., pattern="^(fasta|json|csv)$", description="Export format")
    query: QueryRequest = Field(default_factory=QueryRequest, description="Query filters")
    include_sequences: bool = Field(default=True, description="Include protein sequences")
    include_metadata: bool = Field(default=True, description="Include metadata")
    include_annotations: bool = Field(default=True, description="Include annotations")
    max_sequence_length: Optional[int] = Field(default=None, ge=1, description="Truncate sequences to this length")
    pretty_json: bool = Field(default=True, description="Pretty format JSON output")


class SystemStatus(BaseModel):
    """Response model for system status."""
    database_connected: bool
    interpro_api_available: bool
    uniprot_api_available: bool
    configuration: Dict[str, Any]
    data_summary: Dict[str, Any]


class ErrorResponse(BaseModel):
    """Response model for errors."""
    error: str
    detail: Optional[str] = None
    timestamp: datetime = Field(default_factory=datetime.now)


# Global variables for background tasks
active_collections: Dict[str, Dict[str, Any]] = {}


# FastAPI app instance
app = FastAPI(
    title="Protein Data Collector API",
    description="REST API for TIM barrel protein data collection and analysis",
    version="1.0.0",
    docs_url="/docs",
    redoc_url="/redoc",
    lifespan=lifespan
)


# Dependency to get system configuration
def get_system_config() -> SystemConfig:
    """Get system configuration."""
    return get_config()


# Dependency to get query engine
def get_query_engine() -> QueryEngine:
    """Get query engine instance."""
    return QueryEngine()


# Dependency to get data exporter
def get_data_exporter() -> DataExporter:
    """Get data exporter instance."""
    return DataExporter()


@app.get("/", response_model=Dict[str, str])
async def root():
    """Root endpoint with API information."""
    return {
        "name": "Protein Data Collector API",
        "version": "1.0.0",
        "description": "REST API for TIM barrel protein data collection and analysis",
        "docs": "/docs",
        "health": "/health"
    }


@app.get("/health", response_model=Dict[str, str])
async def health_check():
    """Basic health check endpoint."""
    try:
        # Test database connection
        db_manager = get_database_manager()
        db_connected = db_manager.test_connection()
        
        if db_connected:
            return {"status": "healthy", "database": "connected"}
        else:
            return {"status": "degraded", "database": "disconnected"}
    except Exception as e:
        logger.error(f"Health check failed: {e}")
        return {"status": "unhealthy", "error": str(e)}


@app.get("/health/detailed", response_model=Dict[str, Any])
async def detailed_health_check():
    """Comprehensive health check endpoint."""
    try:
        health_checker = get_health_checker()
        system_health = await health_checker.get_comprehensive_health()
        return system_health.to_dict()
    except Exception as e:
        logger.error(f"Detailed health check failed: {e}")
        return {
            "status": "unhealthy",
            "error": str(e),
            "timestamp": datetime.now().isoformat()
        }


@app.get("/health/metrics", response_model=Dict[str, Any])
async def get_health_metrics():
    """Get performance and system metrics."""
    try:
        health_checker = get_health_checker()
        metrics = await health_checker.get_performance_metrics()
        return metrics.to_dict()
    except Exception as e:
        logger.error(f"Metrics collection failed: {e}")
        return {
            "error": "Metrics collection failed",
            "message": str(e),
            "timestamp": datetime.now().isoformat()
        }


@app.get("/health/readiness")
async def readiness_check():
    """Kubernetes readiness probe endpoint."""
    try:
        health_checker = get_health_checker()
        system_health = await health_checker.get_comprehensive_health()
        
        # Check if critical components are healthy
        critical_components = ["database"]
        for component in system_health.components:
            if component.name in critical_components:
                if component.status.value != "healthy":
                    return PlainTextResponse(
                        content="Not Ready",
                        status_code=503
                    )
        
        return PlainTextResponse(content="Ready", status_code=200)
        
    except Exception as e:
        logger.error(f"Readiness check failed: {e}")
        return PlainTextResponse(
            content="Not Ready",
            status_code=503
        )


@app.get("/health/liveness")
async def liveness_check():
    """Kubernetes liveness probe endpoint."""
    try:
        # Simple check - if the application is running, it's alive
        return PlainTextResponse(content="Alive", status_code=200)
    except Exception as e:
        logger.error(f"Liveness check failed: {e}")
        return PlainTextResponse(
            content="Not Alive",
            status_code=503
        )


@app.get("/status", response_model=SystemStatus)
async def get_system_status(config: SystemConfig = Depends(get_system_config)):
    """Get comprehensive system status."""
    try:
        # Database status
        db_manager = get_database_manager()
        db_connected = db_manager.test_connection()
        
        # API availability (simplified check)
        interpro_available = True  # TODO: Add actual API check
        uniprot_available = True   # TODO: Add actual API check
        
        # Data summary
        if db_connected:
            health_checker = get_health_checker()
            data_summary = await health_checker.get_data_summary()
        else:
            data_summary = {}
        
        # Configuration summary
        config_summary = {
            "database_host": config.database.host,
            "database_port": config.database.port,
            "database_name": config.database.database,
            "interpro_base_url": config.api.interpro_base_url,
            "uniprot_base_url": config.api.uniprot_base_url,
            "max_retries": config.retry.max_retries,
            "log_level": config.logging.level
        }
        
        return SystemStatus(
            database_connected=db_connected,
            interpro_api_available=interpro_available,
            uniprot_api_available=uniprot_available,
            configuration=config_summary,
            data_summary=data_summary
        )
        
    except Exception as e:
        logger.error(f"Status check failed: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Status check failed: {str(e)}"
        )


# Collection endpoints

async def run_collection_background(collection_id: str, request: CollectionRequest):
    """Background task for running data collection."""
    try:
        logger.info(f"Starting collection {collection_id}")
        active_collections[collection_id]["status"] = "running"
        active_collections[collection_id]["start_time"] = datetime.now()
        
        # Run collection
        report = await run_complete_collection(
            progress_file=request.progress_file,
            page_size=request.page_size,
            batch_size=request.batch_size,
            store_data=request.store_data
        )
        
        # Update status
        active_collections[collection_id]["status"] = "completed"
        active_collections[collection_id]["end_time"] = datetime.now()
        active_collections[collection_id]["report"] = report.to_summary_dict()
        
        logger.info(f"Collection {collection_id} completed successfully")
        
    except Exception as e:
        logger.error(f"Collection {collection_id} failed: {e}")
        active_collections[collection_id]["status"] = "failed"
        active_collections[collection_id]["error"] = str(e)
        active_collections[collection_id]["end_time"] = datetime.now()


@app.post("/collection/start", response_model=Dict[str, str])
async def start_collection(
    request: CollectionRequest,
    background_tasks: BackgroundTasks
):
    """Start a new data collection process."""
    try:
        # Generate collection ID
        collection_id = f"collection_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        
        # Initialize collection tracking
        active_collections[collection_id] = {
            "id": collection_id,
            "status": "starting",
            "request": request.model_dump(),
            "created_at": datetime.now()
        }
        
        # Start background task
        background_tasks.add_task(run_collection_background, collection_id, request)
        
        return {
            "collection_id": collection_id,
            "status": "started",
            "message": "Collection started in background"
        }
        
    except Exception as e:
        logger.error(f"Failed to start collection: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to start collection: {str(e)}"
        )


@app.get("/collection/{collection_id}/status", response_model=Dict[str, Any])
async def get_collection_status(collection_id: str):
    """Get status of a specific collection."""
    if collection_id not in active_collections:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Collection {collection_id} not found"
        )
    
    collection_info = active_collections[collection_id].copy()
    
    # Add duration if applicable
    if "start_time" in collection_info:
        start_time = collection_info["start_time"]
        end_time = collection_info.get("end_time", datetime.now())
        collection_info["duration_seconds"] = (end_time - start_time).total_seconds()
    
    return collection_info


@app.get("/collection/list", response_model=List[Dict[str, Any]])
async def list_collections():
    """List all collections."""
    collections = []
    for collection_id, info in active_collections.items():
        collection_summary = {
            "id": collection_id,
            "status": info["status"],
            "created_at": info["created_at"]
        }
        
        if "start_time" in info:
            collection_summary["start_time"] = info["start_time"]
        if "end_time" in info:
            collection_summary["end_time"] = info["end_time"]
        if "error" in info:
            collection_summary["error"] = info["error"]
        
        collections.append(collection_summary)
    
    return sorted(collections, key=lambda x: x["created_at"], reverse=True)


# Query endpoints

@app.post("/query/proteins", response_model=QueryResult)
async def query_proteins(
    request: QueryRequest,
    engine: QueryEngine = Depends(get_query_engine)
):
    """Query proteins with filtering options."""
    try:
        # Convert request to QueryFilters
        filters = QueryFilters(
            pfam_family=request.pfam_family,
            protein_id=request.protein_id,
            organism=request.organism,
            min_sequence_length=request.min_sequence_length,
            max_sequence_length=request.max_sequence_length,
            min_exon_count=request.min_exon_count,
            max_exon_count=request.max_exon_count,
            has_tim_barrel=request.has_tim_barrel,
            tim_barrel_confidence=request.tim_barrel_confidence
        )
        
        # Execute query based on specific criteria
        if request.pfam_family:
            result = engine.search_by_pfam_family(request.pfam_family, request.include_isoforms)
        elif request.protein_id:
            result = engine.get_protein_isoforms(request.protein_id)
        elif request.has_tim_barrel:
            criteria = {'has_location': True}
            if request.tim_barrel_confidence:
                criteria['min_confidence'] = request.tim_barrel_confidence
            result = engine.search_by_tim_barrel_features(criteria)
        else:
            result = engine.filter_proteins(filters, request.limit)
        
        return result
        
    except Exception as e:
        logger.error(f"Query failed: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Query failed: {str(e)}"
        )


@app.get("/query/pfam/{pfam_id}", response_model=QueryResult)
async def get_pfam_family(
    pfam_id: str,
    include_isoforms: bool = Query(default=True, description="Include isoform details"),
    engine: QueryEngine = Depends(get_query_engine)
):
    """Get proteins for a specific PFAM family."""
    try:
        result = engine.search_by_pfam_family(pfam_id, include_isoforms)
        
        if result.total_count == 0:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"PFAM family {pfam_id} not found"
            )
        
        return result
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"PFAM query failed: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"PFAM query failed: {str(e)}"
        )


@app.get("/query/protein/{protein_id}", response_model=QueryResult)
async def get_protein_isoforms(
    protein_id: str,
    engine: QueryEngine = Depends(get_query_engine)
):
    """Get all isoforms for a specific protein."""
    try:
        result = engine.get_protein_isoforms(protein_id)
        
        if result.total_count == 0:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Protein {protein_id} not found"
            )
        
        return result
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Protein query failed: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Protein query failed: {str(e)}"
        )


@app.get("/query/statistics", response_model=Dict[str, Any])
async def get_statistics(engine: QueryEngine = Depends(get_query_engine)):
    """Get database statistics."""
    try:
        return engine.get_summary_statistics()
    except Exception as e:
        logger.error(f"Statistics query failed: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Statistics query failed: {str(e)}"
        )


# Export endpoints

@app.post("/export/proteins")
async def export_proteins(
    request: ExportRequest,
    engine: QueryEngine = Depends(get_query_engine),
    exporter: DataExporter = Depends(get_data_exporter)
):
    """Export protein data in specified format."""
    try:
        # Convert query request to filters
        filters = QueryFilters(
            pfam_family=request.query.pfam_family,
            protein_id=request.query.protein_id,
            organism=request.query.organism,
            min_sequence_length=request.query.min_sequence_length,
            max_sequence_length=request.query.max_sequence_length,
            min_exon_count=request.query.min_exon_count,
            max_exon_count=request.query.max_exon_count,
            has_tim_barrel=request.query.has_tim_barrel,
            tim_barrel_confidence=request.query.tim_barrel_confidence
        )
        
        # Query data
        if request.query.pfam_family:
            result = engine.search_by_pfam_family(request.query.pfam_family, True)
            proteins = result.isoforms
        elif request.query.protein_id:
            result = engine.get_protein_isoforms(request.query.protein_id)
            proteins = result.isoforms
        elif request.query.has_tim_barrel:
            criteria = {'has_location': True}
            if request.query.tim_barrel_confidence:
                criteria['min_confidence'] = request.query.tim_barrel_confidence
            result = engine.search_by_tim_barrel_features(criteria)
            proteins = result.isoforms
        else:
            result = engine.filter_proteins(filters, request.query.limit)
            proteins = result.isoforms
        
        if not proteins:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail="No proteins found matching criteria"
            )
        
        # Set up export options
        export_format = ExportFormat(request.format)
        options = ExportOptions(
            format=export_format,
            include_sequences=request.include_sequences,
            include_metadata=request.include_metadata,
            include_annotations=request.include_annotations,
            max_sequence_length=request.max_sequence_length,
            pretty_json=request.pretty_json
        )
        
        # Export data
        exported_data = exporter.export_proteins(proteins, export_format, options)
        
        # Determine content type and filename
        if request.format == "fasta":
            media_type = "text/plain"
            filename = "proteins.fasta"
        elif request.format == "json":
            media_type = "application/json"
            filename = "proteins.json"
        elif request.format == "csv":
            media_type = "text/csv"
            filename = "proteins.csv"
        else:
            media_type = "text/plain"
            filename = "proteins.txt"
        
        # Return as streaming response
        return StreamingResponse(
            io.StringIO(exported_data),
            media_type=media_type,
            headers={"Content-Disposition": f"attachment; filename={filename}"}
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Export failed: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Export failed: {str(e)}"
        )


@app.get("/alerts", response_model=List[Dict[str, Any]])
async def get_active_alerts():
    """Get all active alerts."""
    try:
        alert_manager = get_alert_manager()
        alerts = alert_manager.get_active_alerts()
        return [alert.to_dict() for alert in alerts]
    except Exception as e:
        logger.error(f"Failed to get alerts: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to get alerts: {str(e)}"
        )


@app.get("/alerts/history", response_model=List[Dict[str, Any]])
async def get_alert_history(hours: int = Query(default=24, ge=1, le=168)):
    """Get alert history for the specified number of hours."""
    try:
        alert_manager = get_alert_manager()
        alerts = alert_manager.get_alert_history(hours)
        return [alert.to_dict() for alert in alerts]
    except Exception as e:
        logger.error(f"Failed to get alert history: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to get alert history: {str(e)}"
        )


@app.post("/alerts/{alert_id}/resolve")
async def resolve_alert(alert_id: str):
    """Resolve an active alert."""
    try:
        alert_manager = get_alert_manager()
        await alert_manager.resolve_alert(alert_id)
        return {"message": f"Alert {alert_id} resolved successfully"}
    except Exception as e:
        logger.error(f"Failed to resolve alert {alert_id}: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to resolve alert: {str(e)}"
        )


@app.get("/monitoring/logs", response_model=Dict[str, Any])
async def get_logging_metrics():
    """Get logging metrics and statistics."""
    try:
        from ..logging_config import get_logging_metrics
        metrics = get_logging_metrics()
        return metrics
    except Exception as e:
        logger.error(f"Failed to get logging metrics: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to get logging metrics: {str(e)}"
        )



# Error handlers

@app.exception_handler(HTTPException)
async def http_exception_handler(request, exc):
    """Handle HTTP exceptions."""
    from fastapi.responses import JSONResponse
    return JSONResponse(
        status_code=exc.status_code,
        content={
            "error": exc.detail,
            "detail": f"HTTP {exc.status_code}",
            "timestamp": datetime.now().isoformat()
        }
    )


@app.exception_handler(Exception)
async def general_exception_handler(request, exc):
    """Handle general exceptions."""
    from fastapi.responses import JSONResponse
    logger.error(f"Unhandled exception: {exc}")
    return JSONResponse(
        status_code=500,
        content={
            "error": "Internal server error",
            "detail": str(exc),
            "timestamp": datetime.now().isoformat()
        }
    )


if __name__ == "__main__":
    import uvicorn
    
    # Load configuration
    config = get_config()
    
    # Run the API server
    uvicorn.run(
        "protein_data_collector.api.rest_api:app",
        host="0.0.0.0",
        port=8000,
        reload=False,
        log_level=config.logging.level.lower()
    )