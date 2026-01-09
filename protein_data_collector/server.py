"""
Server runner for the Protein Data Collector REST API.

This module provides a simple way to start the FastAPI server with proper configuration.
"""

import uvicorn
import logging
from .config import get_config
from .logging_config import setup_logging


def run_server(
    host: str = "0.0.0.0",
    port: int = 8000,
    reload: bool = False,
    workers: int = 1
):
    """
    Run the FastAPI server.
    
    Args:
        host: Host to bind to
        port: Port to bind to
        reload: Enable auto-reload for development
        workers: Number of worker processes
    """
    # Load configuration and setup logging
    config = get_config()
    setup_logging(config.logging)
    
    logger = logging.getLogger(__name__)
    logger.info(f"Starting Protein Data Collector API server on {host}:{port}")
    
    # Run the server
    uvicorn.run(
        "protein_data_collector.api.rest_api:app",
        host=host,
        port=port,
        reload=reload,
        workers=workers,
        log_level=config.logging.level.lower(),
        access_log=True
    )


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Protein Data Collector API Server")
    parser.add_argument("--host", default="0.0.0.0", help="Host to bind to")
    parser.add_argument("--port", type=int, default=8000, help="Port to bind to")
    parser.add_argument("--reload", action="store_true", help="Enable auto-reload")
    parser.add_argument("--workers", type=int, default=1, help="Number of workers")
    
    args = parser.parse_args()
    
    run_server(
        host=args.host,
        port=args.port,
        reload=args.reload,
        workers=args.workers
    )