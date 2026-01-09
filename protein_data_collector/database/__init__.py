"""
Database layer for the Protein Data Collector system.

This package provides database schema definitions, connection management,
migration utilities, and storage operations for protein data.
"""

from .schema import Base, TIMBarrelEntry, InterProProtein, Protein
from .connection import (
    DatabaseManager, 
    get_database_manager, 
    set_database_manager,
    get_db_session, 
    get_db_transaction
)
from .migrations import (
    Migration, 
    MigrationManager, 
    get_migration_manager,
    migrate_database, 
    get_migration_status
)
from .storage import (
    DatabaseStorage,
    StorageStats,
    StorageResult,
    store_all_entities,
    get_database_statistics
)

__all__ = [
    # Schema
    "Base",
    "TIMBarrelEntry", 
    "InterProProtein", 
    "Protein",
    
    # Connection management
    "DatabaseManager",
    "get_database_manager",
    "set_database_manager", 
    "get_db_session",
    "get_db_transaction",
    
    # Migrations
    "Migration",
    "MigrationManager",
    "get_migration_manager",
    "migrate_database",
    "get_migration_status",
    
    # Storage operations
    "DatabaseStorage",
    "StorageStats",
    "StorageResult",
    "store_all_entities",
    "get_database_statistics",
]