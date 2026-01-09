"""
Database migration system for schema updates.

This module provides utilities for managing database schema changes,
including version tracking, migration execution, and rollback capabilities.
"""

import logging
from datetime import datetime
from typing import List, Dict, Any, Optional
from sqlalchemy import Column, String, DateTime, Integer, Text, MetaData, Table
from sqlalchemy.sql import func
from sqlalchemy.exc import SQLAlchemyError

from .connection import DatabaseManager, get_database_manager
from .schema import Base

logger = logging.getLogger(__name__)


class Migration:
    """Represents a single database migration."""
    
    def __init__(self, version: str, description: str):
        self.version = version
        self.description = description
        self.timestamp = datetime.now()
    
    def up(self, db_manager: DatabaseManager) -> None:
        """Apply the migration (implement in subclasses)."""
        raise NotImplementedError("Migration.up() must be implemented")
    
    def down(self, db_manager: DatabaseManager) -> None:
        """Rollback the migration (implement in subclasses)."""
        raise NotImplementedError("Migration.down() must be implemented")
    
    def __repr__(self) -> str:
        return f"<Migration(version='{self.version}', description='{self.description}')>"


class InitialSchemaMigration(Migration):
    """Initial schema creation migration."""
    
    def __init__(self):
        super().__init__("001", "Create initial schema with pfam_families, interpro_proteins, and proteins tables")
    
    def up(self, db_manager: DatabaseManager) -> None:
        """Create all tables defined in the schema."""
        try:
            Base.metadata.create_all(db_manager.engine)
            logger.info("Created initial database schema")
        except SQLAlchemyError as e:
            logger.error(f"Failed to create initial schema: {e}")
            raise
    
    def down(self, db_manager: DatabaseManager) -> None:
        """Drop all tables."""
        try:
            Base.metadata.drop_all(db_manager.engine)
            logger.info("Dropped initial database schema")
        except SQLAlchemyError as e:
            logger.error(f"Failed to drop initial schema: {e}")
            raise


class MigrationManager:
    """Manages database migrations and version tracking."""
    
    def __init__(self, db_manager: Optional[DatabaseManager] = None):
        self.db_manager = db_manager or get_database_manager()
        self.migrations: List[Migration] = [
            InitialSchemaMigration(),
        ]
        self._ensure_migration_table()
    
    def _ensure_migration_table(self) -> None:
        """Ensure the migration tracking table exists."""
        metadata = MetaData()
        migration_table = Table(
            'schema_migrations',
            metadata,
            Column('version', String(50), primary_key=True),
            Column('description', String(255), nullable=False),
            Column('applied_at', DateTime, default=func.now()),
        )
        
        try:
            metadata.create_all(self.db_manager.engine, tables=[migration_table])
            logger.debug("Migration tracking table ensured")
        except SQLAlchemyError as e:
            logger.error(f"Failed to create migration table: {e}")
            raise
    
    def get_applied_migrations(self) -> List[str]:
        """Get list of applied migration versions."""
        try:
            with self.db_manager.get_session() as session:
                result = session.execute(
                    "SELECT version FROM schema_migrations ORDER BY version"
                )
                return [row[0] for row in result.fetchall()]
        except SQLAlchemyError as e:
            logger.error(f"Failed to get applied migrations: {e}")
            return []
    
    def get_pending_migrations(self) -> List[Migration]:
        """Get list of migrations that haven't been applied."""
        applied_versions = set(self.get_applied_migrations())
        return [m for m in self.migrations if m.version not in applied_versions]
    
    def apply_migration(self, migration: Migration) -> None:
        """Apply a single migration."""
        logger.info(f"Applying migration {migration.version}: {migration.description}")
        
        try:
            with self.db_manager.get_transaction() as session:
                # Apply the migration
                migration.up(self.db_manager)
                
                # Record the migration as applied
                session.execute(
                    "INSERT INTO schema_migrations (version, description, applied_at) VALUES (%s, %s, %s)",
                    (migration.version, migration.description, datetime.now())
                )
                
            logger.info(f"Successfully applied migration {migration.version}")
            
        except Exception as e:
            logger.error(f"Failed to apply migration {migration.version}: {e}")
            raise
    
    def rollback_migration(self, migration: Migration) -> None:
        """Rollback a single migration."""
        logger.info(f"Rolling back migration {migration.version}: {migration.description}")
        
        try:
            with self.db_manager.get_transaction() as session:
                # Rollback the migration
                migration.down(self.db_manager)
                
                # Remove the migration record
                session.execute(
                    "DELETE FROM schema_migrations WHERE version = %s",
                    (migration.version,)
                )
                
            logger.info(f"Successfully rolled back migration {migration.version}")
            
        except Exception as e:
            logger.error(f"Failed to rollback migration {migration.version}: {e}")
            raise
    
    def migrate_up(self, target_version: Optional[str] = None) -> None:
        """Apply all pending migrations up to target version."""
        pending = self.get_pending_migrations()
        
        if target_version:
            pending = [m for m in pending if m.version <= target_version]
        
        if not pending:
            logger.info("No pending migrations to apply")
            return
        
        logger.info(f"Applying {len(pending)} pending migrations")
        
        for migration in pending:
            self.apply_migration(migration)
        
        logger.info("All migrations applied successfully")
    
    def migrate_down(self, target_version: str) -> None:
        """Rollback migrations down to target version."""
        applied_versions = self.get_applied_migrations()
        migrations_to_rollback = [
            m for m in reversed(self.migrations) 
            if m.version in applied_versions and m.version > target_version
        ]
        
        if not migrations_to_rollback:
            logger.info("No migrations to rollback")
            return
        
        logger.info(f"Rolling back {len(migrations_to_rollback)} migrations")
        
        for migration in migrations_to_rollback:
            self.rollback_migration(migration)
        
        logger.info("Migrations rolled back successfully")
    
    def get_migration_status(self) -> Dict[str, Any]:
        """Get current migration status."""
        applied = self.get_applied_migrations()
        pending = self.get_pending_migrations()
        
        return {
            "applied_migrations": applied,
            "pending_migrations": [m.version for m in pending],
            "current_version": applied[-1] if applied else None,
            "latest_version": self.migrations[-1].version if self.migrations else None,
            "is_up_to_date": len(pending) == 0
        }


# Global migration manager instance
_migration_manager: Optional[MigrationManager] = None


def get_migration_manager() -> MigrationManager:
    """Get the global migration manager instance."""
    global _migration_manager
    if _migration_manager is None:
        _migration_manager = MigrationManager()
    return _migration_manager


def migrate_database() -> None:
    """Apply all pending database migrations."""
    manager = get_migration_manager()
    manager.migrate_up()


def get_migration_status() -> Dict[str, Any]:
    """Get current database migration status."""
    manager = get_migration_manager()
    return manager.get_migration_status()