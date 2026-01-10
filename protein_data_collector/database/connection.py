"""
Database connection management and session handling.

This module provides database connection pooling, transaction management,
session lifecycle management, and performance monitoring for the Protein Data Collector system.
"""

import logging
import time
from contextlib import contextmanager
from typing import Generator, Optional
from sqlalchemy import create_engine, Engine, event
from sqlalchemy.orm import sessionmaker, Session
from sqlalchemy.pool import QueuePool
from sqlalchemy.exc import SQLAlchemyError

from ..config import DatabaseConfig, get_config
from ..performance import get_performance_monitor, record_database_performance
from .schema import Base

logger = logging.getLogger(__name__)


class DatabaseManager:
    """Manages database connections, sessions, transactions, and performance monitoring."""
    
    def __init__(self, config: Optional[DatabaseConfig] = None):
        """Initialize database manager with configuration."""
        self.config = config or get_config().database
        self._engine: Optional[Engine] = None
        self._session_factory: Optional[sessionmaker] = None
        self.performance_monitor = get_performance_monitor()
        
        # Performance tracking
        self._query_start_times = {}
    
    def _get_connect_args(self) -> dict:
        """Get database-specific connection arguments."""
        if self.config.type == "sqlite":
            return {
                "check_same_thread": False,
                "timeout": 20
            }
        elif "mysql" in self.config.connection_url:
            return {
                "charset": "utf8mb4",
                "autocommit": False
            }
        else:
            return {}
    
    def _setup_performance_monitoring(self, engine: Engine) -> None:
        """Set up SQLAlchemy event listeners for performance monitoring."""
        
        @event.listens_for(engine, "before_cursor_execute")
        def before_cursor_execute(conn, cursor, statement, parameters, context, executemany):
            """Record query start time."""
            context._query_start_time = time.time()
        
        @event.listens_for(engine, "after_cursor_execute")
        def after_cursor_execute(conn, cursor, statement, parameters, context, executemany):
            """Record query completion and performance metrics."""
            if hasattr(context, '_query_start_time'):
                duration_ms = (time.time() - context._query_start_time) * 1000
                
                # Extract query type and table from statement
                query_type = statement.strip().split()[0].upper()
                
                # Try to extract table name (simplified approach)
                table_name = "unknown"
                statement_lower = statement.lower()
                if "from " in statement_lower:
                    parts = statement_lower.split("from ")[1].split()
                    if parts:
                        table_name = parts[0].strip('`"[]')
                elif "into " in statement_lower:
                    parts = statement_lower.split("into ")[1].split()
                    if parts:
                        table_name = parts[0].strip('`"[]')
                elif "update " in statement_lower:
                    parts = statement_lower.split("update ")[1].split()
                    if parts:
                        table_name = parts[0].strip('`"[]')
                
                # Record performance metrics
                self.performance_monitor.record_database_query(
                    query_type=query_type,
                    table=table_name,
                    duration_ms=duration_ms,
                    rows_affected=cursor.rowcount if hasattr(cursor, 'rowcount') else None
                )
    
    @property
    def engine(self) -> Engine:
        """Get or create the database engine with connection pooling and performance monitoring."""
        if self._engine is None:
            system_config = get_config()
            
            if self.config.type == "sqlite":
                # SQLite-specific configuration
                self._engine = create_engine(
                    self.config.connection_url,
                    echo=False,  # Set to True for SQL debugging
                    future=True,  # Use SQLAlchemy 2.0 style
                    connect_args=self._get_connect_args()
                )
                
                # Enable SQLite optimizations
                @event.listens_for(self._engine, "connect")
                def set_sqlite_pragma(dbapi_connection, connection_record):
                    cursor = dbapi_connection.cursor()
                    cursor.execute("PRAGMA foreign_keys=ON")
                    cursor.execute("PRAGMA journal_mode=WAL")
                    cursor.execute("PRAGMA synchronous=NORMAL")
                    cursor.execute("PRAGMA cache_size=10000")
                    cursor.execute("PRAGMA temp_store=memory")
                    cursor.close()
                
            else:
                # MySQL configuration
                pool_size = system_config.collection.connection_pool_size if system_config.collection.enable_connection_pooling else self.config.pool_size
                
                self._engine = create_engine(
                    self.config.connection_url,
                    poolclass=QueuePool,
                    pool_size=pool_size,
                    pool_recycle=self.config.pool_recycle,
                    pool_pre_ping=True,  # Verify connections before use
                    pool_timeout=30,  # Timeout for getting connection from pool
                    max_overflow=pool_size // 2,  # Allow some overflow connections
                    echo=False,  # Set to True for SQL debugging
                    future=True,  # Use SQLAlchemy 2.0 style
                    connect_args=self._get_connect_args()
                )
            
            # Set up performance monitoring
            if system_config.collection.enable_performance_monitoring:
                self._setup_performance_monitoring(self._engine)
            
            db_info = f"SQLite: {self.config.path}" if self.config.type == "sqlite" else f"{self.config.host}:{self.config.port}/{self.config.database}"
            logger.info(f"Created {self.config.type.upper()} database engine for {db_info}")
        return self._engine
    
    @property
    def session_factory(self) -> sessionmaker:
        """Get or create the session factory."""
        if self._session_factory is None:
            self._session_factory = sessionmaker(
                bind=self.engine,
                expire_on_commit=False,
                future=True
            )
        return self._session_factory
    
    def create_tables(self) -> None:
        """Create all database tables defined in the schema."""
        try:
            Base.metadata.create_all(self.engine)
            logger.info("Successfully created database tables")
        except SQLAlchemyError as e:
            logger.error(f"Failed to create database tables: {e}")
            raise
    
    def drop_tables(self) -> None:
        """Drop all database tables (use with caution)."""
        try:
            Base.metadata.drop_all(self.engine)
            logger.info("Successfully dropped database tables")
        except SQLAlchemyError as e:
            logger.error(f"Failed to drop database tables: {e}")
            raise
    
    @contextmanager
    def get_connection(self):
        """
        Get a raw database connection for direct SQL operations.
        
        Usage:
            with db_manager.get_connection() as conn:
                cursor = conn.cursor()
                cursor.execute("SELECT 1")
        """
        connection = self.engine.raw_connection()
        try:
            yield connection
        except Exception as e:
            connection.rollback()
            logger.error(f"Database connection error, rolling back: {e}")
            raise
        finally:
            connection.close()
    
    @contextmanager
    def get_session(self) -> Generator[Session, None, None]:
        """
        Get a database session with automatic transaction management.
        
        Usage:
            with db_manager.get_session() as session:
                # Use session for database operations
                session.add(entity)
                session.commit()  # Explicit commit
        """
        session = self.session_factory()
        try:
            yield session
        except Exception as e:
            session.rollback()
            logger.error(f"Database session error, rolling back: {e}")
            raise
        finally:
            session.close()
    
    @contextmanager
    def get_transaction(self) -> Generator[Session, None, None]:
        """
        Get a database session with automatic transaction commit/rollback.
        
        Usage:
            with db_manager.get_transaction() as session:
                # Use session for database operations
                session.add(entity)
                # Automatic commit on success, rollback on exception
        """
        session = self.session_factory()
        try:
            yield session
            session.commit()
            logger.debug("Database transaction committed successfully")
        except Exception as e:
            session.rollback()
            logger.error(f"Database transaction failed, rolling back: {e}")
            raise
        finally:
            session.close()
    
    def test_connection(self) -> bool:
        """Test database connectivity."""
        try:
            from sqlalchemy import text
            with self.get_session() as session:
                session.execute(text("SELECT 1"))
            logger.info("Database connection test successful")
            return True
        except Exception as e:
            logger.error(f"Database connection test failed: {e}")
            return False
    
    def get_performance_metrics(self) -> dict:
        """Get database performance metrics."""
        pool = self.engine.pool
        
        return {
            "connection_pool": {
                "size": pool.size(),
                "checked_in": pool.checkedin(),
                "checked_out": pool.checkedout(),
                "overflow": pool.overflow(),
                "invalid": pool.invalid()
            },
            "database_config": {
                "host": self.config.host,
                "port": self.config.port,
                "database": self.config.database,
                "pool_size": self.config.pool_size,
                "pool_recycle": self.config.pool_recycle
            }
        }
    
    def close(self) -> None:
        """Close database connections and clean up resources."""
        if self._engine:
            self._engine.dispose()
            logger.info("Database engine disposed")


# Global database manager instance
_db_manager: Optional[DatabaseManager] = None


def get_database_manager() -> DatabaseManager:
    """Get the global database manager instance."""
    global _db_manager
    if _db_manager is None:
        _db_manager = DatabaseManager()
    return _db_manager


def set_database_manager(manager: DatabaseManager) -> None:
    """Set the global database manager instance."""
    global _db_manager
    _db_manager = manager


@contextmanager
def get_db_session() -> Generator[Session, None, None]:
    """Convenience function to get a database session."""
    with get_database_manager().get_session() as session:
        yield session


@contextmanager
def get_db_transaction() -> Generator[Session, None, None]:
    """Convenience function to get a database transaction."""
    with get_database_manager().get_transaction() as session:
        yield session