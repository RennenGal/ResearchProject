"""
Configuration management for the Protein Data Collector system.

This module provides configuration classes and utilities for managing
API endpoints, database settings, retry policies, and other system parameters.
"""

import os
from dataclasses import dataclass, field
from typing import Optional
import json


@dataclass
class RetryConfig:
    """Configuration for retry behavior across all external API calls."""
    max_retries: int = 3
    initial_delay: float = 1.0
    backoff_multiplier: float = 2.0
    max_delay: float = 60.0


@dataclass
class DatabaseConfig:
    """Database connection and configuration settings."""
    type: str = "sqlite"  # "mysql" or "sqlite"
    path: str = "db/protein_data.db"  # SQLite database file path
    host: str = "localhost"
    port: int = 3306
    database: str = "protein_data"
    username: str = "protein_user"
    password: str = ""
    pool_size: int = 10
    pool_recycle: int = 3600
    
    @property
    def connection_url(self) -> str:
        """Generate database connection URL."""
        if self.type == "sqlite":
            return f"sqlite:///{self.path}"
        else:
            return f"mysql+pymysql://{self.username}:{self.password}@{self.host}:{self.port}/{self.database}"


@dataclass
class RateLimitingConfig:
    """Rate limiting configuration for external APIs."""
    # InterPro API rate limiting
    interpro_requests_per_second: float = 10.0
    interpro_burst_limit: int = 50
    interpro_burst_window_seconds: int = 60
    
    # UniProt API rate limiting
    uniprot_requests_per_second: float = 5.0
    uniprot_burst_limit: int = 25
    uniprot_burst_window_seconds: int = 60
    

    
    # Exponential backoff for rate limit violations
    violation_initial_delay: float = 1.0
    violation_backoff_multiplier: float = 2.0
    violation_max_delay: float = 300.0  # 5 minutes max
    
    # Monitoring and reporting
    soft_limit_threshold: float = 0.8  # Warn at 80% of rate limit
    enable_monitoring: bool = True
    enable_reporting: bool = True


@dataclass
class APIConfig:
    """External API endpoints and configuration."""
    interpro_base_url: str = "https://www.ebi.ac.uk/interpro/api/"
    uniprot_base_url: str = "https://rest.uniprot.org/"
    
    # Request timeout settings
    request_timeout: int = 30
    connection_timeout: int = 10



@dataclass
class CollectionConfig:
    """Data collection behavior configuration."""
    batch_size: int = 100
    cache_ttl_hours: int = 24
    enable_progress_tracking: bool = True
    resume_interrupted_collections: bool = True
    
    # Performance optimization settings
    enable_response_caching: bool = True
    cache_max_size: int = 10000
    cache_max_memory_mb: int = 500
    enable_connection_pooling: bool = True
    connection_pool_size: int = 20
    enable_performance_monitoring: bool = True


@dataclass
class LoggingConfig:
    """Logging system configuration."""
    level: str = "INFO"
    format: str = "json"
    log_file: Optional[str] = None
    max_file_size_mb: int = 100
    backup_count: int = 5
    structured: bool = True


@dataclass
class SystemConfig:
    """Main system configuration combining all subsystem configs."""
    database: DatabaseConfig = field(default_factory=DatabaseConfig)
    api: APIConfig = field(default_factory=APIConfig)
    retry: RetryConfig = field(default_factory=RetryConfig)
    rate_limiting: RateLimitingConfig = field(default_factory=RateLimitingConfig)
    collection: CollectionConfig = field(default_factory=CollectionConfig)
    logging: LoggingConfig = field(default_factory=LoggingConfig)
    
    @classmethod
    def from_env(cls) -> "SystemConfig":
        """Create configuration from environment variables."""
        config = cls()
        
        # Database configuration from environment
        if os.getenv("DB_HOST"):
            config.database.host = os.getenv("DB_HOST")
        if os.getenv("DB_PORT"):
            config.database.port = int(os.getenv("DB_PORT"))
        if os.getenv("DB_NAME"):
            config.database.database = os.getenv("DB_NAME")
        if os.getenv("DB_USER"):
            config.database.username = os.getenv("DB_USER")
        if os.getenv("DB_PASSWORD"):
            config.database.password = os.getenv("DB_PASSWORD")
            
        # API configuration from environment
        if os.getenv("INTERPRO_BASE_URL"):
            config.api.interpro_base_url = os.getenv("INTERPRO_BASE_URL")
        if os.getenv("UNIPROT_BASE_URL"):
            config.api.uniprot_base_url = os.getenv("UNIPROT_BASE_URL")
            
        # Rate limiting configuration from environment
        if os.getenv("INTERPRO_REQUESTS_PER_SECOND"):
            config.rate_limiting.interpro_requests_per_second = float(os.getenv("INTERPRO_REQUESTS_PER_SECOND"))
        if os.getenv("UNIPROT_REQUESTS_PER_SECOND"):
            config.rate_limiting.uniprot_requests_per_second = float(os.getenv("UNIPROT_REQUESTS_PER_SECOND"))
        if os.getenv("RATE_LIMIT_MONITORING"):
            config.rate_limiting.enable_monitoring = os.getenv("RATE_LIMIT_MONITORING").lower() == "true"
            
        # Retry configuration from environment
        if os.getenv("MAX_RETRIES"):
            config.retry.max_retries = int(os.getenv("MAX_RETRIES"))
        if os.getenv("INITIAL_DELAY"):
            config.retry.initial_delay = float(os.getenv("INITIAL_DELAY"))
        if os.getenv("BACKOFF_MULTIPLIER"):
            config.retry.backoff_multiplier = float(os.getenv("BACKOFF_MULTIPLIER"))

            
        # Logging configuration from environment
        if os.getenv("LOG_LEVEL"):
            config.logging.level = os.getenv("LOG_LEVEL")
        if os.getenv("LOG_FILE"):
            config.logging.log_file = os.getenv("LOG_FILE")
            
        return config
    
    @classmethod
    def from_file(cls, config_path: str) -> "SystemConfig":
        """Load configuration from JSON file."""
        with open(config_path, 'r') as f:
            config_data = json.load(f)
        
        config = cls()
        
        # Update configuration with file data
        if "database" in config_data:
            db_config = config_data["database"]
            for key, value in db_config.items():
                if hasattr(config.database, key):
                    setattr(config.database, key, value)
                    
        if "api" in config_data:
            api_config = config_data["api"]
            for key, value in api_config.items():
                if hasattr(config.api, key):
                    setattr(config.api, key, value)
                    
        if "rate_limiting" in config_data:
            rate_limiting_config = config_data["rate_limiting"]
            for key, value in rate_limiting_config.items():
                if hasattr(config.rate_limiting, key):
                    setattr(config.rate_limiting, key, value)
                    
        if "retry" in config_data:
            retry_config = config_data["retry"]
            for key, value in retry_config.items():
                if hasattr(config.retry, key):
                    setattr(config.retry, key, value)

                    
        if "collection" in config_data:
            collection_config = config_data["collection"]
            for key, value in collection_config.items():
                if hasattr(config.collection, key):
                    setattr(config.collection, key, value)
                    
        if "logging" in config_data:
            logging_config = config_data["logging"]
            for key, value in logging_config.items():
                if hasattr(config.logging, key):
                    setattr(config.logging, key, value)
        
        return config
    
    def to_file(self, config_path: str) -> None:
        """Save configuration to JSON file."""
        config_data = {
            "database": {
                "type": self.database.type,
                "path": self.database.path,
                "host": self.database.host,
                "port": self.database.port,
                "database": self.database.database,
                "username": self.database.username,
                "pool_size": self.database.pool_size,
                "pool_recycle": self.database.pool_recycle
            },
            "api": {
                "interpro_base_url": self.api.interpro_base_url,
                "uniprot_base_url": self.api.uniprot_base_url,
                "request_timeout": self.api.request_timeout,
                "connection_timeout": self.api.connection_timeout
            },
            "rate_limiting": {
                "interpro_requests_per_second": self.rate_limiting.interpro_requests_per_second,
                "interpro_burst_limit": self.rate_limiting.interpro_burst_limit,
                "interpro_burst_window_seconds": self.rate_limiting.interpro_burst_window_seconds,
                "uniprot_requests_per_second": self.rate_limiting.uniprot_requests_per_second,
                "uniprot_burst_limit": self.rate_limiting.uniprot_burst_limit,
                "uniprot_burst_window_seconds": self.rate_limiting.uniprot_burst_window_seconds,
                "violation_initial_delay": self.rate_limiting.violation_initial_delay,
                "violation_backoff_multiplier": self.rate_limiting.violation_backoff_multiplier,
                "violation_max_delay": self.rate_limiting.violation_max_delay,
                "soft_limit_threshold": self.rate_limiting.soft_limit_threshold,
                "enable_monitoring": self.rate_limiting.enable_monitoring,
                "enable_reporting": self.rate_limiting.enable_reporting
            },
            "retry": {
                "max_retries": self.retry.max_retries,
                "initial_delay": self.retry.initial_delay,
                "backoff_multiplier": self.retry.backoff_multiplier,
                "max_delay": self.retry.max_delay

            },
            "collection": {
                "batch_size": self.collection.batch_size,
                "cache_ttl_hours": self.collection.cache_ttl_hours,
                "enable_progress_tracking": self.collection.enable_progress_tracking,
                "resume_interrupted_collections": self.collection.resume_interrupted_collections
            },
            "logging": {
                "level": self.logging.level,
                "format": self.logging.format,
                "log_file": self.logging.log_file,
                "max_file_size_mb": self.logging.max_file_size_mb,
                "backup_count": self.logging.backup_count
            }
        }
        
        with open(config_path, 'w') as f:
            json.dump(config_data, f, indent=2)


# Global configuration instance
_config: Optional[SystemConfig] = None


def get_config() -> SystemConfig:
    """Get the global configuration instance."""
    global _config
    if _config is None:
        _config = SystemConfig.from_env()
    return _config


def set_config(config: SystemConfig) -> None:
    """Set the global configuration instance."""
    global _config
    _config = config


def load_config_from_file(config_path: str) -> SystemConfig:
    """Load and set configuration from file."""
    config = SystemConfig.from_file(config_path)
    set_config(config)
    return config