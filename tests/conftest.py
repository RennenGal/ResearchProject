"""
Pytest configuration and shared fixtures for the test suite.
"""

import pytest
import tempfile
import os
from pathlib import Path
from protein_data_collector.config import SystemConfig, DatabaseConfig, APIConfig, RetryConfig, RateLimitingConfig, set_config
from protein_data_collector.database.connection import set_database_manager, DatabaseManager


@pytest.fixture
def temp_config_file():
    """Create a temporary configuration file for testing."""
    config_data = {
        "database": {
            "host": "localhost",
            "port": 3306,
            "database": "test_protein_data",
            "username": "test_user",
            "password": "test_password"
        },
        "api": {
            "interpro_base_url": "https://www.ebi.ac.uk/interpro/api/",
            "uniprot_base_url": "https://rest.uniprot.org/"
        },
        "rate_limiting": {
            "interpro_requests_per_second": 10.0,
            "uniprot_requests_per_second": 5.0
        },
        "retry": {
            "max_retries": 2,
            "initial_delay": 0.1,
            "backoff_multiplier": 1.5,
            "max_delay": 5.0
        },
        "logging": {
            "level": "DEBUG",
            "format": "json"
        }
    }
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        import json
        json.dump(config_data, f)
        temp_file = f.name
    
    yield temp_file
    
    # Cleanup
    os.unlink(temp_file)


@pytest.fixture
def test_config():
    """Create a test configuration instance."""
    return SystemConfig(
        database=DatabaseConfig(
            host="localhost",
            port=3306,
            database="test_protein_data",
            username="test_user",
            password="test_password"
        ),
        api=APIConfig(
            interpro_base_url="https://www.ebi.ac.uk/interpro/api/",
            uniprot_base_url="https://rest.uniprot.org/"
        ),
        rate_limiting=RateLimitingConfig(
            interpro_requests_per_second=10.0,
            uniprot_requests_per_second=5.0
        ),
        retry=RetryConfig(
            max_retries=2,
            initial_delay=0.1,
            backoff_multiplier=1.5,
            max_delay=5.0
        )
    )


@pytest.fixture
def mock_env_vars(monkeypatch):
    """Set up mock environment variables for testing."""
    env_vars = {
        "DB_HOST": "test_host",
        "DB_PORT": "3307",
        "DB_NAME": "test_db",
        "DB_USER": "test_user",
        "DB_PASSWORD": "test_pass",
        "MAX_RETRIES": "5",
        "LOG_LEVEL": "DEBUG"
    }
    
    for key, value in env_vars.items():
        monkeypatch.setenv(key, value)
    
    return env_vars


@pytest.fixture(scope="session")
def test_data_dir():
    """Get the test data directory path."""
    return Path(__file__).parent / "data"


@pytest.fixture(autouse=True)
def setup_test_config():
    """Automatically set up test configuration for all tests."""
    # Create test configuration
    test_config = SystemConfig(
        database=DatabaseConfig(
            host="localhost",
            port=3306,
            database="test_protein_data",
            username="test_user",
            password="test_password"
        ),
        api=APIConfig(
            interpro_base_url="https://www.ebi.ac.uk/interpro/api/",
            uniprot_base_url="https://rest.uniprot.org/"
        ),
        rate_limiting=RateLimitingConfig(
            interpro_requests_per_second=10.0,
            uniprot_requests_per_second=5.0
        ),
        retry=RetryConfig(
            max_retries=2,
            initial_delay=0.1,
            backoff_multiplier=1.5,
            max_delay=5.0
        )
    )
    
    # Set as global config
    set_config(test_config)
    
    # Create and set test database manager
    db_manager = DatabaseManager(test_config.database)
    set_database_manager(db_manager)
    
    yield test_config
    
    # Cleanup - reset to None
    set_config(None)
    set_database_manager(None)