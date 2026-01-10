"""
Unit tests for the configuration management system.
"""

import pytest
import os
import json
import tempfile
from protein_data_collector.config import (
    SystemConfig, DatabaseConfig, APIConfig, RetryConfig,
    get_config, set_config, load_config_from_file
)


class TestDatabaseConfig:
    """Test database configuration functionality."""
    
    def test_default_values(self):
        """Test default database configuration values."""
        config = DatabaseConfig()
        assert config.host == "localhost"
        assert config.port == 3306
        assert config.database == "protein_data"
        assert config.username == "protein_user"
        assert config.password == ""
        assert config.pool_size == 10
        assert config.pool_recycle == 3600
    
    def test_connection_url_generation(self):
        """Test MySQL connection URL generation."""
        config = DatabaseConfig(
            type="mysql",
            host="testhost",
            port=3307,
            database="testdb",
            username="testuser",
            password="testpass"
        )
        expected_url = "mysql+pymysql://testuser:testpass@testhost:3307/testdb"
        assert config.connection_url == expected_url
    
    def test_connection_url_no_password(self):
        """Test connection URL generation without password."""
        config = DatabaseConfig(
            type="mysql",
            host="testhost",
            username="testuser",
            password=""
        )
        expected_url = "mysql+pymysql://testuser:@testhost:3306/protein_data"
        assert config.connection_url == expected_url
    
    def test_sqlite_connection_url(self):
        """Test SQLite connection URL generation."""
        config = DatabaseConfig(
            type="sqlite",
            path="test_database.db"
        )
        expected_url = "sqlite:///test_database.db"
        assert config.connection_url == expected_url


class TestRetryConfig:
    """Test retry configuration functionality."""
    
    def test_default_values(self):
        """Test default retry configuration values."""
        config = RetryConfig()
        assert config.max_retries == 3
        assert config.initial_delay == 1.0
        assert config.backoff_multiplier == 2.0
        assert config.max_delay == 60.0
    
    def test_custom_values(self):
        """Test custom retry configuration values."""
        config = RetryConfig(
            max_retries=5,
            initial_delay=0.5,
            backoff_multiplier=1.5,
            max_delay=30.0
        )
        assert config.max_retries == 5
        assert config.initial_delay == 0.5
        assert config.backoff_multiplier == 1.5
        assert config.max_delay == 30.0


class TestAPIConfig:
    """Test API configuration functionality."""
    
    def test_default_values(self):
        """Test default API configuration values."""
        config = APIConfig()
        assert config.interpro_base_url == "https://www.ebi.ac.uk/interpro/api/"
        assert config.uniprot_base_url == "https://rest.uniprot.org/"
        assert config.request_timeout == 30
        assert config.connection_timeout == 10


class TestSystemConfig:
    """Test system configuration functionality."""
    
    def test_default_initialization(self):
        """Test default system configuration initialization."""
        config = SystemConfig()
        assert isinstance(config.database, DatabaseConfig)
        assert isinstance(config.api, APIConfig)
        assert isinstance(config.retry, RetryConfig)
        assert config.database.host == "localhost"
        assert config.api.interpro_base_url == "https://www.ebi.ac.uk/interpro/api/"
        assert config.retry.max_retries == 3
    
    def test_from_env(self, mock_env_vars):
        """Test configuration loading from environment variables."""
        config = SystemConfig.from_env()
        
        assert config.database.host == "test_host"
        assert config.database.port == 3307
        assert config.database.database == "test_db"
        assert config.database.username == "test_user"
        assert config.database.password == "test_pass"
        assert config.retry.max_retries == 5
        assert config.logging.level == "DEBUG"
    
    def test_from_file(self, temp_config_file):
        """Test configuration loading from JSON file."""
        config = SystemConfig.from_file(temp_config_file)
        
        assert config.database.host == "localhost"
        assert config.database.database == "test_protein_data"
        assert config.database.username == "test_user"
        assert config.database.password == "test_password"
        assert config.retry.max_retries == 2
        assert config.retry.initial_delay == 0.1
        assert config.logging.level == "DEBUG"
    
    def test_to_file(self, test_config):
        """Test configuration saving to JSON file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            temp_file = f.name
        
        try:
            test_config.to_file(temp_file)
            
            # Verify file was created and contains expected data
            with open(temp_file, 'r') as f:
                saved_data = json.load(f)
            
            assert saved_data["database"]["host"] == "localhost"
            assert saved_data["database"]["database"] == "test_protein_data"
            assert saved_data["retry"]["max_retries"] == 2
            assert saved_data["api"]["interpro_base_url"] == "https://www.ebi.ac.uk/interpro/api/"
        finally:
            os.unlink(temp_file)


class TestGlobalConfig:
    """Test global configuration management."""
    
    def test_get_config_default(self):
        """Test getting default global configuration."""
        # Reset global config
        set_config(None)
        config = get_config()
        assert isinstance(config, SystemConfig)
        assert config.database.host == "localhost"
    
    def test_set_and_get_config(self, test_config):
        """Test setting and getting global configuration."""
        set_config(test_config)
        retrieved_config = get_config()
        assert retrieved_config is test_config
        assert retrieved_config.database.database == "test_protein_data"
    
    def test_load_config_from_file(self, temp_config_file):
        """Test loading global configuration from file."""
        config = load_config_from_file(temp_config_file)
        
        # Verify it was set as global config
        global_config = get_config()
        assert global_config is config
        assert global_config.database.database == "test_protein_data"