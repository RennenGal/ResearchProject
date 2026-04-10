"""Configuration for the protein data collector."""

import os
from dataclasses import dataclass


@dataclass
class Config:
    db_path: str = "db/protein_data.db"
    interpro_base_url: str = "https://www.ebi.ac.uk/interpro/api"
    uniprot_base_url: str = "https://rest.uniprot.org/uniprotkb"
    request_timeout: int = 30
    request_delay: float = 0.1   # seconds between API calls
    max_retries: int = 3
    retry_initial_delay: float = 1.0
    retry_max_delay: float = 60.0
    batch_size: int = 25


_config = Config()


def get_config() -> Config:
    return _config


def load_config_from_env() -> Config:
    """Override defaults from environment variables."""
    config = Config()
    if v := os.getenv("DB_PATH"):
        config.db_path = v
    if v := os.getenv("REQUEST_DELAY"):
        config.request_delay = float(v)
    if v := os.getenv("MAX_RETRIES"):
        config.max_retries = int(v)
    global _config
    _config = config
    return config
