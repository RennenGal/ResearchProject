"""Configuration for the protein data collector."""

import os
from dataclasses import dataclass
from typing import Dict


@dataclass
class DomainConfig:
    """Domain-specific settings."""
    display_name: str                   # e.g. "TIM barrel"
    interpro_annotation: str            # term for InterPro annotation= query (exact match)
    entries_table: str                  # e.g. "entries"
    table_prefix: str                   # e.g. "" — used to derive organism table names
    accession_col: str = ""             # DB column name for domain accession on proteins table
    location_col: str = ""              # DB column name for domain location on isoforms table
    sequence_col: str = ""              # DB column name for domain sequence on isoforms table
    interpro_search: str = ""           # term for InterPro search= query (text match, fallback)
    cathgene3d_search: str = ""         # search= query for CATH Gene3D (catches structurally-classified entries with no IPR parent)
    extra_accessions: tuple = ()        # additional accessions to always include; never removed by cleanup


@dataclass
class OrganismConfig:
    """Organism-specific settings used by the collection pipeline."""
    display_name: str           # e.g. "Homo sapiens"
    taxon_id: int               # NCBI taxonomy ID used in InterPro queries
    organism_suffix: str        # e.g. "" or "_mus_musculus"

    def protein_table(self, domain: DomainConfig) -> str:
        return f"{domain.table_prefix}proteins{self.organism_suffix}"

    def isoform_table(self, domain: DomainConfig) -> str:
        return f"{domain.table_prefix}isoforms{self.organism_suffix}"

    def affected_isoforms_table(self, domain: DomainConfig) -> str:
        return f"{domain.table_prefix}affected_isoforms{self.organism_suffix}"


DOMAINS: Dict[str, DomainConfig] = {
    "tim_barrel": DomainConfig(
        display_name="TIM barrel",
        interpro_annotation="TIM barrel",
        entries_table="entries",
        table_prefix="",
        accession_col="tim_barrel_accession",
        location_col="tim_barrel_location",
        sequence_col="tim_barrel_sequence",
        cathgene3d_search="3.20.20",
        extra_accessions=(
            "IPR011060",
        ),
    ),
}

ORGANISMS: Dict[str, OrganismConfig] = {
    "homo_sapiens": OrganismConfig(
        display_name="Homo sapiens",
        taxon_id=9606,
        organism_suffix="",
    ),
}


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
