"""Configuration for the protein data collector."""

import os
from dataclasses import dataclass
from typing import Dict


@dataclass
class DomainConfig:
    """Domain-specific settings (e.g. TIM barrel, beta propeller)."""
    display_name: str                   # e.g. "TIM barrel"
    interpro_annotation: str            # term for InterPro annotation= query (exact match)
    entries_table: str                  # e.g. "tb_entries"
    table_prefix: str                   # e.g. "tb_" — used to derive organism tables
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


# Registry of supported domains.  Keys are the CLI --domain argument values.
DOMAINS: Dict[str, DomainConfig] = {
    "tim_barrel": DomainConfig(
        display_name="TIM barrel",
        interpro_annotation="TIM barrel",
        entries_table="tb_entries",
        table_prefix="tb_",
        cathgene3d_search="3.20.20",    # CATH TIM barrel superfamilies (3.20.20.x); catches structurally-confirmed entries with no IPR parent
        extra_accessions=(
            "IPR011060",   # Ribulose-phosphate binding barrel (SSF51366 parent) — not returned by annotation= query
        ),
    ),
    "beta_propeller": DomainConfig(
        display_name="Beta propeller",
        interpro_annotation="",           # no uniform annotation= term in InterPro
        entries_table="bp_entries",
        table_prefix="bp_",
        interpro_search="propeller",      # text search catches named beta-propeller entries
        cathgene3d_search="2.130",        # CATH beta propeller superfamilies (2.130.x.x); catches structurally-confirmed entries with no Pfam entry (e.g. integrins, RCC1)
        extra_accessions=(                # WD40 superfamilies — common beta propeller not named "propeller"
            "PF00400",    # WD domain, G-beta repeat (WD40)
            "IPR001680",  # WD40 repeat
            "IPR036322",  # WD40-repeat-containing domain superfamily
            "IPR015943",  # WD40/YVTN repeat-like-containing domain superfamily
        ),
    ),
}

# Registry of supported organisms.  Keys are the CLI --organism argument values.
ORGANISMS: Dict[str, OrganismConfig] = {
    "homo_sapiens": OrganismConfig(
        display_name="Homo sapiens",
        taxon_id=9606,
        organism_suffix="",
    ),
    "mus_musculus": OrganismConfig(
        display_name="Mus musculus",
        taxon_id=10090,
        organism_suffix="_mus_musculus",
    ),
    "rattus_norvegicus": OrganismConfig(
        display_name="Rattus norvegicus",
        taxon_id=10116,
        organism_suffix="_rattus_norvegicus",
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
