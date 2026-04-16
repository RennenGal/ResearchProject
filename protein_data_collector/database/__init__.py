from .connection import ensure_db, get_connection
from .storage import (
    get_all_domain_entries, get_all_proteins, get_counts,
    get_isoforms_for_protein, get_proteins_without_isoforms,
    upsert_domain_entries, upsert_domain_entry,
    upsert_isoform, upsert_isoforms, upsert_protein, upsert_proteins,
)
