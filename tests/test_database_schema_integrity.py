"""
Property-based tests for database schema integrity.

Feature: protein-data-collector, Property 7: Database Storage Completeness
Validates: Requirements 4.1, 4.2, 4.3, 4.5
"""

import pytest
from datetime import datetime
from contextlib import contextmanager
from hypothesis import given, strategies as st, settings, HealthCheck
from sqlalchemy.exc import IntegrityError

from protein_data_collector.database import (
    Base, PfamFamily, InterProProtein, Protein, DatabaseManager
)
from protein_data_collector.config import DatabaseConfig


class _TestDatabaseConfig(DatabaseConfig):
    """Test database configuration that uses SQLite in-memory database."""
    
    @property
    def connection_url(self) -> str:
        """Generate SQLite in-memory connection URL for testing."""
        return "sqlite:///:memory:"


@contextmanager
def get_test_db_manager():
    """Create a test database manager with in-memory SQLite."""
    test_config = _TestDatabaseConfig(
        host="",
        port=0,
        database="",
        username="",
        password="",
        pool_size=1,
        pool_recycle=3600
    )
    manager = DatabaseManager(test_config)
    # Create all tables
    Base.metadata.create_all(manager.engine)
    try:
        yield manager
    finally:
        # Cleanup
        manager.close()


# Hypothesis strategies for generating test data
pfam_accession_strategy = st.text(
    alphabet=st.characters(whitelist_categories=("Lu", "Ll", "Nd")),
    min_size=5,
    max_size=20
).filter(lambda x: x.strip() and not x.isspace())

protein_name_strategy = st.text(
    alphabet=st.characters(whitelist_categories=("Lu", "Ll", "Nd", "Zs")),
    min_size=1,
    max_size=255
).filter(lambda x: x.strip() and not x.isspace())

description_strategy = st.text(
    alphabet=st.characters(whitelist_categories=("Lu", "Ll", "Nd", "Zs", "Po")),
    min_size=1,
    max_size=1000
).filter(lambda x: x.strip() and not x.isspace())

# Valid amino acid characters for protein sequences
amino_acids = "ACDEFGHIKLMNPQRSTVWY"
protein_sequence_strategy = st.text(
    alphabet=amino_acids,
    min_size=10,
    max_size=1000
)

json_data_strategy = st.dictionaries(
    keys=st.text(min_size=1, max_size=20),
    values=st.one_of(
        st.text(min_size=1, max_size=50),
        st.integers(min_value=1, max_value=1000),
        st.lists(st.integers(min_value=1, max_value=100), min_size=1, max_size=5)
    ),
    min_size=1,
    max_size=5
)


class TestDatabaseSchemaIntegrity:
    """Property-based tests for database schema integrity."""
    
    @given(
        accession=pfam_accession_strategy,
        name=protein_name_strategy,
        description=description_strategy,
        tim_barrel_annotation=description_strategy
    )
    @settings(max_examples=100, suppress_health_check=[HealthCheck.function_scoped_fixture])
    def test_pfam_family_storage_completeness(
        self, accession, name, description, tim_barrel_annotation
    ):
        """
        Feature: protein-data-collector, Property 7: Database Storage Completeness
        
        For any successfully validated PFAM family entity, the system should store 
        all required fields in the local database with proper relational integrity maintained.
        """
        with get_test_db_manager() as test_db_manager:
            with test_db_manager.get_transaction() as session:
                # Create and store PFAM family
                pfam_family = PfamFamily(
                    accession=accession,
                    name=name,
                    description=description,
                    tim_barrel_annotation=tim_barrel_annotation
                )
                session.add(pfam_family)
                session.flush()  # Ensure it's written to database
                
                # Verify all required fields are stored
                stored_family = session.query(PfamFamily).filter_by(accession=accession).first()
                assert stored_family is not None
                assert stored_family.accession == accession
                assert stored_family.name == name
                assert stored_family.description == description
                assert stored_family.tim_barrel_annotation == tim_barrel_annotation
                assert stored_family.created_at is not None
                assert isinstance(stored_family.created_at, datetime)