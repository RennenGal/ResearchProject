"""Tests for database storage and retrieval against the real schema."""

import pytest

from protein_data_collector.database.connection import get_connection
from protein_data_collector.database.schema import init_db
from protein_data_collector.database.storage import (
    get_all_proteins,
    get_all_tim_barrel_entries,
    get_counts,
    get_isoforms_for_protein,
    get_proteins_without_isoforms,
    upsert_isoform,
    upsert_protein,
    upsert_tim_barrel_entry,
)
from protein_data_collector.models.entities import Isoform, Protein, TIMBarrelEntry


@pytest.fixture
def db(tmp_path):
    """In-memory-style SQLite database initialised with the new schema."""
    import sqlite3
    path = str(tmp_path / "test.db")
    conn = sqlite3.connect(path)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    init_db(conn)
    conn.close()
    return path


@pytest.fixture
def entry():
    return TIMBarrelEntry(
        accession="PF00394",
        entry_type="pfam",
        name="TIM barrel",
        tim_barrel_annotation="TIM barrel domain",
    )


@pytest.fixture
def protein():
    return Protein(
        uniprot_id="P04637",
        tim_barrel_accession="PF00394",
        protein_name="Cellular tumor antigen p53",
        gene_name="TP53",
        organism="Homo sapiens",
        reviewed=True,
        annotation_score=5,
    )


_SEQ300 = "ACDEFGHIKLMNPQRSTVWY" * 15  # 300 aa — above fragment threshold


@pytest.fixture
def isoform():
    return Isoform(
        isoform_id="P04637-1",
        uniprot_id="P04637",
        is_canonical=True,
        sequence=_SEQ300,
        sequence_length=300,
        tim_barrel_location={"domain_id": "PF00394", "start": 50, "end": 280,
                             "length": 231, "source": "interpro_api"},
        ensembl_gene_id="ENSG00000141510",
    )


@pytest.fixture
def fragment_isoform():
    short_seq = "ACDEFGHIKLMNPQRSTVWY" * 5  # 100 aa — below fragment threshold
    return Isoform(
        isoform_id="P04637-2",
        uniprot_id="P04637",
        is_canonical=False,
        sequence=short_seq,
        sequence_length=100,
    )


class TestTIMBarrelEntries:
    def test_insert_and_retrieve(self, db, entry):
        with get_connection(db) as conn:
            upsert_tim_barrel_entry(conn, entry)
            conn.commit()
            rows = get_all_tim_barrel_entries(conn)
        assert len(rows) == 1
        assert rows[0]["accession"] == "PF00394"
        assert rows[0]["entry_type"] == "pfam"

    def test_upsert_is_idempotent(self, db, entry):
        with get_connection(db) as conn:
            upsert_tim_barrel_entry(conn, entry)
            upsert_tim_barrel_entry(conn, entry)
            conn.commit()
            assert len(get_all_tim_barrel_entries(conn)) == 1


class TestProteins:
    def test_insert_and_retrieve(self, db, entry, protein):
        with get_connection(db) as conn:
            upsert_tim_barrel_entry(conn, entry)
            conn.commit()
            upsert_protein(conn, protein)
            conn.commit()
            rows = get_all_proteins(conn)
        assert len(rows) == 1
        assert rows[0]["uniprot_id"] == "P04637"
        assert rows[0]["reviewed"] == 1

    def test_foreign_key_enforced(self, db, protein):
        with get_connection(db) as conn:
            import sqlite3
            with pytest.raises(sqlite3.IntegrityError):
                upsert_protein(conn, protein)
                conn.commit()


class TestIsoforms:
    def _seed(self, db, entry, protein):
        with get_connection(db) as conn:
            upsert_tim_barrel_entry(conn, entry)
            conn.commit()
            upsert_protein(conn, protein)
            conn.commit()

    def test_insert_and_retrieve(self, db, entry, protein, isoform):
        self._seed(db, entry, protein)
        with get_connection(db) as conn:
            upsert_isoform(conn, isoform)
            conn.commit()
            rows = get_isoforms_for_protein(conn, "P04637")
        assert len(rows) == 1
        assert rows[0]["isoform_id"] == "P04637-1"
        assert rows[0]["is_canonical"] == 1

    def test_tim_barrel_location_round_trip(self, db, entry, protein, isoform):
        self._seed(db, entry, protein)
        with get_connection(db) as conn:
            upsert_isoform(conn, isoform)
            conn.commit()
            rows = get_isoforms_for_protein(conn, "P04637")
        import json
        loc = json.loads(rows[0]["tim_barrel_location"])
        assert loc["start"] == 50
        assert loc["source"] == "interpro_api"

    def test_tim_barrel_sequence_stored(self, db, entry, protein, isoform):
        self._seed(db, entry, protein)
        with get_connection(db) as conn:
            upsert_isoform(conn, isoform)
            conn.commit()
            rows = get_isoforms_for_protein(conn, "P04637")
        assert rows[0]["tim_barrel_sequence"] == _SEQ300[49:280]
        assert rows[0]["is_fragment"] == 0

    def test_fragment_flagged(self, db, entry, protein, fragment_isoform):
        self._seed(db, entry, protein)
        with get_connection(db) as conn:
            upsert_isoform(conn, fragment_isoform)
            conn.commit()
            rows = get_isoforms_for_protein(conn, "P04637")
        assert rows[0]["is_fragment"] == 1
        assert rows[0]["tim_barrel_sequence"] is None

    def test_proteins_without_isoforms(self, db, entry, protein):
        self._seed(db, entry, protein)
        with get_connection(db) as conn:
            remaining = get_proteins_without_isoforms(conn)
        assert "P04637" in remaining

    def test_proteins_without_isoforms_clears_after_insert(self, db, entry, protein, isoform):
        self._seed(db, entry, protein)
        with get_connection(db) as conn:
            upsert_isoform(conn, isoform)
            conn.commit()
            remaining = get_proteins_without_isoforms(conn)
        assert "P04637" not in remaining


class TestCounts:
    def test_counts(self, db, entry, protein, isoform):
        with get_connection(db) as conn:
            upsert_tim_barrel_entry(conn, entry)
            conn.commit()
            upsert_protein(conn, protein)
            conn.commit()
            upsert_isoform(conn, isoform)
            conn.commit()
            counts = get_counts(conn)
        assert counts["tim_barrel_entries"] == 1
        assert counts["proteins"] == 1
        assert counts["isoforms"] == 1
        assert counts["alternative_isoforms"] == 0
