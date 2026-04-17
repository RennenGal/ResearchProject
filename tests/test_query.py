"""Tests for the QueryEngine against a seeded in-memory database."""

import json
import pytest

from protein_data_collector.database.connection import get_connection
from protein_data_collector.database.schema import init_db
from protein_data_collector.database.storage import (
    upsert_domain_entry, upsert_isoform, upsert_protein,
)
from protein_data_collector.models.entities import Isoform, Protein, TIMBarrelEntry
from protein_data_collector.query.engine import QueryEngine
from protein_data_collector.query.export import to_csv, to_fasta, to_json


@pytest.fixture
def db(tmp_path):
    path = str(tmp_path / "test.db")
    import sqlite3
    conn = sqlite3.connect(path)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    init_db(conn)
    conn.close()
    return path


@pytest.fixture
def seeded_db(db):
    """Database with one family, one protein, two isoforms (canonical + alternative)."""
    entry = TIMBarrelEntry(accession="PF00394", entry_type="pfam",
                           name="TIM barrel", domain_annotation="TIM barrel")
    protein = Protein(uniprot_id="P04637", tim_barrel_accession="PF00394",
                      protein_name="p53", gene_name="TP53", reviewed=True)
    canonical = Isoform(
        isoform_id="P04637-1", uniprot_id="P04637", is_canonical=True,
        sequence="ACDEFGHIKLMNPQRSTVWY", sequence_length=20,
        tim_barrel_location={"start": 1, "end": 18, "source": "interpro_api"},
        ensembl_gene_id="ENSG00000141510",
    )
    alt = Isoform(
        isoform_id="P04637-2", uniprot_id="P04637", is_canonical=False,
        sequence="ACDEFGHIKLMNPQRST", sequence_length=17,
        splice_variants=[{"featureId": "VSP_006535",
                          "location": {"start": {"value": 19}, "end": {"value": 20}},
                          "description": "WY -> deleted in isoform 2"}],
    )
    with get_connection(db) as conn:
        upsert_domain_entry(conn, entry)
        conn.commit()
        upsert_protein(conn, protein)
        conn.commit()
        upsert_isoform(conn, canonical)
        upsert_isoform(conn, alt)
        conn.commit()
    return db


class TestQueryEngine:
    def test_get_all_families(self, seeded_db):
        q = QueryEngine(seeded_db)
        families = q.get_all_families()
        assert len(families) == 1
        assert families[0]["accession"] == "PF00394"

    def test_get_all_proteins(self, seeded_db):
        q = QueryEngine(seeded_db)
        proteins = q.get_all_proteins()
        assert len(proteins) == 1
        assert proteins[0]["gene_name"] == "TP53"

    def test_get_proteins_by_family(self, seeded_db):
        q = QueryEngine(seeded_db)
        proteins = q.get_proteins_by_family("PF00394")
        assert len(proteins) == 1

    def test_get_proteins_by_family_missing(self, seeded_db):
        q = QueryEngine(seeded_db)
        assert q.get_proteins_by_family("PF99999") == []

    def test_get_isoforms_for_protein(self, seeded_db):
        q = QueryEngine(seeded_db)
        isoforms = q.get_isoforms_for_protein("P04637")
        assert len(isoforms) == 2
        # Canonical first
        assert isoforms[0]["isoform_id"] == "P04637-1"
        assert isoforms[0]["is_canonical"] == 1

    def test_tim_barrel_location_deserialized(self, seeded_db):
        q = QueryEngine(seeded_db)
        isoforms = q.get_isoforms_for_protein("P04637")
        loc = isoforms[0]["tim_barrel_location"]
        assert isinstance(loc, dict)
        assert loc["start"] == 1

    def test_splice_variants_deserialized(self, seeded_db):
        q = QueryEngine(seeded_db)
        isoforms = q.get_isoforms_for_protein("P04637")
        alt = next(i for i in isoforms if i["isoform_id"] == "P04637-2")
        sv = alt["splice_variants"]
        assert isinstance(sv, list)
        assert sv[0]["featureId"] == "VSP_006535"

    def test_proteins_with_alternative_isoforms(self, seeded_db):
        q = QueryEngine(seeded_db)
        proteins = q.get_proteins_with_alternative_isoforms()
        assert len(proteins) == 1
        assert proteins[0]["isoform_count"] == 1

    def test_isoforms_with_splice_variants(self, seeded_db):
        q = QueryEngine(seeded_db)
        isoforms = q.get_isoforms_with_splice_variants()
        assert len(isoforms) == 1
        assert isoforms[0]["isoform_id"] == "P04637-2"

    def test_summary(self, seeded_db):
        q = QueryEngine(seeded_db)
        s = q.summary()
        assert s["proteins"] == 1
        assert s["isoforms"] == 2
        assert s["alternative_isoforms"] == 1


class TestExport:
    @pytest.fixture
    def isoforms(self, seeded_db):
        return QueryEngine(seeded_db).get_all_isoforms()

    def test_fasta_has_two_sequences(self, isoforms):
        fasta = to_fasta(isoforms)
        assert fasta.count(">") == 2

    def test_fasta_header_format(self, isoforms):
        fasta = to_fasta(isoforms)
        assert ">P04637-1 uniprot_id=P04637 canonical=1" in fasta

    def test_fasta_sequence_wrapped_at_60(self):
        long_iso = [{"isoform_id": "X-1", "uniprot_id": "X",
                     "sequence": "A" * 120, "is_canonical": 1}]
        fasta = to_fasta(long_iso)
        lines = fasta.strip().split("\n")
        seq_lines = [l for l in lines if not l.startswith(">")]
        assert all(len(l) <= 60 for l in seq_lines)

    def test_csv_has_header_and_rows(self, isoforms):
        csv_text = to_csv(isoforms)
        lines = csv_text.strip().split("\n")
        assert "isoform_id" in lines[0]
        assert len(lines) == 3  # header + 2 isoforms

    def test_csv_has_tim_barrel_flag(self, isoforms):
        csv_text = to_csv(isoforms)
        assert "has_tim_barrel" in csv_text

    def test_json_round_trip(self, isoforms):
        import json
        result = json.loads(to_json(isoforms))
        assert len(result) == 2
