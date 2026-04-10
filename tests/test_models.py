"""Unit tests for Pydantic data models."""

import pytest
from pydantic import ValidationError

from protein_data_collector.models.entities import Isoform, Protein, TIMBarrelEntry


class TestTIMBarrelEntry:
    def test_valid_pfam(self):
        e = TIMBarrelEntry(accession="PF00394", entry_type="pfam",
                           name="TIM", tim_barrel_annotation="TIM barrel")
        assert e.accession == "PF00394"

    def test_valid_interpro(self):
        e = TIMBarrelEntry(accession="IPR013785", entry_type="interpro",
                           name="TIM", tim_barrel_annotation="TIM barrel")
        assert e.entry_type == "interpro"

    def test_bad_accession(self):
        with pytest.raises(ValidationError, match="Accession must start"):
            TIMBarrelEntry(accession="XY12345", entry_type="pfam",
                           name="X", tim_barrel_annotation="X")

    def test_bad_entry_type(self):
        with pytest.raises(ValidationError, match="entry_type must be"):
            TIMBarrelEntry(accession="PF00001", entry_type="unknown",
                           name="X", tim_barrel_annotation="X")


class TestIsoform:
    _BASE = dict(
        isoform_id="P04637-1",
        uniprot_id="P04637",
        is_canonical=True,
        sequence="ACDEFGHIKLMNPQRSTVWY",
        sequence_length=20,
    )

    def test_valid(self):
        iso = Isoform(**self._BASE)
        assert iso.sequence_length == 20

    def test_invalid_amino_acid(self):
        # '1' and '?' are never valid in any IUPAC scheme
        with pytest.raises(ValidationError, match="Invalid amino acid"):
            Isoform(**{**self._BASE, "sequence": "ACDEF1??", "sequence_length": 8})

    def test_length_mismatch(self):
        with pytest.raises(ValidationError, match="sequence_length"):
            Isoform(**{**self._BASE, "sequence_length": 999})

    def test_tim_barrel_bounds_ok(self):
        iso = Isoform(**{**self._BASE,
                         "tim_barrel_location": {"start": 1, "end": 15}})
        assert iso.tim_barrel_location["end"] == 15

    def test_tim_barrel_end_exceeds_length(self):
        with pytest.raises(ValidationError, match="exceeds sequence length"):
            Isoform(**{**self._BASE,
                       "tim_barrel_location": {"start": 1, "end": 999}})

    def test_tim_barrel_start_gte_end(self):
        with pytest.raises(ValidationError, match="start .* >= end"):
            Isoform(**{**self._BASE,
                       "tim_barrel_location": {"start": 10, "end": 5}})

    def test_tim_barrel_missing_coords_allowed(self):
        # Partial location with zeros is treated as unknown
        iso = Isoform(**{**self._BASE,
                         "tim_barrel_location": {"start": 0, "end": 0}})
        assert iso.tim_barrel_location is not None

    def test_optional_fields_default_none(self):
        iso = Isoform(**self._BASE)
        assert iso.exon_annotations is None
        assert iso.splice_variants is None
        assert iso.tim_barrel_location is None
        assert iso.ensembl_gene_id is None
