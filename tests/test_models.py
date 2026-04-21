"""Unit tests for Pydantic data models."""

import pytest
from pydantic import ValidationError

from protein_data_collector.models.entities import Isoform, Protein, TIMBarrelEntry, _FRAGMENT_LENGTH_THRESHOLD


class TestTIMBarrelEntry:
    def test_valid_pfam(self):
        e = TIMBarrelEntry(accession="PF00394", entry_type="pfam",
                           name="TIM", domain_annotation="TIM barrel")
        assert e.accession == "PF00394"

    def test_valid_interpro(self):
        e = TIMBarrelEntry(accession="IPR013785", entry_type="interpro",
                           name="TIM", domain_annotation="TIM barrel")
        assert e.entry_type == "interpro"

    def test_valid_cathgene3d(self):
        e = TIMBarrelEntry(accession="G3DSA:3.20.20.80", entry_type="cathgene3d",
                           name="Glycosidases", domain_annotation="3.20.20")
        assert e.entry_type == "cathgene3d"

    def test_bad_accession(self):
        with pytest.raises(ValidationError, match="Accession must start"):
            TIMBarrelEntry(accession="XY12345", entry_type="pfam",
                           name="X", domain_annotation="X")

    def test_bad_entry_type(self):
        with pytest.raises(ValidationError, match="entry_type must be"):
            TIMBarrelEntry(accession="PF00001", entry_type="unknown",
                           name="X", domain_annotation="X")


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
        assert iso.ensembl_transcript_id is None

    def test_is_fragment_auto_set_for_short_sequence(self):
        short_seq = "ACDEFGHIKL" * 5   # 50 aa — below threshold
        iso = Isoform(isoform_id="P00001-1", uniprot_id="P00001",
                      sequence=short_seq, sequence_length=50)
        assert iso.is_fragment is True

    def test_is_fragment_false_for_full_length(self):
        # _BASE sequence is 20 aa (a test value), so set explicitly above threshold
        long_seq = "ACDEFGHIKLMNPQRSTVWY" * 15  # 300 aa
        iso = Isoform(isoform_id="P00001-1", uniprot_id="P00001",
                      sequence=long_seq, sequence_length=300)
        assert iso.is_fragment is False

    def test_tim_barrel_sequence_auto_sliced(self):
        long_seq = "ACDEFGHIKLMNPQRSTVWY" * 15  # 300 aa
        loc = {"domain_id": "IPR000001", "start": 11, "end": 30, "length": 20, "source": "interpro_api"}
        iso = Isoform(isoform_id="P00001-1", uniprot_id="P00001",
                      sequence=long_seq, sequence_length=300,
                      tim_barrel_location=loc)
        assert iso.tim_barrel_sequence == long_seq[10:30]
        assert len(iso.tim_barrel_sequence) == 20

    def test_tim_barrel_sequence_none_for_fragment(self):
        short_seq = "ACDEFGHIKL" * 5  # 50 aa
        loc = {"domain_id": "IPR000001", "start": 1, "end": 50, "length": 50, "source": "interpro_api"}
        iso = Isoform(isoform_id="P00001-1", uniprot_id="P00001",
                      sequence=short_seq, sequence_length=50,
                      tim_barrel_location=loc)
        assert iso.is_fragment is True
        assert iso.tim_barrel_sequence is None

    def test_fragment_threshold_value(self):
        assert _FRAGMENT_LENGTH_THRESHOLD == 200
