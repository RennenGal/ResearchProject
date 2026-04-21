"""
Tests for UniProt response parsing — the research-critical extraction logic.

All tests use the minimal sample response from conftest.py (no network calls).
"""

import pytest

from protein_data_collector.collector.uniprot_collector import (
    _extract_all_splice_features,
    _extract_alphafold_id,
    _extract_ensembl_transcript_id,
    _parse_alternative_products,
)
from tests.conftest import UNIPROT_P53_RESPONSE


class TestParseAlternativeProducts:
    def test_returns_two_isoforms(self):
        isoforms = _parse_alternative_products(UNIPROT_P53_RESPONSE)
        assert len(isoforms) == 2

    def test_canonical_identified(self):
        isoforms = _parse_alternative_products(UNIPROT_P53_RESPONSE)
        canonical = next(i for i in isoforms if i["isoform_id"] == "P04637-1")
        assert canonical["sequence_status"] == "Displayed"

    def test_alternative_has_vsp_id(self):
        isoforms = _parse_alternative_products(UNIPROT_P53_RESPONSE)
        alt = next(i for i in isoforms if i["isoform_id"] == "P04637-2")
        assert "VSP_006535" in alt["sequence_ids"]

    def test_no_alternative_products(self):
        data = {**UNIPROT_P53_RESPONSE, "comments": []}
        assert _parse_alternative_products(data) == []


class TestExtractSpliceFeatures:
    def test_returns_one_feature(self):
        features = _extract_all_splice_features(UNIPROT_P53_RESPONSE)
        assert len(features) == 1

    def test_feature_has_expected_fields(self):
        features = _extract_all_splice_features(UNIPROT_P53_RESPONSE)
        f = features[0]
        assert f["featureId"] == "VSP_006535"
        assert f["location"]["start"]["value"] == 1
        assert f["location"]["end"]["value"] == 10
        assert "isoform 2" in f["description"]

    def test_non_splice_features_excluded(self):
        data = {
            **UNIPROT_P53_RESPONSE,
            "features": [
                {"type": "Domain", "featureId": "DOM_001",
                 "location": {}, "description": "TIM barrel"},
                {"type": "Alternative sequence", "featureId": "VSP_999",
                 "location": {}, "description": "X -> Y in isoform 2"},
            ],
        }
        features = _extract_all_splice_features(data)
        assert len(features) == 1
        assert features[0]["featureId"] == "VSP_999"

    def test_empty_features(self):
        data = {**UNIPROT_P53_RESPONSE, "features": []}
        assert _extract_all_splice_features(data) == []


class TestExtractCrossReferences:
    def test_ensembl_transcript_id(self):
        gene_id = _extract_ensembl_transcript_id(UNIPROT_P53_RESPONSE)
        assert gene_id == "ENSG00000141510"

    def test_ensembl_missing(self):
        data = {**UNIPROT_P53_RESPONSE, "uniProtKBCrossReferences": []}
        assert _extract_ensembl_transcript_id(data) is None

    def test_alphafold_id(self):
        af_id = _extract_alphafold_id(UNIPROT_P53_RESPONSE)
        assert af_id == "P04637"

    def test_alphafold_missing(self):
        data = {**UNIPROT_P53_RESPONSE, "uniProtKBCrossReferences": []}
        assert _extract_alphafold_id(data) is None


class TestSpliceVariantLinkage:
    """Verify that alternative isoforms are linked to the right splice features."""

    def test_alt_isoform_gets_matching_feature(self):
        isoforms = _parse_alternative_products(UNIPROT_P53_RESPONSE)
        features = _extract_all_splice_features(UNIPROT_P53_RESPONSE)

        alt = next(i for i in isoforms if i["isoform_id"] == "P04637-2")
        vsp_ids = set(alt["sequence_ids"])
        linked = [f for f in features if f["featureId"] in vsp_ids]
        assert len(linked) == 1
        assert linked[0]["featureId"] == "VSP_006535"

    def test_canonical_has_no_vsp_ids(self):
        isoforms = _parse_alternative_products(UNIPROT_P53_RESPONSE)
        canonical = next(i for i in isoforms if i["isoform_id"] == "P04637-1")
        assert canonical["sequence_ids"] == []
