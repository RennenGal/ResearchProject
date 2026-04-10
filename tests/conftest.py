"""Shared fixtures for the test suite."""

import json
import pytest

from protein_data_collector.models.entities import Isoform, Protein, TIMBarrelEntry


# ---------------------------------------------------------------------------
# Minimal sample UniProt API response (structure mirrors the real API)
# ---------------------------------------------------------------------------

UNIPROT_P53_RESPONSE = {
    "primaryAccession": "P04637",
    "uniProtkbId": "P53_HUMAN",
    "entryType": "UniProtKB reviewed (Swiss-Prot)",
    "annotationScore": 5,
    "proteinExistence": "1: Evidence at protein level",
    "organism": {"scientificName": "Homo sapiens", "taxonId": 9606},
    "proteinDescription": {
        "recommendedName": {
            "fullName": {"value": "Cellular tumor antigen p53"}
        }
    },
    "genes": [{"geneName": {"value": "TP53"}}],
    "reviewed": True,
    "sequence": {
        "value": "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDP",
        "length": 59,
    },
    "comments": [
        {
            "commentType": "ALTERNATIVE PRODUCTS",
            "isoforms": [
                {
                    "isoformIds": ["P04637-1"],
                    "name": {"value": "1"},
                    "isoformSequenceStatus": "Displayed",
                    "sequenceIds": [],
                },
                {
                    "isoformIds": ["P04637-2"],
                    "name": {"value": "2"},
                    "isoformSequenceStatus": "Described",
                    "sequenceIds": ["VSP_006535"],
                },
            ],
        }
    ],
    "features": [
        {
            "type": "Alternative sequence",
            "featureId": "VSP_006535",
            "location": {
                "start": {"value": 1, "modifier": "EXACT"},
                "end": {"value": 10, "modifier": "EXACT"},
            },
            "description": "MEEPQSDPSV -> MK in isoform 2",
        }
    ],
    "uniProtKBCrossReferences": [
        {
            "database": "Ensembl",
            "id": "ENSG00000141510",
            "properties": [
                {"key": "ProteinId", "value": "ENSP00000269305"},
                {"key": "TranscriptId", "value": "ENST00000269305"},
            ],
        },
        {
            "database": "AlphaFoldDB",
            "id": "P04637",
        },
    ],
}


# ---------------------------------------------------------------------------
# Model fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def tim_barrel_entry():
    return TIMBarrelEntry(
        accession="PF00394",
        entry_type="pfam",
        name="TIM barrel",
        tim_barrel_annotation="TIM barrel",
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
    )


@pytest.fixture
def canonical_isoform():
    return Isoform(
        isoform_id="P04637-1",
        uniprot_id="P04637",
        is_canonical=True,
        sequence="MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDP",
        sequence_length=59,
    )
