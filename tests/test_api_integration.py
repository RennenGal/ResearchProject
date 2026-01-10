"""
Integration tests for REST API endpoints.

Tests end-to-end collection workflow, query functionality, export capabilities,
and error handling through the API.
"""

import pytest
import asyncio
import json
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock
from fastapi.testclient import TestClient

from protein_data_collector.api.rest_api import app, get_query_engine
from protein_data_collector.models.entities import TIMBarrelEntryModel, InterProProteinModel, ProteinModel
from protein_data_collector.collector.data_collector import CollectionReport, CollectionProgress


@pytest.fixture
def client():
    """Create test client for API."""
    return TestClient(app)


@pytest.fixture
def sample_pfam_family():
    """Sample PFAM family for testing."""
    return TIMBarrelEntryModel(
        accession="PF00121",
        entry_type="pfam",
        name="TIM",
        description="Triosephosphate isomerase",
        tim_barrel_annotation="TIM barrel fold"
    )


@pytest.fixture
def sample_interpro_protein():
    """Sample InterPro protein for testing."""
    return InterProProteinModel(
        uniprot_id="P60174",
        pfam_accession="PF00121",
        name="Triosephosphate isomerase",
        organism="Homo sapiens"
    )


@pytest.fixture
def sample_protein_isoform():
    """Sample protein isoform for testing."""
    sequence = "MSKIAKIGEIKDATPSDFVVKASYDLGADVIITGVPVKDGVLYDGKDVTIPAKGDKDVKFVVGVNVLADAVKVTLGPKGRNVVLDKSFGAPTITKDGVSVAREIELEDKFENMGAQMVKEVASKANDAAGDGTTTATVLAQAIITEGLKAVAAGMNPMDLKRGIDKAVTAAVEELKALSVPCSDSKAIAQVGTISANSDETVGKLIAEAMDKVGKEGVITVEDGTGLQDELDVVEGMQFDRGYLSPYFINKPETGAVELESPFILLADKKISNIREMLPVLEAVAKAGKPLLIIAEDVEGEALATLVVNTMRGIVKVAAVKAPGFGDRRKAMLQDIATLTGGTVISEEIGMELEKATLEDLGQAKRVVINKDTTTIIDGVGEEAAIQGRVAQIRQQIEEATSDYDREKLQERVAKLAGGVAVIKVGAATEVEMKEKKARVEDALHATRAAVEEGVVAGGGVALIRVASKLADLRGQNEDQNVGIKVALRAMEAPLRQIVLNCGEEPSVVANTVKFFRRGNRQIKAEEQAERSVAHAAVVAGVQKDVLQVVQKQHPIFEREGNNLYCEVPINFATRQVYSLIRPNENPAHKSQLVWMACHSAAFEDLRVSSFIRGTKVVPRGKLSTRGVQIASNENMETMESSTLELRSRYYDLDVVGDVVCGTGFVDLMVKELQRIRFGGDKSYICATPGYPLVTPVDMSTGEREMTKGVTSADFTNFDPRGLLPESLDYWTYPGSLTTPPLLECVTWIVLKEPISVSSEQVLKFRKLNFNGEGEPEELMVDNWRPAQPLKNRQIKASFK"
    return ProteinModel(
        isoform_id="P60174-1",
        parent_protein_id="P60174",
        sequence=sequence,
        sequence_length=len(sequence),
        exon_annotations={"exons": [{"start": 1, "end": len(sequence)}]},
        exon_count=1,
        tim_barrel_location={"start": 10, "end": 240, "confidence": 0.95},
        organism="Homo sapiens",
        name="Triosephosphate isomerase",
        description="Catalyzes the interconversion of dihydroxyacetone phosphate and D-glyceraldehyde 3-phosphate"
    )


class TestAPIBasics:
    """Test basic API functionality."""
    
    def test_root_endpoint(self, client):
        """Test root endpoint returns API information."""
        response = client.get("/")
        assert response.status_code == 200
        
        data = response.json()
        assert "name" in data
        assert "version" in data
        assert "description" in data
        assert data["name"] == "Protein Data Collector API"
    
    def test_health_check(self, client):
        """Test health check endpoint."""
        with patch('protein_data_collector.database.connection.get_database_manager') as mock_db:
            mock_db.return_value.test_connection.return_value = True
            
            response = client.get("/health")
            assert response.status_code == 200
            
            data = response.json()
            assert "status" in data
            assert data["status"] in ["healthy", "degraded", "unhealthy"]
    
    def test_health_check_db_failure(self, client):
        """Test health check with database failure."""
        with patch('protein_data_collector.api.rest_api.get_database_manager') as mock_db:
            mock_db.return_value.test_connection.return_value = False
            
            response = client.get("/health")
            assert response.status_code == 200
            
            data = response.json()
            assert data["status"] == "degraded"
            assert data["database"] == "disconnected"


class TestSystemStatus:
    """Test system status endpoints."""
    
    def test_system_status(self, client):
        """Test system status endpoint."""
        with patch('protein_data_collector.api.rest_api.get_database_manager') as mock_db, \
             patch('protein_data_collector.api.rest_api.get_health_checker') as mock_health:
            
            mock_db.return_value.test_connection.return_value = True
            
            # Mock the async method
            async def mock_get_data_summary():
                return {
                    "pfam_families": 10,
                    "proteins": 100,
                    "isoforms": 150,
                    "tim_barrel_annotations": 120,
                    "tim_barrel_coverage": 80.0
                }
            
            mock_health.return_value.get_data_summary = mock_get_data_summary
            
            response = client.get("/status")
            assert response.status_code == 200
            
            data = response.json()
            assert "database_connected" in data
            assert "configuration" in data
            assert "data_summary" in data
            assert data["database_connected"] is True
    
    def test_statistics_endpoint(self, client):
        """Test statistics endpoint."""
        mock_stats = {
            "pfam_families": 5,
            "proteins": 50,
            "isoforms": 75,
            "tim_barrel_annotations": 60,
            "tim_barrel_coverage": 80.0,
            "sequence_length": {"min": 100, "max": 500, "avg": 250.0},
            "exon_count": {"min": 1, "max": 10, "avg": 3.5}
        }
        
        # Create mock engine
        mock_engine = MagicMock()
        mock_engine.get_summary_statistics.return_value = mock_stats
        
        # Override the dependency
        app.dependency_overrides[get_query_engine] = lambda: mock_engine
        
        try:
            response = client.get("/query/statistics")
            assert response.status_code == 200
            
            data = response.json()
            assert data == mock_stats
        finally:
            # Clean up dependency override
            app.dependency_overrides.clear()


class TestQueryEndpoints:
    """Test protein query endpoints."""
    
    def test_query_proteins_basic(self, client, sample_protein_isoform):
        """Test basic protein query."""
        from protein_data_collector.query.engine import QueryResult
        
        # Create a proper QueryResult object
        mock_result = QueryResult(
            pfam_families=[],
            proteins=[],
            isoforms=[sample_protein_isoform.model_dump()],
            total_count=1,
            query_metadata={"filters": {}}
        )
        
        # Mock the dependency function to return a mock engine
        def mock_get_query_engine():
            mock_engine = MagicMock()
            mock_engine.filter_proteins.return_value = mock_result
            return mock_engine
        
        # Override the dependency
        app.dependency_overrides[get_query_engine] = mock_get_query_engine
        
        try:
            query_request = {
                "organism": "Homo sapiens",
                "min_sequence_length": 100,
                "limit": 10
            }
            
            response = client.post("/query/proteins", json=query_request)
            assert response.status_code == 200
            
            data = response.json()
            assert "isoforms" in data
            assert "total_count" in data
            assert data["total_count"] == 1
        finally:
            # Clean up dependency override
            app.dependency_overrides.clear()
    
    def test_query_pfam_family(self, client, sample_pfam_family, sample_protein_isoform):
        """Test querying by PFAM family."""
        from protein_data_collector.query.engine import QueryResult
        
        mock_result = QueryResult(
            pfam_families=[sample_pfam_family.model_dump()],
            proteins=[],
            isoforms=[sample_protein_isoform.model_dump()],
            total_count=1,
            query_metadata={"pfam_id": "PF00121", "found": True}
        )
        
        # Create mock engine
        mock_engine = MagicMock()
        mock_engine.search_by_pfam_family.return_value = mock_result
        
        # Override the dependency
        app.dependency_overrides[get_query_engine] = lambda: mock_engine
        
        try:
            response = client.get("/query/pfam/PF00121")
            assert response.status_code == 200
            
            data = response.json()
            assert "pfam_families" in data
            assert "isoforms" in data
            assert data["total_count"] == 1
            assert len(data["pfam_families"]) == 1
        finally:
            # Clean up dependency override
            app.dependency_overrides.clear()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])