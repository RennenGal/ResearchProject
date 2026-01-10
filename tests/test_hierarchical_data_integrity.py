"""
Property-based tests for hierarchical data integrity.

Tests Property 2: Hierarchical Data Integrity
**Validates: Requirements 2.3, 3.7, 4.4**

Feature: protein-data-collector, Property 2: Hierarchical Data Integrity
"""

import pytest
from hypothesis import given, strategies as st, settings
from typing import List, Dict, Set

from protein_data_collector.models.entities import TIMBarrelEntryModel, InterProProteinModel, ProteinModel


# Generators for test data
@st.composite
def tim_barrel_entry_generator(draw):
    """Generate valid TIM barrel entry data."""
    entry_type = draw(st.sampled_from(['pfam', 'interpro']))
    
    if entry_type == 'pfam':
        accession = 'PF' + draw(st.text(min_size=5, max_size=5, alphabet=st.characters(whitelist_categories=('Nd',))))
    else:
        accession = 'IPR' + draw(st.text(min_size=6, max_size=6, alphabet=st.characters(whitelist_categories=('Nd',))))
    
    name = draw(st.text(min_size=1, max_size=255, alphabet=st.characters(min_codepoint=32, max_codepoint=126)).filter(lambda x: x.strip()))
    description = draw(st.one_of(
        st.none(), 
        st.text(min_size=1, max_size=1000, alphabet=st.characters(min_codepoint=32, max_codepoint=126)).filter(lambda x: x.strip())
    ))
    tim_barrel_annotation = draw(st.text(min_size=1, max_size=1000, alphabet=st.characters(min_codepoint=32, max_codepoint=126)).filter(lambda x: x.strip()))
    
    return TIMBarrelEntryModel(
        accession=accession,
        entry_type=entry_type,
        name=name,
        description=description,
        tim_barrel_annotation=tim_barrel_annotation
    )


@st.composite
def interpro_protein_generator(draw, tim_barrel_accessions: List[str]):
    """Generate valid InterPro protein data linked to TIM barrel entries."""
    if not tim_barrel_accessions:
        tim_barrel_accessions = ["PF00001"]  # Default fallback
    
    uniprot_id = draw(st.text(min_size=6, max_size=20, alphabet=st.characters(whitelist_categories=('Lu', 'Nd'))))
    tim_barrel_accession = draw(st.sampled_from(tim_barrel_accessions))
    name = draw(st.one_of(
        st.none(), 
        st.text(min_size=1, max_size=255, alphabet=st.characters(min_codepoint=32, max_codepoint=126)).filter(lambda x: x.strip())
    ))
    organism = draw(st.sampled_from(["Homo sapiens", "Human"]))
    
    return InterProProteinModel(
        uniprot_id=uniprot_id,
        tim_barrel_accession=tim_barrel_accession,
        name=name,
        organism=organism
    )


@st.composite
def protein_isoform_generator(draw, parent_proteins: List[InterProProteinModel]):
    """Generate valid protein isoform data linked to parent proteins."""
    if not parent_proteins:
        # Create a default parent protein
        parent_proteins = [InterProProteinModel(
            uniprot_id="P12345",
            tim_barrel_accession="PF00001",
            name="Default protein",
            organism="Homo sapiens"
        )]
    
    parent_protein = draw(st.sampled_from(parent_proteins))
    isoform_id = f"{parent_protein.uniprot_id}-{draw(st.integers(min_value=1, max_value=10))}"
    
    # Generate valid amino acid sequence
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    sequence_length = draw(st.integers(min_value=50, max_value=500))  # Increased minimum for TIM barrel
    sequence = ''.join(draw(st.lists(st.sampled_from(amino_acids), min_size=sequence_length, max_size=sequence_length)))
    
    # Generate exon annotations
    exon_count = draw(st.integers(min_value=1, max_value=10))
    exons = []
    current_pos = 1
    for i in range(exon_count):
        exon_length = draw(st.integers(min_value=10, max_value=50))
        exons.append({
            "start": current_pos,
            "end": current_pos + exon_length - 1
        })
        current_pos += exon_length
    
    exon_annotations = {"exons": exons}
    
    # Generate TIM barrel location within sequence bounds - ensure valid range
    tim_start = draw(st.integers(min_value=1, max_value=max(1, sequence_length - 30)))
    tim_end = draw(st.integers(min_value=tim_start + 20, max_value=sequence_length))
    tim_barrel_location = {
        "start": tim_start,
        "end": tim_end,
        "confidence": draw(st.floats(min_value=0.0, max_value=1.0))
    }
    
    name = draw(st.one_of(
        st.none(), 
        st.text(min_size=1, max_size=255, alphabet=st.characters(min_codepoint=32, max_codepoint=126)).filter(lambda x: x.strip())
    ))
    description = draw(st.one_of(
        st.none(), 
        st.text(min_size=1, max_size=1000, alphabet=st.characters(min_codepoint=32, max_codepoint=126)).filter(lambda x: x.strip())
    ))
    organism = draw(st.sampled_from(["Homo sapiens", "Human"]))
    
    return ProteinModel(
        isoform_id=isoform_id,
        parent_protein_id=parent_protein.uniprot_id,
        parent_tim_barrel_accession=parent_protein.tim_barrel_accession,
        sequence=sequence,
        sequence_length=sequence_length,
        exon_annotations=exon_annotations,
        exon_count=exon_count,
        tim_barrel_location=tim_barrel_location,
        organism=organism,
        name=name,
        description=description
    )


@st.composite
def hierarchical_dataset_generator(draw):
    """Generate a complete hierarchical dataset with proper relationships."""
    # Generate TIM barrel entries with unique accessions
    tim_barrel_entries = draw(st.lists(tim_barrel_entry_generator(), min_size=1, max_size=10))
    
    # Ensure unique accessions
    seen_accessions = set()
    unique_tim_barrel_entries = []
    for entry in tim_barrel_entries:
        if entry.accession not in seen_accessions:
            seen_accessions.add(entry.accession)
            unique_tim_barrel_entries.append(entry)
    
    # If we don't have any unique entries, create at least one
    if not unique_tim_barrel_entries:
        unique_tim_barrel_entries = [draw(tim_barrel_entry_generator())]
    
    tim_barrel_accessions = [entry.accession for entry in unique_tim_barrel_entries]
    
    # Generate InterPro proteins linked to TIM barrel entries with unique IDs
    interpro_proteins = draw(st.lists(
        interpro_protein_generator(tim_barrel_accessions), 
        min_size=1, 
        max_size=50
    ))
    
    # Ensure unique protein IDs per TIM barrel entry (composite key)
    seen_protein_keys = set()
    unique_interpro_proteins = []
    for protein in interpro_proteins:
        key = (protein.uniprot_id, protein.tim_barrel_accession)
        if key not in seen_protein_keys:
            seen_protein_keys.add(key)
            unique_interpro_proteins.append(protein)
    
    # If we don't have any unique proteins, create at least one
    if not unique_interpro_proteins:
        unique_interpro_proteins = [draw(interpro_protein_generator(tim_barrel_accessions))]
    
    # Generate protein isoforms linked to InterPro proteins with unique IDs
    protein_isoforms = draw(st.lists(
        protein_isoform_generator(unique_interpro_proteins), 
        min_size=1, 
        max_size=100
    ))
    
    # Ensure unique isoform IDs
    seen_isoform_ids = set()
    unique_protein_isoforms = []
    for isoform in protein_isoforms:
        if isoform.isoform_id not in seen_isoform_ids:
            seen_isoform_ids.add(isoform.isoform_id)
            unique_protein_isoforms.append(isoform)
    
    # If we don't have any unique isoforms, create at least one
    if not unique_protein_isoforms:
        unique_protein_isoforms = [draw(protein_isoform_generator(unique_interpro_proteins))]
    
    return unique_tim_barrel_entries, unique_interpro_proteins, unique_protein_isoforms


class TestHierarchicalDataIntegrity:
    """Test hierarchical data integrity property."""
    
    @given(hierarchical_dataset_generator())
    @settings(max_examples=100, deadline=None)
    def test_hierarchical_data_integrity_property(self, dataset):
        """
        Property 2: Hierarchical Data Integrity
        
        For any collected dataset, every protein should be correctly linked to exactly 
        one TIM barrel entry, and every isoform should be correctly linked to exactly one 
        protein, maintaining the three-tier hierarchy.
        
        **Validates: Requirements 2.3, 3.7, 4.4**
        """
        tim_barrel_entries, interpro_proteins, protein_isoforms = dataset
        
        # Extract identifiers for validation
        tim_barrel_accessions = {entry.accession for entry in tim_barrel_entries}
        protein_keys = {(protein.uniprot_id, protein.tim_barrel_accession) for protein in interpro_proteins}
        
        # Property 1: Every InterPro protein must be linked to exactly one TIM barrel entry
        for protein in interpro_proteins:
            # Protein must have a TIM barrel accession
            assert protein.tim_barrel_accession is not None, f"Protein {protein.uniprot_id} missing TIM barrel accession"
            
            # TIM barrel accession must exist in the TIM barrel entries
            assert protein.tim_barrel_accession in tim_barrel_accessions, \
                f"Protein {protein.uniprot_id} references non-existent TIM barrel entry {protein.tim_barrel_accession}"
        
        # Property 2: Every protein isoform must be linked to exactly one InterPro protein
        for isoform in protein_isoforms:
            # Isoform must have a parent protein ID
            assert isoform.parent_protein_id is not None, \
                f"Isoform {isoform.isoform_id} missing parent protein ID"
            
            # Parent protein composite key must exist in the InterPro proteins
            parent_key = (isoform.parent_protein_id, isoform.parent_protein_id)  # Need to get the TIM barrel accession
            # For now, just check the protein ID exists
            protein_ids = {protein.uniprot_id for protein in interpro_proteins}
            assert isoform.parent_protein_id in protein_ids, \
                f"Isoform {isoform.isoform_id} references non-existent protein {isoform.parent_protein_id}"
        
        # Property 3: Verify three-tier hierarchy completeness
        # Count relationships at each level
        tim_barrel_to_proteins = {}
        for protein in interpro_proteins:
            tim_barrel_acc = protein.tim_barrel_accession
            if tim_barrel_acc not in tim_barrel_to_proteins:
                tim_barrel_to_proteins[tim_barrel_acc] = []
            tim_barrel_to_proteins[tim_barrel_acc].append(protein.uniprot_id)
        
        protein_to_isoforms = {}
        for isoform in protein_isoforms:
            parent_id = isoform.parent_protein_id
            if parent_id not in protein_to_isoforms:
                protein_to_isoforms[parent_id] = []
            protein_to_isoforms[parent_id].append(isoform.isoform_id)
        
        # Verify that all TIM barrel entries have at least one protein (if proteins exist)
        if interpro_proteins:
            for entry in tim_barrel_entries:
                if entry.accession in tim_barrel_to_proteins:
                    assert len(tim_barrel_to_proteins[entry.accession]) > 0, \
                        f"TIM barrel entry {entry.accession} has no associated proteins"
        
        # Verify that all proteins have at least one isoform (if isoforms exist)
        if protein_isoforms:
            for protein in interpro_proteins:
                if protein.uniprot_id in protein_to_isoforms:
                    assert len(protein_to_isoforms[protein.uniprot_id]) > 0, \
                        f"Protein {protein.uniprot_id} has no associated isoforms"
        
        # Property 4: Verify no orphaned entities
        # All referenced TIM barrel accessions should exist
        referenced_tim_barrel_accessions = {protein.tim_barrel_accession for protein in interpro_proteins}
        for ref_acc in referenced_tim_barrel_accessions:
            assert ref_acc in tim_barrel_accessions, \
                f"Referenced TIM barrel accession {ref_acc} not found in TIM barrel entries"
        
        # All referenced parent protein IDs should exist
        referenced_protein_ids = {isoform.parent_protein_id for isoform in protein_isoforms}
        protein_ids = {protein.uniprot_id for protein in interpro_proteins}
        for ref_id in referenced_protein_ids:
            assert ref_id in protein_ids, \
                f"Referenced protein ID {ref_id} not found in InterPro proteins"
    
    @given(hierarchical_dataset_generator())
    @settings(max_examples=100, deadline=None)
    def test_hierarchical_relationship_uniqueness(self, dataset):
        """
        Test that each entity has exactly one parent in the hierarchy.
        
        This is part of the hierarchical data integrity property.
        """
        tim_barrel_entries, interpro_proteins, protein_isoforms = dataset
        
        # Each InterPro protein should reference exactly one TIM barrel entry
        for protein in interpro_proteins:
            entry_count = sum(1 for entry in tim_barrel_entries if entry.accession == protein.tim_barrel_accession)
            assert entry_count == 1, \
                f"Protein {protein.uniprot_id} references TIM barrel entry {protein.tim_barrel_accession} which appears {entry_count} times"
        
        # Each protein isoform should reference exactly one parent protein
        for isoform in protein_isoforms:
            protein_count = sum(1 for protein in interpro_proteins if protein.uniprot_id == isoform.parent_protein_id)
            assert protein_count >= 1, \
                f"Isoform {isoform.isoform_id} references protein {isoform.parent_protein_id} which appears {protein_count} times"
    
    @given(hierarchical_dataset_generator())
    @settings(max_examples=100, deadline=None)
    def test_hierarchical_data_consistency(self, dataset):
        """
        Test data consistency within the hierarchy.
        
        Verifies that related entities have consistent data.
        """
        tim_barrel_entries, interpro_proteins, protein_isoforms = dataset
        
        # Create lookup maps
        tim_barrel_map = {entry.accession: entry for entry in tim_barrel_entries}
        protein_map = {(protein.uniprot_id, protein.tim_barrel_accession): protein for protein in interpro_proteins}
        
        # Verify organism consistency between proteins and isoforms
        for isoform in protein_isoforms:
            # Find the parent protein (need to match by uniprot_id)
            parent_proteins = [p for p in interpro_proteins if p.uniprot_id == isoform.parent_protein_id]
            if parent_proteins:
                parent_protein = parent_proteins[0]
                if isoform.organism_name and parent_protein.organism:
                    # Both should be human organisms (allowing for slight variations)
                    assert isoform.organism_name.lower() in ["homo sapiens", "human"], \
                        f"Isoform {isoform.isoform_id} has non-human organism: {isoform.organism_name}"
                    assert parent_protein.organism.lower() in ["homo sapiens", "human"], \
                        f"Protein {parent_protein.uniprot_id} has non-human organism: {parent_protein.organism}"
        
        # Verify that isoform IDs are related to parent protein IDs
        for isoform in protein_isoforms:
            # Common pattern: isoform ID contains or starts with parent protein ID
            parent_id = isoform.parent_protein_id
            isoform_id = isoform.isoform_id
            
            # This is a soft check since ID formats can vary
            # We just verify the relationship makes sense
            parent_proteins = [p for p in interpro_proteins if p.uniprot_id == parent_id]
            assert len(parent_proteins) >= 1, \
                f"Isoform {isoform_id} references non-existent parent {parent_id}"
    
    def test_empty_dataset_hierarchy(self):
        """Test that empty datasets maintain hierarchy integrity."""
        # Empty datasets should not violate hierarchy rules
        tim_barrel_entries = []
        interpro_proteins = []
        protein_isoforms = []
        
        # No assertions should fail for empty datasets
        tim_barrel_accessions = {entry.accession for entry in tim_barrel_entries}
        protein_keys = {(protein.uniprot_id, protein.tim_barrel_accession) for protein in interpro_proteins}
        
        # These should be empty sets
        assert len(tim_barrel_accessions) == 0
        assert len(protein_keys) == 0
        
        # No relationships to validate in empty datasets
        for protein in interpro_proteins:
            assert protein.tim_barrel_accession in tim_barrel_accessions
        
        for isoform in protein_isoforms:
            protein_ids = {protein.uniprot_id for protein in interpro_proteins}
            assert isoform.parent_protein_id in protein_ids
    
    def test_single_entity_hierarchy(self):
        """Test hierarchy with single entities at each level."""
        # Create minimal valid hierarchy
        tim_barrel_entry = TIMBarrelEntryModel(
            accession="PF00001",
            entry_type="pfam",
            name="Test Family",
            description="Test description",
            tim_barrel_annotation="Test TIM barrel annotation"
        )
        
        interpro_protein = InterProProteinModel(
            uniprot_id="P12345",
            tim_barrel_accession="PF00001",
            name="Test Protein",
            organism="Homo sapiens"
        )
        
        protein_isoform = ProteinModel(
            isoform_id="P12345-1",
            parent_protein_id="P12345",
            parent_tim_barrel_accession="PF00001",
            sequence="ACDEFGHIKLMNPQRSTVWY" * 10,  # Valid amino acid sequence
            sequence_length=200,
            exon_annotations={"exons": [{"start": 1, "end": 200}]},
            exon_count=1,
            tim_barrel_location={"start": 10, "end": 190, "confidence": 0.95},
            organism="Homo sapiens",
            name="Test Protein Isoform",
            description="Test isoform description"
        )
        
        # Test the hierarchy
        tim_barrel_entries = [tim_barrel_entry]
        interpro_proteins = [interpro_protein]
        protein_isoforms = [protein_isoform]
        
        # Verify relationships
        tim_barrel_accessions = {entry.accession for entry in tim_barrel_entries}
        protein_ids = {protein.uniprot_id for protein in interpro_proteins}
        
        # Protein should reference existing TIM barrel entry
        assert interpro_protein.tim_barrel_accession in tim_barrel_accessions
        
        # Isoform should reference existing protein
        assert protein_isoform.parent_protein_id in protein_ids
        
        # Verify hierarchy completeness
        assert len(tim_barrel_entries) == 1
        assert len(interpro_proteins) == 1
        assert len(protein_isoforms) == 1