"""
Database schema definitions for the Protein Data Collector system.

This module defines SQLAlchemy models for the three-tier hierarchy:
- PFAM families with TIM barrel annotations
- Human proteins belonging to PFAM families
- Protein isoforms with detailed sequence and structural data
"""

from datetime import datetime
from typing import Optional, List
from sqlalchemy import (
    Column, String, Text, Integer, DateTime, ForeignKey, JSON, Index, ForeignKeyConstraint
)
from sqlalchemy.orm import declarative_base, relationship, Mapped
from sqlalchemy.sql import func

Base = declarative_base()


class TIMBarrelEntry(Base):
    """Unified table for both PFAM families and InterPro entries with TIM barrel annotations."""
    
    __tablename__ = "tim_barrel_entries"
    
    accession = Column(String(20), primary_key=True, doc="Entry accession identifier (PFxxxxx or IPRxxxxxx)")
    entry_type = Column(String(20), nullable=False, doc="Entry type: 'pfam' or 'interpro'")
    name = Column(String(255), nullable=False, doc="Entry name")
    description = Column(Text, doc="Entry description")
    interpro_type = Column(String(50), doc="InterPro entry type (Domain, Family, etc.) - only for IPR entries")
    tim_barrel_annotation = Column(Text, nullable=False, doc="TIM barrel structural annotation details")
    member_databases = Column(JSON, doc="Member database signatures as JSON - only for IPR entries")
    interpro_id = Column(String(20), doc="Associated InterPro identifier - only for PFAM entries")
    created_at = Column(DateTime, default=func.now(), doc="Record creation timestamp")
    
    # Relationship to proteins (updated to use new table name)
    proteins: Mapped[List["InterProProtein"]] = relationship(
        "InterProProtein", back_populates="tim_barrel_entry", cascade="all, delete-orphan"
    )
    
    def __repr__(self) -> str:
        return f"<TIMBarrelEntry(accession='{self.accession}', type='{self.entry_type}', name='{self.name}')>"
    
    @property
    def is_pfam(self) -> bool:
        """Check if this is a PFAM entry."""
        return self.entry_type == 'pfam'
    
    @property
    def is_interpro(self) -> bool:
        """Check if this is an InterPro entry."""
        return self.entry_type == 'interpro'


class InterProProtein(Base):
    """Human proteins belonging to TIM barrel entries (from InterPro)."""
    
    __tablename__ = "interpro_proteins"
    
    uniprot_id = Column(String(20), primary_key=True, doc="UniProt protein identifier")
    tim_barrel_accession = Column(
        String(20), 
        ForeignKey("tim_barrel_entries.accession", ondelete="CASCADE"), 
        primary_key=True,
        doc="Associated TIM barrel entry accession (PFAM or InterPro)"
    )
    name = Column(String(255), doc="Protein name")
    organism = Column(String(100), default="Homo sapiens", doc="Source organism")
    created_at = Column(DateTime, default=func.now(), doc="Record creation timestamp")
    
    # Relationships
    tim_barrel_entry: Mapped["TIMBarrelEntry"] = relationship("TIMBarrelEntry", back_populates="proteins")
    isoforms: Mapped[List["Protein"]] = relationship(
        "Protein", back_populates="parent_protein", cascade="all, delete-orphan"
    )
    
    def __repr__(self) -> str:
        return f"<InterProProtein(uniprot_id='{self.uniprot_id}', tim_barrel='{self.tim_barrel_accession}')>"


class Protein(Base):
    """Protein isoforms with detailed annotations (from UniProt)."""
    
    __tablename__ = "proteins"
    
    isoform_id = Column(String(30), primary_key=True, doc="UniProt isoform identifier")
    parent_protein_id = Column(
        String(20), 
        nullable=False,
        doc="Parent protein UniProt ID"
    )
    parent_tim_barrel_accession = Column(
        String(20),
        nullable=False,
        doc="Parent protein TIM barrel accession"
    )
    sequence = Column(Text, nullable=False, doc="Protein amino acid sequence")
    sequence_length = Column(Integer, nullable=False, doc="Length of protein sequence")
    exon_annotations = Column(JSON, doc="Exon annotation data as JSON")
    exon_count = Column(Integer, doc="Number of exons calculated from annotations")
    tim_barrel_location = Column(JSON, doc="TIM barrel location coordinates as JSON")
    organism = Column(String(100), doc="Source organism")
    name = Column(String(255), doc="Protein name")
    description = Column(Text, doc="Protein description")
    created_at = Column(DateTime, default=func.now(), doc="Record creation timestamp")
    
    # Composite foreign key constraint
    __table_args__ = (
        ForeignKeyConstraint(
            ["parent_protein_id", "parent_tim_barrel_accession"],
            ["interpro_proteins.uniprot_id", "interpro_proteins.tim_barrel_accession"],
            ondelete="CASCADE"
        ),
    )
    
    # Relationship
    parent_protein: Mapped["InterProProtein"] = relationship(
        "InterProProtein", 
        back_populates="isoforms"
    )
    
    def __repr__(self) -> str:
        return f"<Protein(isoform_id='{self.isoform_id}', parent='{self.parent_protein_id}')>"


# Create indexes for efficient querying
Index("idx_tim_barrel_entries_type", TIMBarrelEntry.entry_type)
Index("idx_interpro_proteins_tim_barrel", InterProProtein.tim_barrel_accession)
Index("idx_interpro_proteins_uniprot", InterProProtein.uniprot_id)
Index("idx_proteins_parent", Protein.parent_protein_id)
Index("idx_proteins_parent_composite", Protein.parent_protein_id, Protein.parent_tim_barrel_accession)
Index("idx_interpro_proteins_organism", InterProProtein.organism)
Index("idx_proteins_sequence_length", Protein.sequence_length)
Index("idx_proteins_exon_count", Protein.exon_count)