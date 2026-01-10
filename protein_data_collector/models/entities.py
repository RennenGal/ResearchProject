"""
Pydantic models for protein data entities.

This module defines data models for PFAM families, InterPro proteins, and protein isoforms
with comprehensive validation and serialization capabilities.
"""

from datetime import datetime
from typing import Optional, Dict, Any, List
from pydantic import BaseModel, Field, field_validator, model_validator, ConfigDict, model_serializer
import json

from .validation import ProteinSequenceValidator, TIMBarrelLocationValidator


class TIMBarrelEntryModel(BaseModel):
    """Unified Pydantic model for both PFAM families and InterPro entries with TIM barrel annotations."""
    
    accession: str = Field(
        ..., 
        min_length=5, 
        max_length=20,
        description="Entry accession identifier (PFxxxxx or IPRxxxxxx)"
    )
    entry_type: str = Field(
        ...,
        description="Entry type: 'pfam' or 'interpro'"
    )
    name: str = Field(
        ..., 
        min_length=1, 
        max_length=255,
        description="Entry name"
    )
    description: Optional[str] = Field(
        None, 
        max_length=10000,
        description="Entry description"
    )
    interpro_type: Optional[str] = Field(
        None,
        max_length=50,
        description="InterPro entry type (Domain, Family, etc.) - only for IPR entries"
    )
    tim_barrel_annotation: str = Field(
        ..., 
        min_length=1,
        description="TIM barrel structural annotation details"
    )
    member_databases: Dict[str, Any] = Field(
        default_factory=dict,
        description="Member database signatures as dictionary - only for IPR entries"
    )
    interpro_id: Optional[str] = Field(
        None,
        max_length=20,
        description="Associated InterPro identifier - only for PFAM entries"
    )
    created_at: Optional[datetime] = Field(
        default_factory=datetime.now,
        description="Record creation timestamp"
    )
    
    @field_validator('accession')
    @classmethod
    def validate_accession(cls, v):
        """Validate accession format for both PFAM and InterPro."""
        if not v.strip():
            raise ValueError("Accession cannot be empty or whitespace")
        v = v.strip()
        # PFAM accessions start with PF, InterPro with IPR
        if not (v.startswith('PF') or v.startswith('IPR')):
            raise ValueError("Accession must start with 'PF' (PFAM) or 'IPR' (InterPro)")
        return v
    
    @field_validator('entry_type')
    @classmethod
    def validate_entry_type(cls, v):
        """Validate entry type is either pfam or interpro."""
        if v not in ['pfam', 'interpro']:
            raise ValueError("Entry type must be 'pfam' or 'interpro'")
        return v
    
    @field_validator('name', 'description', 'tim_barrel_annotation')
    @classmethod
    def validate_text_fields(cls, v):
        """Validate text fields are not empty or whitespace-only."""
        if v is not None and not v.strip():
            raise ValueError("Text fields cannot be empty or whitespace-only")
        return v.strip() if v else v
    
    @model_validator(mode='after')
    def validate_consistency(self):
        """Validate consistency between accession and entry_type."""
        if self.accession.startswith('PF') and self.entry_type != 'pfam':
            raise ValueError("PFAM accessions (PF*) must have entry_type='pfam'")
        if self.accession.startswith('IPR') and self.entry_type != 'interpro':
            raise ValueError("InterPro accessions (IPR*) must have entry_type='interpro'")
        return self
    
    @property
    def is_pfam(self) -> bool:
        """Check if this is a PFAM entry."""
        return self.entry_type == 'pfam'
    
    @property
    def is_interpro(self) -> bool:
        """Check if this is an InterPro entry."""
        return self.entry_type == 'interpro'
    
    @model_serializer
    def serialize_model(self):
        """Custom serializer to handle datetime objects."""
        data = {
            'accession': self.accession,
            'entry_type': self.entry_type,
            'name': self.name,
            'description': self.description,
            'interpro_type': self.interpro_type,
            'tim_barrel_annotation': self.tim_barrel_annotation,
            'member_databases': self.member_databases,
            'interpro_id': self.interpro_id,
            'created_at': self.created_at.isoformat() if self.created_at else None
        }
        return data
    
    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "accession": "PF00121",
                "entry_type": "pfam",
                "name": "Triosephosphate isomerase",
                "description": "TIM barrel enzyme family",
                "interpro_type": None,
                "tim_barrel_annotation": "Eight-fold alpha/beta barrel structure",
                "member_databases": {},
                "interpro_id": "IPR000652"
            }
        }
    )


class InterProProteinModel(BaseModel):
    """Pydantic model for InterPro protein entities."""
    
    uniprot_id: str = Field(
        ..., 
        min_length=6, 
        max_length=20,
        description="UniProt protein identifier"
    )
    tim_barrel_accession: str = Field(
        ..., 
        min_length=5, 
        max_length=20,
        description="Associated TIM barrel entry accession (PFAM or InterPro)"
    )
    name: Optional[str] = Field(
        None, 
        max_length=255,
        description="Protein name"
    )
    organism: str = Field(
        default="Homo sapiens", 
        max_length=100,
        description="Source organism"
    )
    basic_metadata: Dict[str, Any] = Field(
        default_factory=dict,
        description="Additional protein metadata"
    )
    created_at: Optional[datetime] = Field(
        default_factory=datetime.now,
        description="Record creation timestamp"
    )
    
    @field_validator('uniprot_id')
    @classmethod
    def validate_uniprot_id(cls, v):
        """Validate UniProt ID format."""
        if not v.strip():
            raise ValueError("UniProt ID cannot be empty or whitespace")
        # UniProt IDs are typically 6-10 characters, alphanumeric
        if not v.isalnum():
            raise ValueError("UniProt ID must be alphanumeric")
        return v.strip().upper()
    
    @field_validator('tim_barrel_accession')
    @classmethod
    def validate_tim_barrel_accession(cls, v):
        """Validate TIM barrel entry accession format."""
        if not v.strip():
            raise ValueError("TIM barrel entry accession cannot be empty or whitespace")
        return v.strip()
    
    @field_validator('organism')
    @classmethod
    def validate_organism(cls, v):
        """Validate organism name."""
        if not v.strip():
            raise ValueError("Organism cannot be empty or whitespace")
        return v.strip()
    
    @field_validator('name')
    @classmethod
    def validate_name(cls, v):
        """Validate protein name."""
        if v is not None and not v.strip():
            raise ValueError("Protein name cannot be empty or whitespace-only")
        return v.strip() if v else v
    
    @model_serializer
    def serialize_model(self):
        """Custom serializer to handle datetime objects."""
        data = {
            'uniprot_id': self.uniprot_id,
            'tim_barrel_accession': self.tim_barrel_accession,
            'name': self.name,
            'organism': self.organism,
            'basic_metadata': self.basic_metadata,
            'created_at': self.created_at.isoformat() if self.created_at else None
        }
        return data
    
    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "uniprot_id": "P60174",
                "tim_barrel_accession": "PF00121",
                "name": "Triosephosphate isomerase",
                "organism": "Homo sapiens",
                "basic_metadata": {
                    "gene_name": "TPI1",
                    "protein_length": 249
                }
            }
        }
    )


class ProteinModel(BaseModel):
    """Comprehensive Pydantic model for protein entities with all 67 UniProt fields."""
    
    # Primary identifiers
    isoform_id: str = Field(
        ..., 
        min_length=6, 
        max_length=30,
        description="UniProt isoform identifier"
    )
    parent_protein_id: str = Field(
        ..., 
        min_length=6, 
        max_length=20,
        description="Parent protein UniProt ID"
    )
    parent_tim_barrel_accession: str = Field(
        ...,
        description="Associated TIM barrel entry accession"
    )
    
    # Category field for filtering
    data_category: str = Field(
        default="sequences",
        description="Data category for filtering"
    )
    
    # 1. Names & Taxonomy (14 fields)
    accession: Optional[str] = Field(None, description="UniProt accession")
    entry_name: Optional[str] = Field(None, description="UniProt entry name")
    gene_names: Optional[str] = Field(None, description="Gene names")
    gene_primary: Optional[str] = Field(None, description="Primary gene name")
    gene_synonym: Optional[str] = Field(None, description="Gene synonyms")
    gene_oln: Optional[str] = Field(None, description="Gene OLN")
    gene_orf: Optional[str] = Field(None, description="Gene ORF")
    organism_name: Optional[str] = Field(None, description="Organism name")
    organism_id: Optional[int] = Field(None, description="Organism taxonomy ID")
    protein_name: Optional[str] = Field(None, description="Protein name")
    proteomes: Optional[str] = Field(None, description="Proteomes")
    lineage: Optional[str] = Field(None, description="Taxonomic lineage")
    lineage_ids: Optional[str] = Field(None, description="Taxonomic lineage IDs")
    virus_hosts: Optional[str] = Field(None, description="Virus hosts")
    
    # 2. Sequences (19 fields)
    alternative_products: Optional[str] = Field(None, description="Alternative products")
    alternative_sequence: Optional[str] = Field(None, description="Alternative sequence")
    error_gmodel_pred: Optional[bool] = Field(None, description="Gene model prediction error")
    fragment: Optional[bool] = Field(None, description="Fragment")
    organelle: Optional[str] = Field(None, description="Organelle")
    sequence: str = Field(..., min_length=10, description="Protein amino acid sequence")
    sequence_length: int = Field(..., gt=0, description="Length of protein sequence")
    mass: Optional[float] = Field(None, description="Molecular mass")
    mass_spectrometry: Optional[str] = Field(None, description="Mass spectrometry")
    natural_variant: Optional[str] = Field(None, description="Natural variants")
    non_adjacent_residues: Optional[str] = Field(None, description="Non-adjacent residues")
    non_standard_residue: Optional[str] = Field(None, description="Non-standard residues")
    non_terminal_residue: Optional[str] = Field(None, description="Non-terminal residues")
    polymorphism: Optional[str] = Field(None, description="Polymorphism")
    rna_editing: Optional[str] = Field(None, description="RNA editing")
    sequence_caution: Optional[str] = Field(None, description="Sequence caution")
    sequence_conflict: Optional[str] = Field(None, description="Sequence conflict")
    sequence_uncertainty: Optional[str] = Field(None, description="Sequence uncertainty")
    sequence_version: Optional[int] = Field(None, description="Sequence version")
    
    # 3. Function (16 fields)
    absorption: Optional[str] = Field(None, description="Absorption")
    active_site: Optional[str] = Field(None, description="Active site")
    activity_regulation: Optional[str] = Field(None, description="Activity regulation")
    binding_site: Optional[str] = Field(None, description="Binding site")
    catalytic_activity: Optional[str] = Field(None, description="Catalytic activity")
    cofactor: Optional[str] = Field(None, description="Cofactor")
    dna_binding: Optional[str] = Field(None, description="DNA binding")
    ec_number: Optional[str] = Field(None, description="EC number")
    function_cc: Optional[str] = Field(None, description="Function")
    kinetics: Optional[str] = Field(None, description="Kinetics")
    pathway: Optional[str] = Field(None, description="Pathway")
    ph_dependence: Optional[str] = Field(None, description="pH dependence")
    redox_potential: Optional[str] = Field(None, description="Redox potential")
    rhea_id: Optional[str] = Field(None, description="Rhea ID")
    site: Optional[str] = Field(None, description="Site")
    temp_dependence: Optional[str] = Field(None, description="Temperature dependence")
    
    # 4. Interaction (2 fields)
    interacts_with: Optional[str] = Field(None, description="Interacts with")
    subunit_structure: Optional[str] = Field(None, description="Subunit structure")
    
    # 5. Gene Ontology (5 fields)
    go_biological_process: Optional[str] = Field(None, description="GO biological process")
    go_cellular_component: Optional[str] = Field(None, description="GO cellular component")
    go_molecular_function: Optional[str] = Field(None, description="GO molecular function")
    go_terms: Optional[str] = Field(None, description="GO terms")
    go_ids: Optional[str] = Field(None, description="GO IDs")
    
    # 6. Structure (4 fields)
    structure_3d: Optional[str] = Field(None, description="3D structure")
    beta_strand: Optional[str] = Field(None, description="Beta strand")
    helix: Optional[str] = Field(None, description="Helix")
    turn: Optional[str] = Field(None, description="Turn")
    
    # 7. Date Information (4 fields)
    date_created: Optional[datetime] = Field(None, description="Date created")
    date_modified: Optional[datetime] = Field(None, description="Date modified")
    date_sequence_modified: Optional[datetime] = Field(None, description="Date sequence modified")
    entry_version: Optional[int] = Field(None, description="Entry version")
    
    # 8. Family & Domains (9 fields)
    coiled_coil: Optional[str] = Field(None, description="Coiled coil")
    compositional_bias: Optional[str] = Field(None, description="Compositional bias")
    domain_cc: Optional[str] = Field(None, description="Domain")
    domain_ft: Optional[str] = Field(None, description="Domain feature")
    motif: Optional[str] = Field(None, description="Motif")
    protein_families: Optional[str] = Field(None, description="Protein families")
    region: Optional[str] = Field(None, description="Region")
    repeat_region: Optional[str] = Field(None, description="Repeat region")
    zinc_finger: Optional[str] = Field(None, description="Zinc finger")
    
    # 9. 3D Structure Databases (7 fields)
    xref_alphafolddb: Optional[str] = Field(None, description="AlphaFoldDB cross-reference")
    xref_bmrb: Optional[str] = Field(None, description="BMRB cross-reference")
    xref_pcddb: Optional[str] = Field(None, description="PCDDB cross-reference")
    xref_pdb: Optional[str] = Field(None, description="PDB cross-reference")
    xref_pdbsum: Optional[str] = Field(None, description="PDBsum cross-reference")
    xref_sasbdb: Optional[str] = Field(None, description="SASBDB cross-reference")
    xref_smr: Optional[str] = Field(None, description="SMR cross-reference")
    
    # Custom fields for TIM barrel research
    exon_annotations: Dict[str, Any] = Field(
        default_factory=dict,
        description="Exon annotation data as dictionary"
    )
    exon_count: Optional[int] = Field(
        None, 
        ge=0,
        description="Number of exons calculated from annotations"
    )
    tim_barrel_location: Optional[Dict[str, Any]] = Field(
        None,
        description="TIM barrel location coordinates as dictionary"
    )
    
    # Quality indicators
    protein_existence: Optional[str] = Field(None, description="Protein existence")
    reviewed: Optional[bool] = Field(None, description="Reviewed status")
    annotation_score: Optional[int] = Field(None, description="Annotation score")
    
    # Metadata
    created_at: Optional[datetime] = Field(
        default_factory=datetime.now,
        description="Record creation timestamp"
    )
    updated_at: Optional[datetime] = Field(
        default_factory=datetime.now,
        description="Record update timestamp"
    )
    
    @field_validator('sequence')
    @classmethod
    def validate_sequence(cls, v):
        """Validate protein sequence contains only valid amino acids."""
        validator = ProteinSequenceValidator(allow_extended=True)
        result = validator.validate(v)
        if not result.is_valid:
            raise ValueError(f"Invalid protein sequence: {result.error_message}")
        return v.upper().strip()
    
    @field_validator('sequence_length')
    @classmethod
    def validate_sequence_length_matches(cls, v, info):
        """Validate sequence_length matches actual sequence length."""
        if 'sequence' in info.data and info.data['sequence']:
            actual_length = len(info.data['sequence'].strip())
            if v != actual_length:
                raise ValueError(f"Sequence length {v} does not match actual sequence length {actual_length}")
        return v
    
    @field_validator('tim_barrel_location')
    @classmethod
    def validate_tim_barrel_location(cls, v, info):
        """Validate TIM barrel location coordinates."""
        if v is not None and v and 'sequence_length' in info.data:  # Only validate non-empty locations
            validator = TIMBarrelLocationValidator()
            result = validator.validate(v, info.data['sequence_length'])
            if not result.is_valid:
                raise ValueError(f"Invalid TIM barrel location: {result.error_message}")
        return v
    
    @field_validator('exon_count')
    @classmethod
    def validate_exon_count_matches_annotations(cls, v, info):
        """Validate exon count matches exon annotations."""
        if v is not None and 'exon_annotations' in info.data:
            exon_data = info.data['exon_annotations']
            if isinstance(exon_data, dict) and 'exons' in exon_data:
                actual_count = len(exon_data['exons'])
                if v != actual_count:
                    raise ValueError(f"Exon count {v} does not match annotations count {actual_count}")
        return v
    
    @model_validator(mode='after')
    def validate_consistency(self):
        """Validate overall data consistency."""
        # Ensure isoform_id is related to parent_protein_id
        isoform_id = self.isoform_id or ''
        parent_id = self.parent_protein_id or ''
        
        if isoform_id and parent_id:
            # Isoform IDs often contain the parent protein ID
            if not (parent_id in isoform_id or isoform_id.startswith(parent_id)):
                # This is a warning, not an error, as formats can vary
                pass
        
        return self
    
    def to_json_dict(self) -> Dict[str, Any]:
        """Convert to dictionary suitable for JSON storage."""
        data = self.model_dump()
        # Convert datetime to ISO string
        for field in ['created_at', 'updated_at', 'date_created', 'date_modified', 'date_sequence_modified']:
            if data.get(field):
                data[field] = data[field].isoformat()
        return data
    
    @classmethod
    def from_json_dict(cls, data: Dict[str, Any]) -> 'ProteinModel':
        """Create instance from JSON dictionary."""
        # Convert ISO strings back to datetime
        for field in ['created_at', 'updated_at', 'date_created', 'date_modified', 'date_sequence_modified']:
            if field in data and isinstance(data[field], str):
                try:
                    data[field] = datetime.fromisoformat(data[field])
                except ValueError:
                    data[field] = None
        return cls(**data)
    
    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "isoform_id": "P60174-1",
                "parent_protein_id": "P60174",
                "parent_tim_barrel_accession": "IPR000652",
                "data_category": "sequences",
                "accession": "P60174",
                "entry_name": "TPIS_HUMAN",
                "protein_name": "Triosephosphate isomerase",
                "gene_primary": "TPI1",
                "organism_name": "Homo sapiens",
                "sequence": "MAPSRKFFVGGNWKMNGRKQSLGELIGTLNAAKVPADTEVVCAPPTAYIDFARQKLDPKIAVAAQNCYKVTNGAFTGEISPGMIKDCGATWVVLGHSERRHVFGESDELIGQKVAHALAEGLGVIACIGEKLDEREAGITEKVVFEQTKAIADNVKDWSKVVLAYEPVWAIGTGKTATPQQAQEVHEKLRGWLKSNVSDAVAQSTRIIYGGSVTGATCKELASQPDVDGFLVGGASLKPEFVDIINAKQ",
                "sequence_length": 249,
                "mass": 26669.0,
                "ec_number": "5.3.1.1",
                "function_cc": "Catalyzes the interconversion of dihydroxyacetone phosphate and D-glyceraldehyde 3-phosphate",
                "xref_pdb": "1HTI;1N55;1TPH;1TRI;2JK2;2VEI;4POC;5TRI;6TIM;7TIM",
                "xref_alphafolddb": "P60174",
                "protein_existence": "1: Evidence at protein level",
                "reviewed": True,
                "annotation_score": 5
            }
        }
    )