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
    """Pydantic model for protein isoform entities with detailed annotations."""
    
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
    sequence: str = Field(
        ..., 
        min_length=10,
        description="Protein amino acid sequence"
    )
    sequence_length: int = Field(
        ..., 
        gt=0,
        description="Length of protein sequence"
    )
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
    organism: Optional[str] = Field(
        None, 
        max_length=100,
        description="Source organism"
    )
    name: Optional[str] = Field(
        None, 
        max_length=255,
        description="Protein name"
    )
    description: Optional[str] = Field(
        None,
        description="Protein description"
    )
    created_at: Optional[datetime] = Field(
        default_factory=datetime.now,
        description="Record creation timestamp"
    )
    
    @field_validator('isoform_id')
    @classmethod
    def validate_isoform_id(cls, v):
        """Validate isoform ID format."""
        if not v.strip():
            raise ValueError("Isoform ID cannot be empty or whitespace")
        return v.strip()
    
    @field_validator('parent_protein_id')
    @classmethod
    def validate_parent_protein_id(cls, v):
        """Validate parent protein ID format."""
        if not v.strip():
            raise ValueError("Parent protein ID cannot be empty or whitespace")
        if not v.isalnum():
            raise ValueError("Parent protein ID must be alphanumeric")
        return v.strip().upper()
    
    @field_validator('sequence')
    @classmethod
    def validate_sequence(cls, v):
        """Validate protein sequence contains only valid amino acids."""
        validator = ProteinSequenceValidator()
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
        if v is not None and 'sequence_length' in info.data:
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
    
    @field_validator('name', 'description', 'organism')
    @classmethod
    def validate_text_fields(cls, v):
        """Validate text fields are not empty or whitespace-only."""
        if v is not None and not v.strip():
            raise ValueError("Text fields cannot be empty or whitespace-only")
        return v.strip() if v else v
    
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
        data = self.dict()
        # Convert datetime to ISO string
        if data.get('created_at'):
            data['created_at'] = data['created_at'].isoformat()
        return data
    
    @classmethod
    def from_json_dict(cls, data: Dict[str, Any]) -> 'ProteinModel':
        """Create instance from JSON dictionary."""
        if 'created_at' in data and isinstance(data['created_at'], str):
            data['created_at'] = datetime.fromisoformat(data['created_at'])
        return cls(**data)
    
    @model_serializer
    def serialize_model(self):
        """Custom serializer to handle datetime objects."""
        data = {
            'isoform_id': self.isoform_id,
            'parent_protein_id': self.parent_protein_id,
            'sequence': self.sequence,
            'sequence_length': self.sequence_length,
            'exon_annotations': self.exon_annotations,
            'exon_count': self.exon_count,
            'tim_barrel_location': self.tim_barrel_location,
            'organism': self.organism,
            'name': self.name,
            'description': self.description,
            'created_at': self.created_at.isoformat() if self.created_at else None
        }
        return data
    
    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "isoform_id": "P60174-1",
                "parent_protein_id": "P60174",
                "sequence": "MSKIAKIGEHTPSALAIMENANVLARYASICQQNGIVPIVEPEILPDGDHDLKRCQYVTEKVLAAVYKALSDHHIYLEGTLLKPNMVTPGHACTQKFSHEEIAMATVTALRRTVPPAVTGITFLSGGQSEEEASINLNAINKCPLLKPWALTFSYGRALQASALKAWGGKKENLKAAQEEYVKRALANSLACQGKYTPSGQAGAAASESLFVSNHAY",
                "sequence_length": 249,
                "exon_annotations": {
                    "exons": [
                        {"start": 1, "end": 50},
                        {"start": 51, "end": 150},
                        {"start": 151, "end": 249}
                    ]
                },
                "exon_count": 3,
                "tim_barrel_location": {
                    "start": 10,
                    "end": 240,
                    "confidence": 0.95
                },
                "organism": "Homo sapiens",
                "name": "Triosephosphate isomerase",
                "description": "Catalyzes the interconversion of dihydroxyacetone phosphate and D-glyceraldehyde 3-phosphate"
            }
        }
    )