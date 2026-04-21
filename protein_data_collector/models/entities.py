"""Pydantic models for the three-tier protein database."""

from datetime import datetime
from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field, field_validator, model_validator


# Standard 20 + IUPAC ambiguity codes used by UniProt (B, Z, X, U, O, J)
_VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBZXUOJ")


class TIMBarrelEntry(BaseModel):
    accession: str
    entry_type: str  # 'pfam' | 'interpro' | 'cathgene3d'
    name: str
    description: Optional[str] = None
    domain_annotation: str = ""
    created_at: Optional[datetime] = None

    @field_validator("accession")
    @classmethod
    def validate_accession(cls, v: str) -> str:
        v = v.strip()
        if not (v.startswith("PF") or v.startswith("IPR") or v.startswith("G3DSA:")):
            raise ValueError("Accession must start with 'PF', 'IPR', or 'G3DSA:'")
        return v

    @field_validator("entry_type")
    @classmethod
    def validate_entry_type(cls, v: str) -> str:
        if v not in ("pfam", "interpro", "cathgene3d"):
            raise ValueError("entry_type must be 'pfam', 'interpro', or 'cathgene3d'")
        return v


class Protein(BaseModel):
    uniprot_id: str
    tim_barrel_accession: str
    protein_name: Optional[str] = None
    gene_name: Optional[str] = None
    organism: str = "Homo sapiens"
    reviewed: Optional[bool] = None
    protein_existence: Optional[str] = None
    annotation_score: Optional[int] = None
    canonical_uniprot_id: Optional[str] = None  # None = canonical; set = redundant entry pointing to its canonical
    created_at: Optional[datetime] = None


# Minimum length to contain a full TIM barrel domain (~250 residues typical).
# Sequences below this threshold are UniProt fragments, not full-length proteins.
_FRAGMENT_LENGTH_THRESHOLD = 200


class Isoform(BaseModel):
    isoform_id: str                          # e.g. P04637-1
    uniprot_id: str
    is_canonical: bool = False
    sequence: str
    sequence_length: int
    is_fragment: bool = False                # True if sequence_length < _FRAGMENT_LENGTH_THRESHOLD
    exon_count: Optional[int] = None
    exon_annotations: Optional[List[Dict[str, Any]]] = None   # [{start, end}, ...]
    splice_variants: Optional[List[Dict[str, Any]]] = None    # UniProt Alternative-sequence features
    tim_barrel_location: Optional[Dict[str, Any]] = None      # {domain_id, start, end, length, source}
    tim_barrel_sequence: Optional[str] = None                 # sequence[start-1:end]; None for fragments or missing location
    ensembl_transcript_id: Optional[str] = None
    alphafold_id: Optional[str] = None
    created_at: Optional[datetime] = None

    @field_validator("sequence")
    @classmethod
    def validate_sequence(cls, v: str) -> str:
        invalid = set(v.upper()) - _VALID_AA
        if invalid:
            raise ValueError(f"Invalid amino acid characters: {sorted(invalid)}")
        return v

    @model_validator(mode="after")
    def check_length(self) -> "Isoform":
        if self.sequence and len(self.sequence) != self.sequence_length:
            raise ValueError(
                f"sequence_length {self.sequence_length} != actual length {len(self.sequence)}"
            )
        return self

    @model_validator(mode="after")
    def check_tim_barrel_bounds(self) -> "Isoform":
        loc = self.tim_barrel_location
        if loc and self.sequence_length:
            start, end = loc.get("start", 0), loc.get("end", 0)
            if start > 0 and end > 0:
                if start >= end:
                    raise ValueError(f"TIM barrel start {start} >= end {end}")
                if end > self.sequence_length:
                    raise ValueError(
                        f"TIM barrel end {end} exceeds sequence length {self.sequence_length}"
                    )
        return self

    @model_validator(mode="after")
    def compute_derived_fields(self) -> "Isoform":
        # is_fragment: sequence too short to contain a full TIM barrel domain
        if not self.is_fragment:
            self.is_fragment = self.sequence_length < _FRAGMENT_LENGTH_THRESHOLD

        # tim_barrel_sequence: slice from location if available and not a fragment
        if self.tim_barrel_sequence is None and self.tim_barrel_location and not self.is_fragment:
            start = self.tim_barrel_location.get("start")
            end = self.tim_barrel_location.get("end")
            if start and end and self.sequence:
                self.tim_barrel_sequence = self.sequence[start - 1:end]

        return self
