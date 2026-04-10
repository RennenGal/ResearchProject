"""Pydantic models for the three-tier protein database."""

from datetime import datetime
from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field, field_validator, model_validator


# Standard 20 + IUPAC ambiguity codes used by UniProt (B, Z, X, U, O, J)
_VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBZXUOJ")


class TIMBarrelEntry(BaseModel):
    accession: str
    entry_type: str  # 'pfam' | 'interpro'
    name: str
    description: Optional[str] = None
    tim_barrel_annotation: str
    created_at: Optional[datetime] = None

    @field_validator("accession")
    @classmethod
    def validate_accession(cls, v: str) -> str:
        v = v.strip()
        if not (v.startswith("PF") or v.startswith("IPR")):
            raise ValueError("Accession must start with 'PF' or 'IPR'")
        return v

    @field_validator("entry_type")
    @classmethod
    def validate_entry_type(cls, v: str) -> str:
        if v not in ("pfam", "interpro"):
            raise ValueError("entry_type must be 'pfam' or 'interpro'")
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
    created_at: Optional[datetime] = None


class Isoform(BaseModel):
    isoform_id: str                          # e.g. P04637-1
    uniprot_id: str
    is_canonical: bool = False
    sequence: str
    sequence_length: int
    exon_count: Optional[int] = None
    exon_annotations: Optional[List[Dict[str, Any]]] = None   # [{start, end}, ...]
    splice_variants: Optional[List[Dict[str, Any]]] = None    # UniProt Alternative-sequence features
    tim_barrel_location: Optional[Dict[str, Any]] = None      # {start, end, source}
    ensembl_gene_id: Optional[str] = None
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
