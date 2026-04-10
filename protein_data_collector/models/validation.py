"""Standalone validation helpers (used outside Pydantic when needed)."""

from typing import Any, Dict, Optional


# Standard 20 + IUPAC ambiguity codes used by UniProt (B, Z, X, U, O, J)
_VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBZXUOJ")


def validate_sequence(sequence: str) -> None:
    """Raise ValueError if sequence contains non-amino-acid characters."""
    invalid = set(sequence.upper()) - _VALID_AA
    if invalid:
        raise ValueError(f"Invalid amino acid characters: {sorted(invalid)}")


def validate_tim_barrel_location(location: Dict[str, Any], sequence_length: int) -> None:
    """Raise ValueError if TIM barrel coordinates are out of bounds or invalid."""
    start = location.get("start", 0)
    end = location.get("end", 0)
    if start <= 0 or end <= 0:
        return  # partial / missing coordinates — not an error at this stage
    if start >= end:
        raise ValueError(f"TIM barrel start {start} >= end {end}")
    if end > sequence_length:
        raise ValueError(
            f"TIM barrel end {end} exceeds sequence length {sequence_length}"
        )
