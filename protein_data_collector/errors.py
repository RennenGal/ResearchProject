"""Exception hierarchy for the protein data collector."""


class ProteinDataError(Exception):
    """Base exception for all collector errors."""


class APIError(ProteinDataError):
    """HTTP or API-level error from an external service."""

    def __init__(self, message: str, status_code: int = None):
        super().__init__(message)
        self.status_code = status_code


class NetworkError(ProteinDataError):
    """Connection or timeout error."""


class DatabaseError(ProteinDataError):
    """SQLite operation failed."""


class ValidationError(ProteinDataError):
    """Data failed validation (bad sequence, out-of-bounds coordinates, etc.)."""


class DataError(ProteinDataError):
    """Unexpected or missing data in an API response."""
