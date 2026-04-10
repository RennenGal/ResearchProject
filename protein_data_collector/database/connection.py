"""SQLite connection management."""

import sqlite3
from contextlib import contextmanager
from pathlib import Path
from typing import Generator

from .schema import init_db
from ..config import get_config


@contextmanager
def get_connection(db_path: str = None) -> Generator[sqlite3.Connection, None, None]:
    """Yield a sqlite3 connection and close it on exit."""
    path = db_path or get_config().db_path
    conn = sqlite3.connect(path)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    try:
        yield conn
    finally:
        conn.close()


def ensure_db(db_path: str = None) -> None:
    """Create the database file and initialize schema if needed."""
    path = db_path or get_config().db_path
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with get_connection(path) as conn:
        init_db(conn)
