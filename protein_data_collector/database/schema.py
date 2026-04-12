"""Database schema — SQL DDL for the three-tier hierarchy."""

import sqlite3


_CREATE_TABLES = """
CREATE TABLE IF NOT EXISTS tim_barrel_entries (
    accession            TEXT PRIMARY KEY,
    entry_type           TEXT NOT NULL CHECK (entry_type IN ('pfam', 'interpro')),
    name                 TEXT NOT NULL,
    description          TEXT,
    tim_barrel_annotation TEXT NOT NULL,
    created_at           DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS proteins (
    uniprot_id           TEXT PRIMARY KEY,
    tim_barrel_accession TEXT NOT NULL,
    protein_name         TEXT,
    gene_name            TEXT,
    organism             TEXT DEFAULT 'Homo sapiens',
    reviewed             INTEGER,       -- 1 = Swiss-Prot, 0 = TrEMBL
    protein_existence    TEXT,
    annotation_score     INTEGER,
    created_at           DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (tim_barrel_accession)
        REFERENCES tim_barrel_entries(accession) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS isoforms (
    isoform_id           TEXT PRIMARY KEY,   -- e.g. P04637-1
    uniprot_id           TEXT NOT NULL,
    is_canonical         INTEGER NOT NULL DEFAULT 0,
    sequence             TEXT NOT NULL,
    sequence_length      INTEGER NOT NULL,
    is_fragment          INTEGER NOT NULL DEFAULT 0,  -- 1 if sequence_length < 200 (cannot contain full TIM barrel)
    exon_count           INTEGER,
    exon_annotations     TEXT,   -- JSON: [{start, end}, ...] in protein coordinates
    splice_variants      TEXT,   -- JSON: UniProt Alternative-sequence features for this isoform
    tim_barrel_location  TEXT,   -- JSON: {domain_id, start, end, length, source}
    tim_barrel_sequence  TEXT,   -- subsequence sequence[start-1:end] from tim_barrel_location; NULL if no location or is_fragment
    ensembl_gene_id      TEXT,
    alphafold_id         TEXT,
    created_at           DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (uniprot_id)
        REFERENCES proteins(uniprot_id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_proteins_tim_barrel  ON proteins(tim_barrel_accession);
CREATE INDEX IF NOT EXISTS idx_isoforms_uniprot     ON isoforms(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_isoforms_canonical   ON isoforms(uniprot_id, is_canonical);
CREATE INDEX IF NOT EXISTS idx_isoforms_length      ON isoforms(sequence_length);
"""


def init_db(conn: sqlite3.Connection) -> None:
    """Create tables and indexes if they do not already exist."""
    conn.executescript(_CREATE_TABLES)
    conn.commit()
