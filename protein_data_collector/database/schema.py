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
    canonical_uniprot_id TEXT,          -- NULL = this IS the canonical; non-NULL = redundant, points to canonical
    created_at           DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (tim_barrel_accession)
        REFERENCES tim_barrel_entries(accession) ON DELETE CASCADE,
    FOREIGN KEY (canonical_uniprot_id)
        REFERENCES proteins(uniprot_id)
);

-- Archive: all collected isoforms including redundant entries from duplicate proteins.
-- Redundant proteins are those where proteins.canonical_uniprot_id IS NOT NULL.
CREATE TABLE IF NOT EXISTS isoforms_with_duplicates (
    isoform_id           TEXT PRIMARY KEY,   -- e.g. P04637-1
    uniprot_id           TEXT NOT NULL,
    is_canonical         INTEGER NOT NULL DEFAULT 0,
    sequence             TEXT NOT NULL,
    sequence_length      INTEGER NOT NULL,
    is_fragment          INTEGER NOT NULL DEFAULT 0,
    exon_count           INTEGER,
    exon_annotations     TEXT,   -- JSON: [{start, end}, ...] in protein coordinates
    splice_variants      TEXT,   -- JSON: UniProt Alternative-sequence features for this isoform
    tim_barrel_location  TEXT,   -- JSON: {domain_id, start, end, length, source}
    tim_barrel_sequence  TEXT,   -- subsequence sequence[start-1:end]; NULL if no location or is_fragment
    ensembl_gene_id      TEXT,
    alphafold_id         TEXT,
    created_at           DATETIME DEFAULT CURRENT_TIMESTAMP
);

-- Primary working table: isoforms belonging only to canonical (non-redundant) proteins.
-- Populated from isoforms_with_duplicates WHERE proteins.canonical_uniprot_id IS NULL.
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

CREATE INDEX IF NOT EXISTS idx_proteins_tim_barrel      ON proteins(tim_barrel_accession);
CREATE INDEX IF NOT EXISTS idx_proteins_canonical_id    ON proteins(canonical_uniprot_id);
CREATE INDEX IF NOT EXISTS idx_isoforms_uniprot         ON isoforms(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_isoforms_canonical       ON isoforms(uniprot_id, is_canonical);
CREATE INDEX IF NOT EXISTS idx_isoforms_length          ON isoforms(sequence_length);

-- Analysis table: alternative isoforms where the TIM barrel is partially affected by AS.
-- Populated by scripts/build_tim_barrel_isoforms.py via sliding-window ungapped alignment
-- of the canonical TIM barrel sequence against each alternative isoform.
-- Included:  12.5% <= identity < 95%   (at least one beta-alpha motif present, meaningful AS effect)
-- Excluded:  identity >= 95%  (TIM barrel effectively unchanged by AS)
-- Excluded:  identity < 12.5%  (< 1 beta-alpha motif, effectively gone)
CREATE TABLE IF NOT EXISTS tim_barrel_isoforms (
    isoform_id                     TEXT PRIMARY KEY,
    uniprot_id                     TEXT NOT NULL,
    is_canonical                   INTEGER NOT NULL DEFAULT 0,
    sequence                       TEXT NOT NULL,
    sequence_length                INTEGER NOT NULL,
    is_fragment                    INTEGER NOT NULL DEFAULT 0,
    exon_count                     INTEGER,
    exon_annotations               TEXT,   -- JSON
    splice_variants                TEXT,   -- JSON
    -- Alignment-derived TIM barrel position in this isoform
    tim_barrel_location            TEXT,   -- JSON: {start, end, length, source:"local_alignment"}
    tim_barrel_sequence            TEXT,   -- isoform[start-1:end] at alignment position
    -- Canonical reference used for alignment
    canonical_tim_barrel_location  TEXT,   -- JSON: original canonical location
    canonical_tim_barrel_sequence  TEXT,   -- canonical TIM barrel sequence (the query)
    -- Alignment results
    identity_percentage            REAL NOT NULL,  -- alignment_score / tim_barrel_length * 100
    alignment_score                INTEGER NOT NULL,
    ensembl_gene_id                TEXT,
    alphafold_id                   TEXT,
    created_at                     DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (uniprot_id) REFERENCES proteins(uniprot_id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_tb_isoforms_uniprot     ON tim_barrel_isoforms(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_tb_isoforms_identity    ON tim_barrel_isoforms(identity_percentage);

-- Guard: reject any isoform insert whose protein has been marked redundant.
CREATE TRIGGER IF NOT EXISTS trg_block_redundant_isoform
BEFORE INSERT ON isoforms
BEGIN
    SELECT CASE
        WHEN (
            SELECT canonical_uniprot_id
            FROM proteins
            WHERE uniprot_id = NEW.uniprot_id
        ) IS NOT NULL
        THEN RAISE(ABORT, 'Isoform rejected: protein is redundant — insert under its canonical_uniprot_id instead')
    END;
END;
"""


def init_db(conn: sqlite3.Connection) -> None:
    """Create tables and indexes if they do not already exist."""
    conn.executescript(_CREATE_TABLES)
    conn.commit()
