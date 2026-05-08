"""Database schema — SQL DDL for TIM barrel / Homo sapiens tables."""

import sqlite3


_CREATE_TABLES = """
-- ============================================================
-- Domain entries
-- ============================================================

CREATE TABLE IF NOT EXISTS entries (
    accession         TEXT PRIMARY KEY,
    entry_type        TEXT NOT NULL CHECK (entry_type IN ('pfam', 'interpro', 'cathgene3d')),
    name              TEXT NOT NULL,
    description       TEXT,
    domain_annotation TEXT NOT NULL,
    created_at        DATETIME DEFAULT CURRENT_TIMESTAMP
);

-- ============================================================
-- TIM barrel — Homo sapiens
-- ============================================================

CREATE TABLE IF NOT EXISTS proteins (
    uniprot_id           TEXT PRIMARY KEY,
    tim_barrel_accession TEXT NOT NULL,
    protein_name         TEXT,
    gene_name            TEXT,
    organism             TEXT,
    reviewed             INTEGER,
    protein_existence    TEXT,
    annotation_score     INTEGER,
    canonical_uniprot_id TEXT,
    created_at           DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (tim_barrel_accession) REFERENCES entries(accession) ON DELETE CASCADE,
    FOREIGN KEY (canonical_uniprot_id) REFERENCES proteins(uniprot_id)
);

CREATE TABLE IF NOT EXISTS isoforms (
    isoform_id          TEXT PRIMARY KEY,
    uniprot_id          TEXT NOT NULL,
    is_canonical        INTEGER NOT NULL DEFAULT 0,
    sequence            TEXT NOT NULL,
    sequence_length     INTEGER NOT NULL,
    is_fragment         INTEGER NOT NULL DEFAULT 0,
    gene_name           TEXT,
    exon_count          INTEGER,
    exon_annotations    TEXT,
    splice_variants     TEXT,
    tim_barrel_location TEXT,
    tim_barrel_sequence TEXT,
    ensembl_transcript_id     TEXT,
    alphafold_id        TEXT,
    created_at          DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (uniprot_id) REFERENCES proteins(uniprot_id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS affected_isoforms (
    isoform_id                TEXT PRIMARY KEY,
    uniprot_id                TEXT NOT NULL,
    is_canonical              INTEGER NOT NULL DEFAULT 0,
    sequence                  TEXT NOT NULL,
    sequence_length           INTEGER NOT NULL,
    is_fragment               INTEGER NOT NULL DEFAULT 0,
    gene_name                 TEXT,
    exon_count                INTEGER,
    exon_annotations          TEXT,
    splice_variants           TEXT,
    domain_location           TEXT,
    domain_sequence           TEXT,
    canonical_domain_location TEXT,
    canonical_domain_sequence TEXT,
    identity_percentage               REAL NOT NULL,
    alignment_score                   INTEGER NOT NULL,
    exon_boundary_in_domain           INTEGER NOT NULL DEFAULT 0,
    exon_boundaries_in_domain_count   INTEGER NOT NULL DEFAULT 0,
    ensembl_transcript_id             TEXT,
    alphafold_id                      TEXT,
    vsp_domain_events                 TEXT,
    detection_method                  TEXT,
    disrupted_motifs                  TEXT,
    created_at                        DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (uniprot_id) REFERENCES proteins(uniprot_id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_proteins_canonical ON proteins(canonical_uniprot_id);
CREATE INDEX IF NOT EXISTS idx_proteins_entry     ON proteins(tim_barrel_accession);
CREATE INDEX IF NOT EXISTS idx_isoforms_uniprot   ON isoforms(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_isoforms_canonical ON isoforms(uniprot_id, is_canonical);
CREATE INDEX IF NOT EXISTS idx_isoforms_length    ON isoforms(sequence_length);
CREATE INDEX IF NOT EXISTS idx_affected_uniprot   ON affected_isoforms(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_affected_identity  ON affected_isoforms(identity_percentage);

CREATE TRIGGER IF NOT EXISTS trg_block_redundant
BEFORE INSERT ON isoforms
BEGIN
    SELECT CASE
        WHEN (SELECT canonical_uniprot_id FROM proteins WHERE uniprot_id = NEW.uniprot_id) IS NOT NULL
        THEN RAISE(ABORT, 'Isoform rejected: protein is redundant')
    END;
END;

-- ============================================================
-- Ensembl transcript expansion — TIM barrel, Homo sapiens
-- ============================================================

CREATE TABLE IF NOT EXISTS ensembl_transcripts (
    enst_id              TEXT PRIMARY KEY,
    ensg_id              TEXT,
    ensp_id              TEXT,
    uniprot_id           TEXT NOT NULL,
    gene_name            TEXT,
    sequence             TEXT,
    sequence_length      INTEGER,
    is_fragment          INTEGER NOT NULL DEFAULT 0,
    is_mane_select       INTEGER NOT NULL DEFAULT 0,
    biotype              TEXT,
    duplicate_isoform_id TEXT,
    duplicate_enst_id    TEXT,
    exon_annotations     TEXT,
    created_at           DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (uniprot_id) REFERENCES proteins(uniprot_id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS ensembl_affected (
    id                        INTEGER PRIMARY KEY AUTOINCREMENT,
    enst_id                   TEXT NOT NULL,
    uniprot_id                TEXT NOT NULL,
    gene_name                 TEXT,
    domain_location           TEXT,
    domain_sequence           TEXT,
    canonical_domain_location TEXT,
    canonical_domain_sequence TEXT,
    alignment_identity        REAL NOT NULL,
    alignment_score           INTEGER NOT NULL,
    insertion_detected                INTEGER NOT NULL DEFAULT 0,
    exon_boundary_in_domain           INTEGER NOT NULL DEFAULT 0,
    exon_boundaries_in_domain_count   INTEGER NOT NULL DEFAULT 0,
    disrupted_motifs                  TEXT,
    created_at                        DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (enst_id)    REFERENCES ensembl_transcripts(enst_id) ON DELETE CASCADE,
    FOREIGN KEY (uniprot_id) REFERENCES proteins(uniprot_id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_enst_uniprot   ON ensembl_transcripts(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_enst_gene      ON ensembl_transcripts(ensg_id);
CREATE INDEX IF NOT EXISTS idx_enst_dup       ON ensembl_transcripts(duplicate_isoform_id);
CREATE INDEX IF NOT EXISTS idx_enst_aff_enst  ON ensembl_affected(enst_id);
CREATE INDEX IF NOT EXISTS idx_enst_aff_uid   ON ensembl_affected(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_enst_aff_ident ON ensembl_affected(alignment_identity);

-- ============================================================
-- TIM barrel canonical analysis — Homo sapiens
-- One row per canonical protein, used for motif annotation.
-- exon_annotations: [{exon, start, end}] in full-sequence coords (1-based, inclusive)
-- motif_annotations: [{motif, start, end, beta_start, beta_end, alpha_start, alpha_end}] x8
-- ============================================================

CREATE TABLE IF NOT EXISTS canonical_analysis (
    uniprot_id            TEXT PRIMARY KEY,
    gene_name             TEXT,
    sequence              TEXT NOT NULL,
    domain_start          INTEGER,
    domain_end            INTEGER,
    domain_sequence       TEXT,
    exon_annotations      TEXT,
    motif_annotations     TEXT,
    dssp_source           TEXT,
    hmmer_source          TEXT,
    hmmer_annotations     TEXT,
    pdb_motif_annotations TEXT,
    pdb_source            TEXT,
    created_at            DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (uniprot_id) REFERENCES proteins(uniprot_id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_can_gene ON canonical_analysis(gene_name);
"""


def init_db(conn: sqlite3.Connection) -> None:
    """Create all tables, indexes, and triggers if they do not already exist."""
    conn.executescript(_CREATE_TABLES)
    conn.commit()
