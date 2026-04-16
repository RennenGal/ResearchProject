"""Database schema — SQL DDL for all domain/organism table sets.

Naming convention
-----------------
  {domain_prefix}_entries                          — domain family entries (one per domain)
  {domain_prefix}_proteins[_{organism}]            — collected proteins
  {domain_prefix}_isoforms[_{organism}]            — collected isoforms
  {domain_prefix}_affected_isoforms[_{organism}]   — AS-affected domain isoforms (analysis)

Domain prefixes : tb (TIM barrel), bp (beta propeller)
Organisms       : (none = Homo sapiens), _mus_musculus, _rattus_norvegicus
"""

import sqlite3


_CREATE_TABLES = """
-- ============================================================
-- Domain entries (one table per domain, shared across organisms)
-- ============================================================

CREATE TABLE IF NOT EXISTS tb_entries (
    accession         TEXT PRIMARY KEY,
    entry_type        TEXT NOT NULL CHECK (entry_type IN ('pfam', 'interpro')),
    name              TEXT NOT NULL,
    description       TEXT,
    domain_annotation TEXT NOT NULL,
    created_at        DATETIME DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS bp_entries (
    accession         TEXT PRIMARY KEY,
    entry_type        TEXT NOT NULL CHECK (entry_type IN ('pfam', 'interpro')),
    name              TEXT NOT NULL,
    description       TEXT,
    domain_annotation TEXT NOT NULL,
    created_at        DATETIME DEFAULT CURRENT_TIMESTAMP
);

-- ============================================================
-- TIM barrel — Homo sapiens
-- ============================================================

CREATE TABLE IF NOT EXISTS tb_proteins (
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
    FOREIGN KEY (tim_barrel_accession) REFERENCES tb_entries(accession) ON DELETE CASCADE,
    FOREIGN KEY (canonical_uniprot_id) REFERENCES tb_proteins(uniprot_id)
);

CREATE TABLE IF NOT EXISTS tb_isoforms (
    isoform_id          TEXT PRIMARY KEY,
    uniprot_id          TEXT NOT NULL,
    is_canonical        INTEGER NOT NULL DEFAULT 0,
    sequence            TEXT NOT NULL,
    sequence_length     INTEGER NOT NULL,
    is_fragment         INTEGER NOT NULL DEFAULT 0,
    exon_count          INTEGER,
    exon_annotations    TEXT,
    splice_variants     TEXT,
    tim_barrel_location TEXT,
    tim_barrel_sequence TEXT,
    ensembl_gene_id     TEXT,
    alphafold_id        TEXT,
    created_at          DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (uniprot_id) REFERENCES tb_proteins(uniprot_id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS tb_affected_isoforms (
    isoform_id                TEXT PRIMARY KEY,
    uniprot_id                TEXT NOT NULL,
    is_canonical              INTEGER NOT NULL DEFAULT 0,
    sequence                  TEXT NOT NULL,
    sequence_length           INTEGER NOT NULL,
    is_fragment               INTEGER NOT NULL DEFAULT 0,
    exon_count                INTEGER,
    exon_annotations          TEXT,
    splice_variants           TEXT,
    domain_location           TEXT,
    domain_sequence           TEXT,
    canonical_domain_location TEXT,
    canonical_domain_sequence TEXT,
    identity_percentage       REAL NOT NULL,
    alignment_score           INTEGER NOT NULL,
    ensembl_gene_id           TEXT,
    alphafold_id              TEXT,
    created_at                DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (uniprot_id) REFERENCES tb_proteins(uniprot_id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_tb_proteins_canonical ON tb_proteins(canonical_uniprot_id);
CREATE INDEX IF NOT EXISTS idx_tb_proteins_entry     ON tb_proteins(tim_barrel_accession);
CREATE INDEX IF NOT EXISTS idx_tb_isoforms_uniprot   ON tb_isoforms(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_tb_isoforms_canonical ON tb_isoforms(uniprot_id, is_canonical);
CREATE INDEX IF NOT EXISTS idx_tb_isoforms_length    ON tb_isoforms(sequence_length);
CREATE INDEX IF NOT EXISTS idx_tb_affected_uniprot   ON tb_affected_isoforms(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_tb_affected_identity  ON tb_affected_isoforms(identity_percentage);

CREATE TRIGGER IF NOT EXISTS trg_block_redundant_tb
BEFORE INSERT ON tb_isoforms
BEGIN
    SELECT CASE
        WHEN (SELECT canonical_uniprot_id FROM tb_proteins WHERE uniprot_id = NEW.uniprot_id) IS NOT NULL
        THEN RAISE(ABORT, 'Isoform rejected: protein is redundant')
    END;
END;

-- ============================================================
-- TIM barrel — Mus musculus
-- ============================================================

CREATE TABLE IF NOT EXISTS tb_proteins_mus_musculus (
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
    FOREIGN KEY (tim_barrel_accession) REFERENCES tb_entries(accession) ON DELETE CASCADE,
    FOREIGN KEY (canonical_uniprot_id) REFERENCES tb_proteins_mus_musculus(uniprot_id)
);

CREATE TABLE IF NOT EXISTS tb_isoforms_mus_musculus (
    isoform_id          TEXT PRIMARY KEY,
    uniprot_id          TEXT NOT NULL,
    is_canonical        INTEGER NOT NULL DEFAULT 0,
    sequence            TEXT NOT NULL,
    sequence_length     INTEGER NOT NULL,
    is_fragment         INTEGER NOT NULL DEFAULT 0,
    exon_count          INTEGER,
    exon_annotations    TEXT,
    splice_variants     TEXT,
    tim_barrel_location TEXT,
    tim_barrel_sequence TEXT,
    ensembl_gene_id     TEXT,
    alphafold_id        TEXT,
    created_at          DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (uniprot_id) REFERENCES tb_proteins_mus_musculus(uniprot_id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS tb_affected_isoforms_mus_musculus (
    isoform_id                TEXT PRIMARY KEY,
    uniprot_id                TEXT NOT NULL,
    is_canonical              INTEGER NOT NULL DEFAULT 0,
    sequence                  TEXT NOT NULL,
    sequence_length           INTEGER NOT NULL,
    is_fragment               INTEGER NOT NULL DEFAULT 0,
    exon_count                INTEGER,
    exon_annotations          TEXT,
    splice_variants           TEXT,
    domain_location           TEXT,
    domain_sequence           TEXT,
    canonical_domain_location TEXT,
    canonical_domain_sequence TEXT,
    identity_percentage       REAL NOT NULL,
    alignment_score           INTEGER NOT NULL,
    ensembl_gene_id           TEXT,
    alphafold_id              TEXT,
    created_at                DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (uniprot_id) REFERENCES tb_proteins_mus_musculus(uniprot_id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_tb_mm_proteins_canonical ON tb_proteins_mus_musculus(canonical_uniprot_id);
CREATE INDEX IF NOT EXISTS idx_tb_mm_proteins_entry     ON tb_proteins_mus_musculus(tim_barrel_accession);
CREATE INDEX IF NOT EXISTS idx_tb_mm_isoforms_uniprot   ON tb_isoforms_mus_musculus(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_tb_mm_isoforms_canonical ON tb_isoforms_mus_musculus(uniprot_id, is_canonical);
CREATE INDEX IF NOT EXISTS idx_tb_mm_isoforms_length    ON tb_isoforms_mus_musculus(sequence_length);
CREATE INDEX IF NOT EXISTS idx_tb_mm_affected_uniprot   ON tb_affected_isoforms_mus_musculus(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_tb_mm_affected_identity  ON tb_affected_isoforms_mus_musculus(identity_percentage);

CREATE TRIGGER IF NOT EXISTS trg_block_redundant_tb_mm
BEFORE INSERT ON tb_isoforms_mus_musculus
BEGIN
    SELECT CASE
        WHEN (SELECT canonical_uniprot_id FROM tb_proteins_mus_musculus WHERE uniprot_id = NEW.uniprot_id) IS NOT NULL
        THEN RAISE(ABORT, 'Isoform rejected: protein is redundant')
    END;
END;

-- ============================================================
-- TIM barrel — Rattus norvegicus
-- ============================================================

CREATE TABLE IF NOT EXISTS tb_proteins_rattus_norvegicus (
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
    FOREIGN KEY (tim_barrel_accession) REFERENCES tb_entries(accession) ON DELETE CASCADE,
    FOREIGN KEY (canonical_uniprot_id) REFERENCES tb_proteins_rattus_norvegicus(uniprot_id)
);

CREATE TABLE IF NOT EXISTS tb_isoforms_rattus_norvegicus (
    isoform_id          TEXT PRIMARY KEY,
    uniprot_id          TEXT NOT NULL,
    is_canonical        INTEGER NOT NULL DEFAULT 0,
    sequence            TEXT NOT NULL,
    sequence_length     INTEGER NOT NULL,
    is_fragment         INTEGER NOT NULL DEFAULT 0,
    exon_count          INTEGER,
    exon_annotations    TEXT,
    splice_variants     TEXT,
    tim_barrel_location TEXT,
    tim_barrel_sequence TEXT,
    ensembl_gene_id     TEXT,
    alphafold_id        TEXT,
    created_at          DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (uniprot_id) REFERENCES tb_proteins_rattus_norvegicus(uniprot_id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS tb_affected_isoforms_rattus_norvegicus (
    isoform_id                TEXT PRIMARY KEY,
    uniprot_id                TEXT NOT NULL,
    is_canonical              INTEGER NOT NULL DEFAULT 0,
    sequence                  TEXT NOT NULL,
    sequence_length           INTEGER NOT NULL,
    is_fragment               INTEGER NOT NULL DEFAULT 0,
    exon_count                INTEGER,
    exon_annotations          TEXT,
    splice_variants           TEXT,
    domain_location           TEXT,
    domain_sequence           TEXT,
    canonical_domain_location TEXT,
    canonical_domain_sequence TEXT,
    identity_percentage       REAL NOT NULL,
    alignment_score           INTEGER NOT NULL,
    ensembl_gene_id           TEXT,
    alphafold_id              TEXT,
    created_at                DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (uniprot_id) REFERENCES tb_proteins_rattus_norvegicus(uniprot_id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_tb_rn_proteins_canonical ON tb_proteins_rattus_norvegicus(canonical_uniprot_id);
CREATE INDEX IF NOT EXISTS idx_tb_rn_proteins_entry     ON tb_proteins_rattus_norvegicus(tim_barrel_accession);
CREATE INDEX IF NOT EXISTS idx_tb_rn_isoforms_uniprot   ON tb_isoforms_rattus_norvegicus(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_tb_rn_isoforms_canonical ON tb_isoforms_rattus_norvegicus(uniprot_id, is_canonical);
CREATE INDEX IF NOT EXISTS idx_tb_rn_isoforms_length    ON tb_isoforms_rattus_norvegicus(sequence_length);
CREATE INDEX IF NOT EXISTS idx_tb_rn_affected_uniprot   ON tb_affected_isoforms_rattus_norvegicus(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_tb_rn_affected_identity  ON tb_affected_isoforms_rattus_norvegicus(identity_percentage);

CREATE TRIGGER IF NOT EXISTS trg_block_redundant_tb_rn
BEFORE INSERT ON tb_isoforms_rattus_norvegicus
BEGIN
    SELECT CASE
        WHEN (SELECT canonical_uniprot_id FROM tb_proteins_rattus_norvegicus WHERE uniprot_id = NEW.uniprot_id) IS NOT NULL
        THEN RAISE(ABORT, 'Isoform rejected: protein is redundant')
    END;
END;

-- ============================================================
-- Beta propeller — Homo sapiens
-- ============================================================

CREATE TABLE IF NOT EXISTS bp_proteins (
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
    FOREIGN KEY (tim_barrel_accession) REFERENCES bp_entries(accession) ON DELETE CASCADE,
    FOREIGN KEY (canonical_uniprot_id) REFERENCES bp_proteins(uniprot_id)
);

CREATE TABLE IF NOT EXISTS bp_isoforms (
    isoform_id          TEXT PRIMARY KEY,
    uniprot_id          TEXT NOT NULL,
    is_canonical        INTEGER NOT NULL DEFAULT 0,
    sequence            TEXT NOT NULL,
    sequence_length     INTEGER NOT NULL,
    is_fragment         INTEGER NOT NULL DEFAULT 0,
    exon_count          INTEGER,
    exon_annotations    TEXT,
    splice_variants     TEXT,
    tim_barrel_location TEXT,
    tim_barrel_sequence TEXT,
    ensembl_gene_id     TEXT,
    alphafold_id        TEXT,
    created_at          DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (uniprot_id) REFERENCES bp_proteins(uniprot_id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS bp_affected_isoforms (
    isoform_id                TEXT PRIMARY KEY,
    uniprot_id                TEXT NOT NULL,
    is_canonical              INTEGER NOT NULL DEFAULT 0,
    sequence                  TEXT NOT NULL,
    sequence_length           INTEGER NOT NULL,
    is_fragment               INTEGER NOT NULL DEFAULT 0,
    exon_count                INTEGER,
    exon_annotations          TEXT,
    splice_variants           TEXT,
    domain_location           TEXT,
    domain_sequence           TEXT,
    canonical_domain_location TEXT,
    canonical_domain_sequence TEXT,
    identity_percentage       REAL NOT NULL,
    alignment_score           INTEGER NOT NULL,
    ensembl_gene_id           TEXT,
    alphafold_id              TEXT,
    created_at                DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (uniprot_id) REFERENCES bp_proteins(uniprot_id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_bp_proteins_canonical ON bp_proteins(canonical_uniprot_id);
CREATE INDEX IF NOT EXISTS idx_bp_proteins_entry     ON bp_proteins(tim_barrel_accession);
CREATE INDEX IF NOT EXISTS idx_bp_isoforms_uniprot   ON bp_isoforms(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_bp_isoforms_canonical ON bp_isoforms(uniprot_id, is_canonical);
CREATE INDEX IF NOT EXISTS idx_bp_isoforms_length    ON bp_isoforms(sequence_length);
CREATE INDEX IF NOT EXISTS idx_bp_affected_uniprot   ON bp_affected_isoforms(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_bp_affected_identity  ON bp_affected_isoforms(identity_percentage);

CREATE TRIGGER IF NOT EXISTS trg_block_redundant_bp
BEFORE INSERT ON bp_isoforms
BEGIN
    SELECT CASE
        WHEN (SELECT canonical_uniprot_id FROM bp_proteins WHERE uniprot_id = NEW.uniprot_id) IS NOT NULL
        THEN RAISE(ABORT, 'Isoform rejected: protein is redundant')
    END;
END;

-- ============================================================
-- Beta propeller — Mus musculus
-- ============================================================

CREATE TABLE IF NOT EXISTS bp_proteins_mus_musculus (
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
    FOREIGN KEY (tim_barrel_accession) REFERENCES bp_entries(accession) ON DELETE CASCADE,
    FOREIGN KEY (canonical_uniprot_id) REFERENCES bp_proteins_mus_musculus(uniprot_id)
);

CREATE TABLE IF NOT EXISTS bp_isoforms_mus_musculus (
    isoform_id          TEXT PRIMARY KEY,
    uniprot_id          TEXT NOT NULL,
    is_canonical        INTEGER NOT NULL DEFAULT 0,
    sequence            TEXT NOT NULL,
    sequence_length     INTEGER NOT NULL,
    is_fragment         INTEGER NOT NULL DEFAULT 0,
    exon_count          INTEGER,
    exon_annotations    TEXT,
    splice_variants     TEXT,
    tim_barrel_location TEXT,
    tim_barrel_sequence TEXT,
    ensembl_gene_id     TEXT,
    alphafold_id        TEXT,
    created_at          DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (uniprot_id) REFERENCES bp_proteins_mus_musculus(uniprot_id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS bp_affected_isoforms_mus_musculus (
    isoform_id                TEXT PRIMARY KEY,
    uniprot_id                TEXT NOT NULL,
    is_canonical              INTEGER NOT NULL DEFAULT 0,
    sequence                  TEXT NOT NULL,
    sequence_length           INTEGER NOT NULL,
    is_fragment               INTEGER NOT NULL DEFAULT 0,
    exon_count                INTEGER,
    exon_annotations          TEXT,
    splice_variants           TEXT,
    domain_location           TEXT,
    domain_sequence           TEXT,
    canonical_domain_location TEXT,
    canonical_domain_sequence TEXT,
    identity_percentage       REAL NOT NULL,
    alignment_score           INTEGER NOT NULL,
    ensembl_gene_id           TEXT,
    alphafold_id              TEXT,
    created_at                DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (uniprot_id) REFERENCES bp_proteins_mus_musculus(uniprot_id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_bp_mm_proteins_canonical ON bp_proteins_mus_musculus(canonical_uniprot_id);
CREATE INDEX IF NOT EXISTS idx_bp_mm_proteins_entry     ON bp_proteins_mus_musculus(tim_barrel_accession);
CREATE INDEX IF NOT EXISTS idx_bp_mm_isoforms_uniprot   ON bp_isoforms_mus_musculus(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_bp_mm_isoforms_canonical ON bp_isoforms_mus_musculus(uniprot_id, is_canonical);
CREATE INDEX IF NOT EXISTS idx_bp_mm_isoforms_length    ON bp_isoforms_mus_musculus(sequence_length);
CREATE INDEX IF NOT EXISTS idx_bp_mm_affected_uniprot   ON bp_affected_isoforms_mus_musculus(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_bp_mm_affected_identity  ON bp_affected_isoforms_mus_musculus(identity_percentage);

CREATE TRIGGER IF NOT EXISTS trg_block_redundant_bp_mm
BEFORE INSERT ON bp_isoforms_mus_musculus
BEGIN
    SELECT CASE
        WHEN (SELECT canonical_uniprot_id FROM bp_proteins_mus_musculus WHERE uniprot_id = NEW.uniprot_id) IS NOT NULL
        THEN RAISE(ABORT, 'Isoform rejected: protein is redundant')
    END;
END;

-- ============================================================
-- Beta propeller — Rattus norvegicus
-- ============================================================

CREATE TABLE IF NOT EXISTS bp_proteins_rattus_norvegicus (
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
    FOREIGN KEY (tim_barrel_accession) REFERENCES bp_entries(accession) ON DELETE CASCADE,
    FOREIGN KEY (canonical_uniprot_id) REFERENCES bp_proteins_rattus_norvegicus(uniprot_id)
);

CREATE TABLE IF NOT EXISTS bp_isoforms_rattus_norvegicus (
    isoform_id          TEXT PRIMARY KEY,
    uniprot_id          TEXT NOT NULL,
    is_canonical        INTEGER NOT NULL DEFAULT 0,
    sequence            TEXT NOT NULL,
    sequence_length     INTEGER NOT NULL,
    is_fragment         INTEGER NOT NULL DEFAULT 0,
    exon_count          INTEGER,
    exon_annotations    TEXT,
    splice_variants     TEXT,
    tim_barrel_location TEXT,
    tim_barrel_sequence TEXT,
    ensembl_gene_id     TEXT,
    alphafold_id        TEXT,
    created_at          DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (uniprot_id) REFERENCES bp_proteins_rattus_norvegicus(uniprot_id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS bp_affected_isoforms_rattus_norvegicus (
    isoform_id                TEXT PRIMARY KEY,
    uniprot_id                TEXT NOT NULL,
    is_canonical              INTEGER NOT NULL DEFAULT 0,
    sequence                  TEXT NOT NULL,
    sequence_length           INTEGER NOT NULL,
    is_fragment               INTEGER NOT NULL DEFAULT 0,
    exon_count                INTEGER,
    exon_annotations          TEXT,
    splice_variants           TEXT,
    domain_location           TEXT,
    domain_sequence           TEXT,
    canonical_domain_location TEXT,
    canonical_domain_sequence TEXT,
    identity_percentage       REAL NOT NULL,
    alignment_score           INTEGER NOT NULL,
    ensembl_gene_id           TEXT,
    alphafold_id              TEXT,
    created_at                DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (uniprot_id) REFERENCES bp_proteins_rattus_norvegicus(uniprot_id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_bp_rn_proteins_canonical ON bp_proteins_rattus_norvegicus(canonical_uniprot_id);
CREATE INDEX IF NOT EXISTS idx_bp_rn_proteins_entry     ON bp_proteins_rattus_norvegicus(tim_barrel_accession);
CREATE INDEX IF NOT EXISTS idx_bp_rn_isoforms_uniprot   ON bp_isoforms_rattus_norvegicus(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_bp_rn_isoforms_canonical ON bp_isoforms_rattus_norvegicus(uniprot_id, is_canonical);
CREATE INDEX IF NOT EXISTS idx_bp_rn_isoforms_length    ON bp_isoforms_rattus_norvegicus(sequence_length);
CREATE INDEX IF NOT EXISTS idx_bp_rn_affected_uniprot   ON bp_affected_isoforms_rattus_norvegicus(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_bp_rn_affected_identity  ON bp_affected_isoforms_rattus_norvegicus(identity_percentage);

CREATE TRIGGER IF NOT EXISTS trg_block_redundant_bp_rn
BEFORE INSERT ON bp_isoforms_rattus_norvegicus
BEGIN
    SELECT CASE
        WHEN (SELECT canonical_uniprot_id FROM bp_proteins_rattus_norvegicus WHERE uniprot_id = NEW.uniprot_id) IS NOT NULL
        THEN RAISE(ABORT, 'Isoform rejected: protein is redundant')
    END;
END;
"""


def init_db(conn: sqlite3.Connection) -> None:
    """Create all tables, indexes, and triggers if they do not already exist."""
    conn.executescript(_CREATE_TABLES)
    conn.commit()
