-- SQLite Schema for Protein Data Collector
-- Generated from working database: protein_data.db
-- Contains focused schema with 67 UniProt fields across 9 categories
-- Optimized for TIM barrel protein research

CREATE TABLE tim_barrel_entries (
        accession TEXT PRIMARY KEY,
        entry_type TEXT NOT NULL CHECK (entry_type IN ('pfam', 'interpro')),
        name TEXT NOT NULL,
        description TEXT,
        interpro_type TEXT,
        tim_barrel_annotation TEXT NOT NULL,
        member_databases TEXT,
        interpro_id TEXT,
        created_at DATETIME DEFAULT CURRENT_TIMESTAMP
    );
CREATE TABLE interpro_proteins (
        uniprot_id TEXT NOT NULL,
        tim_barrel_accession TEXT NOT NULL,
        name TEXT,
        organism TEXT DEFAULT 'Homo sapiens',
        created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
        PRIMARY KEY (uniprot_id, tim_barrel_accession),
        FOREIGN KEY (tim_barrel_accession) REFERENCES tim_barrel_entries(accession) ON DELETE CASCADE
    );
CREATE TABLE proteins (
        -- Primary identifiers
        isoform_id TEXT PRIMARY KEY,
        parent_protein_id TEXT NOT NULL,
        parent_tim_barrel_accession TEXT NOT NULL,
        
        -- Category field for filtering
        data_category TEXT DEFAULT 'sequences' CHECK (data_category IN ('names_taxonomy', 'sequences', 'function', 'interaction', 'gene_ontology', 'structure', 'dates', 'family_domains', '3d_structure_databases')),
        
        -- 1. Names & Taxonomy (14 fields)
        accession TEXT,
        entry_name TEXT,
        gene_names TEXT,
        gene_primary TEXT,
        gene_synonym TEXT,
        gene_oln TEXT,
        gene_orf TEXT,
        organism_name TEXT,
        organism_id INTEGER,
        protein_name TEXT,
        proteomes TEXT,
        lineage TEXT,
        lineage_ids TEXT,
        virus_hosts TEXT,
        
        -- 2. Sequences (19 fields)
        alternative_products TEXT,
        alternative_sequence TEXT,
        error_gmodel_pred BOOLEAN,
        fragment BOOLEAN,
        organelle TEXT,
        sequence TEXT NOT NULL,
        sequence_length INTEGER NOT NULL,
        mass REAL,
        mass_spectrometry TEXT,
        natural_variant TEXT,
        non_adjacent_residues TEXT,
        non_standard_residue TEXT,
        non_terminal_residue TEXT,
        polymorphism TEXT,
        rna_editing TEXT,
        sequence_caution TEXT,
        sequence_conflict TEXT,
        sequence_uncertainty TEXT,
        sequence_version INTEGER,
        
        -- 3. Function (16 fields)
        absorption TEXT,
        active_site TEXT,
        activity_regulation TEXT,
        binding_site TEXT,
        catalytic_activity TEXT,
        cofactor TEXT,
        dna_binding TEXT,
        ec_number TEXT,
        function_cc TEXT,
        kinetics TEXT,
        pathway TEXT,
        ph_dependence TEXT,
        redox_potential TEXT,
        rhea_id TEXT,
        site TEXT,
        temp_dependence TEXT,
        
        -- 4. Interaction (2 fields)
        interacts_with TEXT,
        subunit_structure TEXT,
        
        -- 5. Gene Ontology (5 fields)
        go_biological_process TEXT,
        go_cellular_component TEXT,
        go_molecular_function TEXT,
        go_terms TEXT,
        go_ids TEXT,
        
        -- 6. Structure (4 fields)
        structure_3d TEXT,
        beta_strand TEXT,
        helix TEXT,
        turn TEXT,
        
        -- 7. Date Information (4 fields)
        date_created DATE,
        date_modified DATE,
        date_sequence_modified DATE,
        entry_version INTEGER,
        
        -- 8. Family & Domains (9 fields)
        coiled_coil TEXT,
        compositional_bias TEXT,
        domain_cc TEXT,
        domain_ft TEXT,
        motif TEXT,
        protein_families TEXT,
        region TEXT,
        repeat_region TEXT,
        zinc_finger TEXT,
        
        -- 9. 3D Structure Databases (7 fields)
        xref_alphafolddb TEXT,
        xref_bmrb TEXT,
        xref_pcddb TEXT,
        xref_pdb TEXT,
        xref_pdbsum TEXT,
        xref_sasbdb TEXT,
        xref_smr TEXT,
        
        -- Custom fields for TIM barrel research
        exon_annotations TEXT,
        exon_count INTEGER,
        tim_barrel_location TEXT,
        
        -- Quality indicators
        protein_existence TEXT,
        reviewed BOOLEAN,
        annotation_score INTEGER,
        
        -- Metadata
        created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
        updated_at DATETIME DEFAULT CURRENT_TIMESTAMP,
        
        -- Foreign key constraint
        FOREIGN KEY (parent_protein_id, parent_tim_barrel_accession) 
            REFERENCES interpro_proteins(uniprot_id, tim_barrel_accession) ON DELETE CASCADE
    );
CREATE TABLE collection_progress (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        collection_id TEXT NOT NULL UNIQUE,
        phase TEXT NOT NULL CHECK (phase IN ('pfam_families', 'interpro_proteins', 'protein_isoforms', 'completed', 'failed')),
        progress_data TEXT,
        started_at DATETIME DEFAULT CURRENT_TIMESTAMP,
        updated_at DATETIME DEFAULT CURRENT_TIMESTAMP,
        completed_at DATETIME
    );
CREATE TABLE sqlite_sequence(name,seq);
CREATE INDEX idx_tim_barrel_type ON tim_barrel_entries(entry_type);
CREATE INDEX idx_interpro_tim_barrel ON interpro_proteins(tim_barrel_accession);
CREATE INDEX idx_proteins_parent ON proteins(parent_protein_id);
CREATE INDEX idx_proteins_organism ON proteins(organism_name);
CREATE INDEX idx_proteins_length ON proteins(sequence_length);
CREATE INDEX idx_proteins_reviewed ON proteins(reviewed);
CREATE INDEX idx_proteins_category ON proteins(data_category);
CREATE VIEW protein_summary AS
    SELECT 
        p.isoform_id,
        p.parent_protein_id,
        p.protein_name,
        p.gene_primary,
        p.organism_name,
        p.sequence_length,
        p.exon_count,
        p.reviewed,
        p.protein_existence,
        p.annotation_score,
        p.data_category,
        p.parent_tim_barrel_accession as tim_barrel_accession,
        tbe.name as tim_barrel_name,
        tbe.entry_type as tim_barrel_type,
        CASE WHEN p.tim_barrel_location IS NOT NULL THEN 1 ELSE 0 END as has_tim_barrel,
        CASE WHEN p.xref_pdb IS NOT NULL AND p.xref_pdb != '' THEN 1 ELSE 0 END as has_3d_structure,
        p.created_at
    FROM proteins p
    JOIN interpro_proteins ip ON p.parent_protein_id = ip.uniprot_id 
        AND p.parent_tim_barrel_accession = ip.tim_barrel_accession
    JOIN tim_barrel_entries tbe ON ip.tim_barrel_accession = tbe.accession
/* protein_summary(isoform_id,parent_protein_id,protein_name,gene_primary,organism_name,sequence_length,exon_count,reviewed,protein_existence,annotation_score,data_category,tim_barrel_accession,tim_barrel_name,tim_barrel_type,has_tim_barrel,has_3d_structure,created_at) */;
CREATE VIEW collection_stats AS
    SELECT 
        (SELECT COUNT(*) FROM tim_barrel_entries) as tim_barrel_entries_count,
        (SELECT COUNT(*) FROM interpro_proteins) as interpro_proteins_count,
        (SELECT COUNT(*) FROM proteins) as protein_isoforms_count,
        (SELECT COUNT(*) FROM proteins WHERE reviewed = 1) as reviewed_proteins_count,
        (SELECT AVG(sequence_length) FROM proteins) as avg_sequence_length
/* collection_stats(tim_barrel_entries_count,interpro_proteins_count,protein_isoforms_count,reviewed_proteins_count,avg_sequence_length) */;
