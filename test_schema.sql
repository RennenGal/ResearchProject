-- Test Database Schema for Protein Data Collector
USE test_protein_data;

-- PFAM families with TIM barrel annotations
CREATE TABLE IF NOT EXISTS pfam_families (
    accession VARCHAR(20) PRIMARY KEY,
    name VARCHAR(255) NOT NULL,
    description TEXT,
    tim_barrel_annotation TEXT NOT NULL,
    interpro_id VARCHAR(20),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    
    INDEX idx_pfam_name (name),
    INDEX idx_pfam_interpro (interpro_id),
    INDEX idx_pfam_created (created_at)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Human proteins belonging to PFAM families (from InterPro)
CREATE TABLE IF NOT EXISTS interpro_proteins (
    uniprot_id VARCHAR(20) PRIMARY KEY,
    pfam_accession VARCHAR(20) NOT NULL,
    name VARCHAR(255),
    organism VARCHAR(100) DEFAULT 'Homo sapiens',
    basic_metadata JSON,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    
    FOREIGN KEY (pfam_accession) REFERENCES pfam_families(accession) ON DELETE CASCADE,
    INDEX idx_interpro_pfam (pfam_accession),
    INDEX idx_interpro_organism (organism),
    INDEX idx_interpro_name (name),
    INDEX idx_interpro_created (created_at)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Protein isoforms with detailed annotations (from UniProt)
CREATE TABLE IF NOT EXISTS proteins (
    isoform_id VARCHAR(30) PRIMARY KEY,
    parent_protein_id VARCHAR(20) NOT NULL,
    sequence TEXT NOT NULL,
    sequence_length INTEGER NOT NULL,
    exon_annotations JSON,
    exon_count INTEGER,
    tim_barrel_location JSON,
    organism VARCHAR(100),
    name VARCHAR(255),
    description TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    
    FOREIGN KEY (parent_protein_id) REFERENCES interpro_proteins(uniprot_id) ON DELETE CASCADE,
    INDEX idx_proteins_parent (parent_protein_id),
    INDEX idx_proteins_organism (organism),
    INDEX idx_proteins_length (sequence_length),
    INDEX idx_proteins_exon_count (exon_count),
    INDEX idx_proteins_name (name),
    INDEX idx_proteins_created (created_at),
    
    -- Full-text search on sequence (for BLAST-like searches)
    FULLTEXT KEY ft_sequence (sequence),
    FULLTEXT KEY ft_description (description)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Collection progress tracking table (optional for research)
CREATE TABLE IF NOT EXISTS collection_progress (
    id INT AUTO_INCREMENT PRIMARY KEY,
    collection_id VARCHAR(50) NOT NULL UNIQUE,
    phase ENUM('pfam_families', 'interpro_proteins', 'protein_isoforms', 'completed', 'failed') NOT NULL,
    progress_data JSON,
    started_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    completed_at TIMESTAMP NULL,
    
    INDEX idx_progress_collection (collection_id),
    INDEX idx_progress_phase (phase),
    INDEX idx_progress_started (started_at)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- API logs table for tracking API usage
CREATE TABLE IF NOT EXISTS api_logs (
    id INT AUTO_INCREMENT PRIMARY KEY,
    api_name VARCHAR(50) NOT NULL,
    endpoint VARCHAR(255) NOT NULL,
    method VARCHAR(10) NOT NULL,
    status_code INT,
    response_time_ms INT,
    error_message TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    
    INDEX idx_api_logs_api (api_name),
    INDEX idx_api_logs_status (status_code),
    INDEX idx_api_logs_created (created_at)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Health metrics table for system monitoring
CREATE TABLE IF NOT EXISTS health_metrics (
    id INT AUTO_INCREMENT PRIMARY KEY,
    metric_name VARCHAR(100) NOT NULL,
    metric_value DECIMAL(10,2) NOT NULL,
    metric_unit VARCHAR(20),
    recorded_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    
    INDEX idx_health_metrics_name (metric_name),
    INDEX idx_health_metrics_recorded (recorded_at)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Create views for common research queries
CREATE OR REPLACE VIEW protein_summary AS
SELECT 
    p.isoform_id,
    p.parent_protein_id,
    p.name as protein_name,
    p.organism,
    p.sequence_length,
    p.exon_count,
    ip.pfam_accession,
    pf.name as pfam_name,
    CASE WHEN p.tim_barrel_location IS NOT NULL THEN 1 ELSE 0 END as has_tim_barrel,
    p.created_at
FROM proteins p
JOIN interpro_proteins ip ON p.parent_protein_id = ip.uniprot_id
JOIN pfam_families pf ON ip.pfam_accession = pf.accession;

-- Create view for collection statistics
CREATE OR REPLACE VIEW collection_stats AS
SELECT 
    (SELECT COUNT(*) FROM pfam_families) as pfam_families_count,
    (SELECT COUNT(*) FROM interpro_proteins) as interpro_proteins_count,
    (SELECT COUNT(*) FROM proteins) as protein_isoforms_count,
    (SELECT COUNT(*) FROM proteins WHERE tim_barrel_location IS NOT NULL) as tim_barrel_proteins_count,
    (SELECT AVG(sequence_length) FROM proteins) as avg_sequence_length,
    (SELECT AVG(exon_count) FROM proteins WHERE exon_count IS NOT NULL) as avg_exon_count;

-- Show completion message
SELECT 'Test database schema created successfully' as status;