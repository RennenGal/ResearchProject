-- Test Database Schema for Protein Data Collector
USE test_protein_data;

-- Unified table for both PFAM families and InterPro entries with TIM barrel annotations
CREATE TABLE IF NOT EXISTS tim_barrel_entries (
    accession VARCHAR(20) PRIMARY KEY,
    entry_type VARCHAR(20) NOT NULL,
    name VARCHAR(255) NOT NULL,
    description TEXT,
    interpro_type VARCHAR(50),
    tim_barrel_annotation TEXT NOT NULL,
    member_databases JSON,
    interpro_id VARCHAR(20),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    
    INDEX idx_tim_barrel_type (entry_type),
    INDEX idx_tim_barrel_name (name),
    INDEX idx_tim_barrel_interpro (interpro_id),
    INDEX idx_tim_barrel_created (created_at)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Human proteins belonging to TIM barrel entries (from InterPro)
-- Uses composite primary key to allow same protein in multiple entries
CREATE TABLE IF NOT EXISTS interpro_proteins (
    uniprot_id VARCHAR(20) NOT NULL,
    tim_barrel_accession VARCHAR(20) NOT NULL,
    name VARCHAR(255),
    organism VARCHAR(100) DEFAULT 'Homo sapiens',
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    
    PRIMARY KEY (uniprot_id, tim_barrel_accession),
    FOREIGN KEY (tim_barrel_accession) REFERENCES tim_barrel_entries(accession) ON DELETE CASCADE,
    INDEX idx_interpro_tim_barrel (tim_barrel_accession),
    INDEX idx_interpro_organism (organism),
    INDEX idx_interpro_name (name),
    INDEX idx_interpro_created (created_at)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

-- Protein isoforms with detailed annotations (from UniProt)
CREATE TABLE IF NOT EXISTS proteins (
    isoform_id VARCHAR(30) PRIMARY KEY,
    parent_protein_id VARCHAR(20) NOT NULL,
    parent_tim_barrel_accession VARCHAR(20) NOT NULL,
    sequence TEXT NOT NULL,
    sequence_length INTEGER NOT NULL,
    exon_annotations JSON,
    exon_count INTEGER,
    tim_barrel_location JSON,
    organism VARCHAR(100),
    name VARCHAR(255),
    description TEXT,
    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    
    FOREIGN KEY (parent_protein_id, parent_tim_barrel_accession) 
        REFERENCES interpro_proteins(uniprot_id, tim_barrel_accession) ON DELETE CASCADE,
    INDEX idx_proteins_parent (parent_protein_id),
    INDEX idx_proteins_parent_composite (parent_protein_id, parent_tim_barrel_accession),
    INDEX idx_proteins_organism (organism),
    INDEX idx_proteins_length (sequence_length),
    INDEX idx_proteins_exon_count (exon_count),
    INDEX idx_proteins_name (name),
    INDEX idx_proteins_created (created_at),
    
    -- Full-text search on sequence (for BLAST-like searches)
    FULLTEXT KEY ft_sequence (sequence),
    FULLTEXT KEY ft_description (description)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;

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
    p.parent_tim_barrel_accession as tim_barrel_accession,
    tbe.name as tim_barrel_name,
    tbe.entry_type as tim_barrel_type,
    CASE WHEN p.tim_barrel_location IS NOT NULL THEN 1 ELSE 0 END as has_tim_barrel,
    p.created_at
FROM proteins p
JOIN interpro_proteins ip ON p.parent_protein_id = ip.uniprot_id 
    AND p.parent_tim_barrel_accession = ip.tim_barrel_accession
JOIN tim_barrel_entries tbe ON ip.tim_barrel_accession = tbe.accession;

-- Create view for collection statistics
CREATE OR REPLACE VIEW collection_stats AS
SELECT 
    (SELECT COUNT(*) FROM tim_barrel_entries) as tim_barrel_entries_count,
    (SELECT COUNT(*) FROM tim_barrel_entries WHERE entry_type = 'pfam') as pfam_entries_count,
    (SELECT COUNT(*) FROM tim_barrel_entries WHERE entry_type = 'interpro') as interpro_entries_count,
    (SELECT COUNT(*) FROM interpro_proteins) as interpro_proteins_count,
    (SELECT COUNT(DISTINCT uniprot_id) FROM interpro_proteins) as unique_proteins_count,
    (SELECT COUNT(*) FROM proteins) as protein_isoforms_count,
    (SELECT COUNT(*) FROM proteins WHERE tim_barrel_location IS NOT NULL) as tim_barrel_proteins_count,
    (SELECT AVG(sequence_length) FROM proteins) as avg_sequence_length,
    (SELECT AVG(exon_count) FROM proteins WHERE exon_count IS NOT NULL) as avg_exon_count;

-- Show completion message
SELECT 'Test database schema created successfully' as status;