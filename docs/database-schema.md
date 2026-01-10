# Database Schema Documentation

## Overview

The Protein Data Collector uses a unified database schema designed to efficiently store and query TIM barrel protein data collected from InterPro and UniProt databases. The schema follows a three-tier hierarchical structure that reflects the biological relationships between protein families, proteins, and their isoforms.

## Schema Architecture

The database schema has been designed with the following principles:

- **Unified Structure**: Single table for both PFAM families and InterPro entries
- **Composite Primary Keys**: Support for proteins belonging to multiple TIM barrel entries
- **Referential Integrity**: Foreign key constraints ensure data consistency
- **Performance Optimization**: Strategic indexes for common query patterns
- **Extensibility**: JSON fields for flexible metadata storage

## Database Tables

### 1. tim_barrel_entries

**Purpose**: Unified storage for both PFAM families and InterPro entries with TIM barrel annotations.

**Structure**:
```sql
CREATE TABLE tim_barrel_entries (
    accession VARCHAR(20) PRIMARY KEY,           -- Entry accession (PF##### or IPR######)
    entry_type VARCHAR(20) NOT NULL,             -- 'pfam' or 'interpro'
    name VARCHAR(255) NOT NULL,                  -- Entry name
    description TEXT,                            -- Detailed description
    interpro_type VARCHAR(50),                   -- InterPro type (Domain, Family, etc.)
    tim_barrel_annotation TEXT NOT NULL,         -- TIM barrel annotation details
    member_databases JSON,                       -- Member databases (InterPro entries only)
    interpro_id VARCHAR(20),                     -- Associated InterPro ID (PFAM entries only)
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    
    -- Indexes for efficient querying
    INDEX idx_tim_barrel_type (entry_type),
    INDEX idx_tim_barrel_name (name),
    INDEX idx_tim_barrel_interpro (interpro_id),
    INDEX idx_tim_barrel_created (created_at)
);
```

**Key Features**:
- **Unified Design**: Stores both PFAM families (PF#####) and InterPro entries (IPR######)
- **Entry Type Discrimination**: `entry_type` field distinguishes between 'pfam' and 'interpro'
- **Flexible Metadata**: JSON fields for complex data structures
- **Cross-References**: Links between PFAM families and their InterPro entries

**Sample Data**:
```sql
-- PFAM family example
INSERT INTO tim_barrel_entries VALUES (
    'PF00113', 'pfam', 'Enolase, C-terminal TIM barrel domain',
    'TIM barrel domain found in enolases', NULL,
    'Eight-fold alpha/beta barrel structure characteristic of TIM barrel proteins',
    '{}', 'IPR000322', NOW()
);

-- InterPro entry example  
INSERT INTO tim_barrel_entries VALUES (
    'IPR013785', 'interpro', 'Aldolase-type TIM barrel',
    'TIM barrel found in aldolase-type enzymes', 'Domain',
    'Aldolase-type TIM barrel domain with specific structural features',
    '{"pfam": ["PF00274"], "smart": ["SM00849"]}', NULL, NOW()
);
```

### 2. interpro_proteins

**Purpose**: Human proteins belonging to TIM barrel entries, collected from InterPro API.

**Structure**:
```sql
CREATE TABLE interpro_proteins (
    uniprot_id VARCHAR(20) NOT NULL,             -- UniProt protein identifier
    tim_barrel_accession VARCHAR(20) NOT NULL,   -- TIM barrel entry accession
    name VARCHAR(255),                           -- Protein name
    organism VARCHAR(100) DEFAULT 'Homo sapiens', -- Source organism
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    
    -- Composite primary key allows same protein in multiple entries
    PRIMARY KEY (uniprot_id, tim_barrel_accession),
    
    -- Foreign key to TIM barrel entries
    FOREIGN KEY (tim_barrel_accession) 
        REFERENCES tim_barrel_entries(accession) ON DELETE CASCADE,
    
    -- Indexes for efficient querying
    INDEX idx_interpro_tim_barrel (tim_barrel_accession),
    INDEX idx_interpro_organism (organism),
    INDEX idx_interpro_name (name),
    INDEX idx_interpro_created (created_at)
);
```

**Key Features**:
- **Composite Primary Key**: `(uniprot_id, tim_barrel_accession)` allows proteins to belong to multiple TIM barrel entries
- **Human Focus**: Primarily stores Homo sapiens proteins
- **Referential Integrity**: Foreign key ensures valid TIM barrel entry references
- **Cascade Deletion**: Proteins are removed when their TIM barrel entry is deleted

**Sample Data**:
```sql
-- Same protein in multiple TIM barrel entries
INSERT INTO interpro_proteins VALUES 
    ('P06733', 'PF00113', 'Alpha-enolase', 'Homo sapiens', NOW()),
    ('P06733', 'IPR000322', 'Alpha-enolase', 'Homo sapiens', NOW());
```

### 3. proteins

**Purpose**: Detailed protein isoform data collected from UniProt, including sequences and structural annotations.

**Structure**:
```sql
CREATE TABLE proteins (
    isoform_id VARCHAR(30) PRIMARY KEY,          -- UniProt isoform identifier
    parent_protein_id VARCHAR(20) NOT NULL,      -- Parent protein UniProt ID
    parent_tim_barrel_accession VARCHAR(20) NOT NULL, -- Parent TIM barrel accession
    sequence TEXT NOT NULL,                      -- Amino acid sequence
    sequence_length INTEGER NOT NULL,            -- Sequence length
    exon_annotations JSON,                       -- Exon structure data
    exon_count INTEGER,                          -- Number of exons
    tim_barrel_location JSON,                    -- TIM barrel coordinates
    organism VARCHAR(100),                       -- Source organism
    name VARCHAR(255),                           -- Protein name
    description TEXT,                            -- Protein description
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    
    -- Composite foreign key to interpro_proteins
    FOREIGN KEY (parent_protein_id, parent_tim_barrel_accession) 
        REFERENCES interpro_proteins(uniprot_id, tim_barrel_accession) ON DELETE CASCADE,
    
    -- Indexes for efficient querying
    INDEX idx_proteins_parent (parent_protein_id),
    INDEX idx_proteins_parent_composite (parent_protein_id, parent_tim_barrel_accession),
    INDEX idx_proteins_organism (organism),
    INDEX idx_proteins_length (sequence_length),
    INDEX idx_proteins_exon_count (exon_count),
    INDEX idx_proteins_name (name),
    INDEX idx_proteins_created (created_at),
    
    -- Full-text search capabilities
    FULLTEXT KEY ft_sequence (sequence),
    FULLTEXT KEY ft_description (description)
);
```

**Key Features**:
- **Detailed Annotations**: Complete protein sequences and structural data
- **JSON Flexibility**: Complex annotations stored as JSON for flexibility
- **Composite Foreign Key**: Links to specific protein-TIM barrel combinations
- **Full-Text Search**: Enables sequence and description searching
- **Performance Indexes**: Optimized for common query patterns

**Sample Data**:
```sql
INSERT INTO proteins VALUES (
    'P06733-1', 'P06733', 'PF00113',
    'MSILKIHAREIFDSRGNPTVEVDLFTSKGLFRAAVPSGASTGIYEALELRDNDKTRYMGKGVSKAVEHINKTIAPALVSKKLNVTEQEKIDKLMIEMDGTENKSKFGANAILGVSLAVCKAGAVEKGVPLYRHIADLAGNSEVILPVPAFNVINGGSHAGNKLAMQEFMILPVGAANFREAMRIGAEVYHNLKNVIKEKYGKDATNVGDEGGFAPNILENKEGLELLKTAIGKAGYTDKVVIGMDVAASEFFRSGKYDLDFKSPDDPSRYISPDQLADLYKSFIKDYPVVSIEDPFDQDDWGAWQKFTASAGIQVVGDDLTVTNPKRIAKAVNEKSCNCLLLKVNQIGSVTESLQACKLAQANGWGVMVSHRSGETEDTFIADLVVGLCTGQIKTGAPCRSERLAKYNQLLRIEEELGSKAKFAGRNFRRNPLAK',
    434,
    '{"exons": [{"start": 1, "end": 150}, {"start": 151, "end": 300}, {"start": 301, "end": 434}]}',
    3,
    '{"start": 50, "end": 400, "confidence": 0.95, "method": "structural_alignment"}',
    'Homo sapiens', 'Alpha-enolase',
    'Multifunctional enzyme that catalyzes the interconversion of 2-phosphoglycerate and phosphoenolpyruvate',
    NOW()
);
```

## Database Views

### 1. protein_summary

**Purpose**: Simplified view combining protein data with TIM barrel entry information.

```sql
CREATE VIEW protein_summary AS
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
```

### 2. collection_stats

**Purpose**: Real-time statistics about the collected data.

```sql
CREATE VIEW collection_stats AS
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
```

## Schema Evolution History

### Migration from Old Schema

The current unified schema represents a significant improvement over the original design:

**Old Schema Issues**:
- Separate `pfam_families` and `interpro_entries` tables (redundant)
- Single primary key in `interpro_proteins` table (prevented proteins from belonging to multiple entries)
- No support for composite relationships

**New Schema Benefits**:
- Unified `tim_barrel_entries` table for both PFAM and InterPro data
- Composite primary keys support many-to-many relationships
- Improved referential integrity with proper foreign key constraints
- Better performance through strategic indexing

**Migration Process**:
1. Created new unified `tim_barrel_entries` table
2. Migrated data from separate `pfam_families` and `interpro_entries` tables
3. Recreated `interpro_proteins` table with composite primary key
4. Updated `proteins` table with composite foreign key constraints
5. Removed redundant old tables

## Data Relationships

### Entity Relationship Diagram

```
tim_barrel_entries (1) ──────── (M) interpro_proteins (1) ──────── (M) proteins
     │                                    │                              │
     │ accession                          │ (uniprot_id,                 │ (parent_protein_id,
     │                                    │  tim_barrel_accession)       │  parent_tim_barrel_accession)
     └────────────────────────────────────┘                              │
                                                                         │ isoform_id
```

### Relationship Details

1. **TIM Barrel Entries → InterPro Proteins**: One-to-Many
   - Each TIM barrel entry can have multiple associated proteins
   - Each protein-entry combination is unique (composite primary key)

2. **InterPro Proteins → Proteins**: One-to-Many
   - Each InterPro protein can have multiple isoforms
   - Isoforms are linked via composite foreign key

3. **Cross-Entry Relationships**: Many-to-Many (via composite keys)
   - Same protein can belong to multiple TIM barrel entries
   - Same protein can have isoforms associated with different entries

## Query Patterns

### Common Queries

**1. Get all proteins for a specific TIM barrel entry:**
```sql
SELECT ip.uniprot_id, ip.name, COUNT(p.isoform_id) as isoform_count
FROM interpro_proteins ip
LEFT JOIN proteins p ON ip.uniprot_id = p.parent_protein_id 
    AND ip.tim_barrel_accession = p.parent_tim_barrel_accession
WHERE ip.tim_barrel_accession = 'PF00113'
GROUP BY ip.uniprot_id, ip.name;
```

**2. Find proteins with TIM barrel structural annotations:**
```sql
SELECT p.isoform_id, p.name, p.sequence_length, 
       JSON_EXTRACT(p.tim_barrel_location, '$.confidence') as confidence
FROM proteins p
WHERE p.tim_barrel_location IS NOT NULL
ORDER BY JSON_EXTRACT(p.tim_barrel_location, '$.confidence') DESC;
```

**3. Get collection statistics by entry type:**
```sql
SELECT tbe.entry_type, 
       COUNT(DISTINCT ip.uniprot_id) as unique_proteins,
       COUNT(p.isoform_id) as total_isoforms,
       AVG(p.sequence_length) as avg_length
FROM tim_barrel_entries tbe
LEFT JOIN interpro_proteins ip ON tbe.accession = ip.tim_barrel_accession
LEFT JOIN proteins p ON ip.uniprot_id = p.parent_protein_id 
    AND ip.tim_barrel_accession = p.parent_tim_barrel_accession
GROUP BY tbe.entry_type;
```

**4. Search proteins by sequence similarity (full-text):**
```sql
SELECT p.isoform_id, p.name, p.sequence_length,
       MATCH(p.sequence) AGAINST('GLYCERALDEHYDE' IN NATURAL LANGUAGE MODE) as relevance
FROM proteins p
WHERE MATCH(p.sequence) AGAINST('GLYCERALDEHYDE' IN NATURAL LANGUAGE MODE)
ORDER BY relevance DESC;
```

## Performance Considerations

### Indexing Strategy

**Primary Indexes**:
- Primary keys provide clustered indexes for fast lookups
- Composite primary key in `interpro_proteins` enables efficient many-to-many queries

**Secondary Indexes**:
- Foreign key indexes for join performance
- Entry type index for filtering by PFAM vs InterPro
- Sequence length and exon count for range queries
- Creation timestamp for temporal queries

**Full-Text Indexes**:
- Protein sequences for similarity searching
- Descriptions for text-based queries

### Query Optimization Tips

1. **Use Composite Keys**: When querying proteins, always include both `parent_protein_id` and `parent_tim_barrel_accession` for optimal performance

2. **Leverage Views**: Use `protein_summary` view for common queries instead of complex joins

3. **JSON Queries**: Use `JSON_EXTRACT()` for efficient JSON field queries

4. **Batch Operations**: Use batch inserts for large data collections

## Data Integrity Constraints

### Foreign Key Constraints

1. **interpro_proteins.tim_barrel_accession** → **tim_barrel_entries.accession**
   - Ensures all proteins reference valid TIM barrel entries
   - CASCADE DELETE removes proteins when entries are deleted

2. **proteins.(parent_protein_id, parent_tim_barrel_accession)** → **interpro_proteins.(uniprot_id, tim_barrel_accession)**
   - Ensures all isoforms reference valid protein-entry combinations
   - CASCADE DELETE removes isoforms when parent proteins are deleted

### Data Validation

**Application-Level Validation**:
- Protein sequences contain only valid amino acids (ACDEFGHIKLMNPQRSTVWY)
- TIM barrel coordinates are within sequence bounds
- Accession formats match expected patterns (PF##### or IPR######)
- Entry types are either 'pfam' or 'interpro'

**Database-Level Constraints**:
- NOT NULL constraints on required fields
- VARCHAR length limits prevent oversized data
- DEFAULT values for common fields (organism, timestamps)

## Backup and Recovery

### Backup Strategy

**Full Backup**:
```sql
mysqldump --single-transaction --routines --triggers protein_collector > backup.sql
```

**Table-Specific Backup**:
```sql
mysqldump --single-transaction protein_collector tim_barrel_entries interpro_proteins proteins > data_backup.sql
```

**Data-Only Backup**:
```sql
mysqldump --no-create-info --single-transaction protein_collector > data_only.sql
```

### Recovery Procedures

**Full Recovery**:
```sql
mysql protein_collector < backup.sql
```

**Selective Recovery**:
```sql
mysql protein_collector < data_backup.sql
```

## Monitoring and Maintenance

### Health Checks

**Data Integrity Check**:
```sql
-- Check for orphaned records
SELECT 'Orphaned proteins' as issue, COUNT(*) as count
FROM interpro_proteins ip
LEFT JOIN tim_barrel_entries tbe ON ip.tim_barrel_accession = tbe.accession
WHERE tbe.accession IS NULL

UNION ALL

SELECT 'Orphaned isoforms' as issue, COUNT(*) as count
FROM proteins p
LEFT JOIN interpro_proteins ip ON p.parent_protein_id = ip.uniprot_id 
    AND p.parent_tim_barrel_accession = ip.tim_barrel_accession
WHERE ip.uniprot_id IS NULL;
```

**Performance Monitoring**:
```sql
-- Check index usage
SHOW INDEX FROM interpro_proteins;
SHOW INDEX FROM proteins;

-- Check table sizes
SELECT 
    table_name,
    ROUND(((data_length + index_length) / 1024 / 1024), 2) AS 'Size (MB)'
FROM information_schema.tables
WHERE table_schema = 'protein_collector'
ORDER BY (data_length + index_length) DESC;
```

### Maintenance Tasks

**Regular Maintenance**:
```sql
-- Optimize tables
OPTIMIZE TABLE tim_barrel_entries, interpro_proteins, proteins;

-- Update statistics
ANALYZE TABLE tim_barrel_entries, interpro_proteins, proteins;

-- Check table integrity
CHECK TABLE tim_barrel_entries, interpro_proteins, proteins;
```

## Current Data Statistics

Based on the latest collection (January 2026):

- **TIM Barrel Entries**: 49 total (18 PFAM families + 31 InterPro entries)
- **Human Proteins**: 407 unique proteins collected
- **Protein Isoforms**: Ready for collection (Phase 3)
- **Collection Success Rate**: 100% (0 parsing errors after bug fixes)

## Future Enhancements

### Planned Improvements

1. **Partitioning**: Consider table partitioning for large-scale data
2. **Archiving**: Implement data archiving for historical collections
3. **Replication**: Set up read replicas for query performance
4. **Monitoring**: Enhanced monitoring with performance metrics
5. **Compression**: Evaluate sequence compression for storage efficiency

### Schema Extensions

1. **Protein Interactions**: Table for protein-protein interactions
2. **Structural Data**: Enhanced structural annotation storage
3. **Evolutionary Data**: Phylogenetic and evolutionary information
4. **Expression Data**: Gene expression and tissue specificity data

---

*Last Updated: January 10, 2026*  
*Schema Version: 2.0 (Unified Structure)*  
*Status: Production Ready*