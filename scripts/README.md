# TIM Barrel Collection Scripts

This directory contains the essential scripts for collecting TIM barrel protein family data from InterPro.

## Available Scripts

### 1. `collect_tim_barrel_entries.py` - Main Collection Script
**Purpose**: Unified script to collect both PFAM families and InterPro entries with TIM barrel annotations.

**Features**:
- **Hybrid Search Strategy**: Combines direct PFAM family search with InterPro entry (IPR) search
- **Comprehensive Coverage**: Finds both protein families (PF) and structural classifications (IPR)
- **Unified Storage**: Stores all entries in a single `tim_barrel_entries` table
- **No Artificial Limits**: Collects all available TIM barrel entries found
- **Deduplication**: Automatically removes duplicate entries

**Usage**:
```bash
# Production run
python scripts/collect_tim_barrel_entries.py

# Dry run (test without storing data)
python scripts/collect_tim_barrel_entries.py --dry-run

# Verbose output
python scripts/collect_tim_barrel_entries.py --verbose

# Custom configuration
python scripts/collect_tim_barrel_entries.py --config config/production.json
```

**Current Results**: Successfully collects **49 total entries** (18 PFAM families + 31 InterPro entries)

### 2. `collect_human_proteins.py` - Human Protein Collection Script
**Purpose**: Collect all Homo sapiens proteins for each TIM barrel entry found in step 1.

**Features**:
- **Comprehensive Protein Collection**: Gets proteins for both PFAM families and InterPro entries
- **Human-Specific Filtering**: Focuses on Homo sapiens proteins only
- **Duplicate Detection**: Avoids storing duplicate proteins across entries
- **Progress Tracking**: Detailed logging and statistics
- **Database Integration**: Stores proteins with proper relationships to TIM barrel entries

**Usage**:
```bash
# Production run (collect all human proteins)
python scripts/collect_human_proteins.py

# Dry run (test without storing data)
python scripts/collect_human_proteins.py --dry-run

# Verbose output
python scripts/collect_human_proteins.py --verbose

# Custom configuration
python scripts/collect_human_proteins.py --config config/production.json
```

**Prerequisites**: Must run `collect_tim_barrel_entries.py` first to populate TIM barrel entries.

### 3. `protein_collection_summary.py` - Protein Collection Status
**Purpose**: Display current status of human protein collection for TIM barrel entries.

**Features**:
- Shows collection progress (percentage of entries processed)
- Lists entries with proteins collected and counts
- Shows top entries by protein count
- Displays example collected proteins
- Provides next steps guidance

**Usage**:
```bash
python scripts/protein_collection_summary.py
```

### 4. `tim_barrel_summary.py` - TIM Barrel Entry Status
**Purpose**: Display current status of TIM barrel entry collection from the database.

**Features**:
- Shows count of PFAM families and InterPro entries
- Lists all collected entries with details
- Displays collection strategy information

**Usage**:
```bash
python scripts/tim_barrel_summary.py
```

## Collection Workflow

The complete data collection process follows this sequence:

### Step 1: Collect TIM Barrel Entries
```bash
python scripts/collect_tim_barrel_entries.py
```
- Finds 49 TIM barrel entries (18 PFAM + 31 InterPro)
- Stores entries in `tim_barrel_entries` table

### Step 2: Collect Human Proteins
```bash
python scripts/collect_human_proteins.py
```
- For each TIM barrel entry, finds all Homo sapiens proteins
- Stores proteins in `interpro_proteins` table with relationships

### Step 3: Check Collection Status
```bash
python scripts/protein_collection_summary.py
```
- Shows progress of human protein collection
- Displays statistics and next steps

### Step 4: View TIM Barrel Entry Status
```bash
python scripts/tim_barrel_summary.py
```
- Shows current collection status and statistics

## Collection Strategy

### TIM Barrel Entry Collection (Step 1)
The unified collection script uses a **comprehensive hybrid search approach**:

#### Phase 1: Direct PFAM Family Search
- Searches for PFAM families using multiple TIM barrel-related terms
- Terms include: "TIM barrel", "TIM-barrel", "triosephosphate", "aldolase", etc.
- Filters results to ensure they are actually TIM barrel-related

#### Phase 2: InterPro Entry Search  
- Searches for InterPro entries (IPR records) with TIM barrel annotations
- Captures structural classifications like IPR013785 (Aldolase-type TIM barrel)
- Includes domain families, homologous superfamilies, and active sites

#### Unified Storage
- All entries stored in single `tim_barrel_entries` table
- Each entry has `entry_type` field: 'pfam' or 'interpro'
- Automatic deduplication by accession number
- Comprehensive metadata preservation

### Human Protein Collection (Step 2)
The protein collection script processes each TIM barrel entry:

#### PFAM Family Processing
- Uses InterPro API endpoint: `/protein/UniProt/entry/pfam/{accession}/`
- Filters by `tax_lineage=Homo sapiens`
- Collects all human proteins belonging to each PFAM family

#### InterPro Entry Processing
- Uses InterPro API endpoint: `/protein/UniProt/entry/interpro/{accession}/`
- Filters by `tax_lineage=Homo sapiens`
- Collects all human proteins annotated with each InterPro entry

#### Protein Storage
- Stores proteins in `interpro_proteins` table
- Links each protein to its TIM barrel entry via `tim_barrel_accession`
- Handles both PFAM and InterPro relationships uniformly
- Automatic duplicate detection across entries

## Database Schema

### TIM Barrel Entries Table
The unified `tim_barrel_entries` table contains:
- `accession`: Primary identifier (PF##### or IPR######)
- `entry_type`: 'pfam' or 'interpro'
- `name`: Entry name
- `description`: Detailed description
- `interpro_type`: Type for InterPro entries (domain, family, etc.)
- `tim_barrel_annotation`: TIM barrel-specific annotation
- `member_databases`: Associated databases (for InterPro entries)
- `interpro_id`: Associated InterPro ID (for PFAM entries)

### Human Proteins Table
The `interpro_proteins` table contains:
- `uniprot_id`: UniProt protein identifier (primary key)
- `tim_barrel_accession`: Foreign key to TIM barrel entry
- `name`: Protein name
- `organism`: Source organism (Homo sapiens)
- `created_at`: Record creation timestamp

### Relationships
- One-to-many: Each TIM barrel entry can have multiple proteins
- Foreign key constraint ensures data integrity
- Cascade delete removes proteins when TIM barrel entry is deleted

## Migration History

This directory previously contained multiple scripts that have been consolidated:

**Removed Scripts** (functionality now in unified scripts):
- `collect_pfam_families.py` - Old separate PFAM collection
- `manual_tim_barrel_additions.py` - Manual addition tool
- `add_interpro_entries_table.py` - Database migration script
- `migrate_to_unified_table.py` - Table consolidation script

The current setup represents the final, clean state after successful migration to the unified approach.