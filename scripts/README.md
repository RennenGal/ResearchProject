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

### 2. `tim_barrel_summary.py` - Status Summary
**Purpose**: Display current status of TIM barrel collection from the database.

**Features**:
- Shows count of PFAM families and InterPro entries
- Lists all collected entries with details
- Displays collection strategy information

**Usage**:
```bash
python scripts/tim_barrel_summary.py
```

## Collection Strategy

The unified collection script uses a **comprehensive hybrid search approach**:

### Phase 1: Direct PFAM Family Search
- Searches for PFAM families using multiple TIM barrel-related terms
- Terms include: "TIM barrel", "TIM-barrel", "triosephosphate", "aldolase", etc.
- Filters results to ensure they are actually TIM barrel-related

### Phase 2: InterPro Entry Search  
- Searches for InterPro entries (IPR records) with TIM barrel annotations
- Captures structural classifications like IPR013785 (Aldolase-type TIM barrel)
- Includes domain families, homologous superfamilies, and active sites

### Unified Storage
- All entries stored in single `tim_barrel_entries` table
- Each entry has `entry_type` field: 'pfam' or 'interpro'
- Automatic deduplication by accession number
- Comprehensive metadata preservation

## Database Schema

The unified `tim_barrel_entries` table contains:
- `accession`: Primary identifier (PF##### or IPR######)
- `entry_type`: 'pfam' or 'interpro'
- `name`: Entry name
- `description`: Detailed description
- `interpro_type`: Type for InterPro entries (domain, family, etc.)
- `tim_barrel_annotation`: TIM barrel-specific annotation
- `member_databases`: Associated databases (for InterPro entries)
- `interpro_id`: Associated InterPro ID (for PFAM entries)

## Migration History

This directory previously contained multiple scripts that have been consolidated:

**Removed Scripts** (functionality now in unified script):
- `collect_pfam_families.py` - Old separate PFAM collection
- `manual_tim_barrel_additions.py` - Manual addition tool
- `add_interpro_entries_table.py` - Database migration script
- `migrate_to_unified_table.py` - Table consolidation script

The current setup represents the final, clean state after successful migration to the unified approach.