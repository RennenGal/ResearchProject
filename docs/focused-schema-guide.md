# Focused UniProt Schema Guide

## 🎯 Overview

This focused schema provides a manageable approach to UniProt data collection, containing **67 essential fields** across **9 selected categories** specifically chosen for TIM barrel protein research.

## 📊 Selected Categories

| Category | Fields | Description | TIM Barrel Relevance |
|----------|--------|-------------|---------------------|
| **Names & Taxonomy** | 14 | Protein identifiers, gene names, organism info | Essential for identification |
| **Sequences** | 19 | Protein sequences, variants, mass | Core data for analysis |
| **Function** | 16 | Catalytic activity, pathways, binding sites | Functional classification |
| **Interaction** | 2 | Protein-protein interactions, subunit structure | Complex formation |
| **Gene Ontology** | 5 | GO terms for biological processes, functions | Standardized annotations |
| **Structure** | 4 | Secondary structure elements (helix, strand, turn) | **Critical for TIM barrels** |
| **Date Information** | 4 | Creation, modification dates, versions | Data provenance |
| **Family & Domains** | 9 | Protein families, domains, motifs, repeats | Structural classification |
| **3D Structure Databases** | 7 | PDB, AlphaFold, and other structure references | Experimental validation |

**Total: 67 fields + 3 quality indicators + 1 category field = 71 columns**

## 🗂️ Category Field for Filtering

The schema includes a `data_category` ENUM field that allows easy filtering by data type:

```sql
-- Filter by category
SELECT * FROM proteins_focused WHERE data_category = 'structure';

-- Get proteins with structural data
SELECT * FROM proteins_focused 
WHERE data_category IN ('structure', '3d_structure_databases');

-- Category-based statistics
SELECT data_category, COUNT(*) as protein_count 
FROM proteins_focused 
GROUP BY data_category;
```

### Category Values:
- `names_taxonomy`
- `sequences` (default)
- `function`
- `interaction`
- `gene_ontology`
- `structure`
- `dates`
- `family_domains`
- `3d_structure_databases`

## 🔧 Implementation

### 1. Database Migration
```bash
# Run the focused migration script
python scripts/migrate_to_focused_schema.py
```

### 2. Field Management
```python
# Use the focused field management module
from scripts.focused_uniprot_fields import get_fields_for_api_request

# Get all focused fields (67 fields)
all_fields = get_fields_for_api_request()

# Get specific categories
structural_fields = get_fields_for_api_request(['structure', '3d_structure_databases'])
functional_fields = get_fields_for_api_request(['function', 'gene_ontology'])

# Get essential TIM barrel fields only
from scripts.focused_uniprot_fields import ESSENTIAL_TIM_BARREL_FIELDS
essential_fields = ','.join(ESSENTIAL_TIM_BARREL_FIELDS)
```

### 3. API Collection Examples
```python
# Example collection with focused fields
import requests
from scripts.focused_uniprot_fields import get_fields_for_api_request

def collect_protein_data(uniprot_id, categories=None):
    """Collect protein data with focused fields."""
    fields = get_fields_for_api_request(categories)
    
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
    params = {
        'fields': fields,
        'format': 'json'
    }
    
    response = requests.get(url, params=params)
    return response.json()

# Collect structural data only
structural_data = collect_protein_data('P60174', ['structure', '3d_structure_databases'])

# Collect functional data only  
functional_data = collect_protein_data('P60174', ['function', 'gene_ontology'])

# Collect all focused data
complete_data = collect_protein_data('P60174')
```

## 📈 Research Applications

### 1. TIM Barrel Structure Analysis
```sql
-- Find proteins with secondary structure annotations
SELECT 
    isoform_id,
    protein_name,
    gene_primary,
    CASE WHEN helix IS NOT NULL THEN 'Yes' ELSE 'No' END as has_helix_data,
    CASE WHEN beta_strand IS NOT NULL THEN 'Yes' ELSE 'No' END as has_strand_data,
    CASE WHEN turn IS NOT NULL THEN 'Yes' ELSE 'No' END as has_turn_data
FROM proteins_focused 
WHERE tim_barrel_location IS NOT NULL
  AND (helix IS NOT NULL OR beta_strand IS NOT NULL OR turn IS NOT NULL);
```

### 2. Functional Classification
```sql
-- Analyze TIM barrel proteins by function
SELECT 
    ec_number,
    COUNT(*) as protein_count,
    AVG(sequence_length) as avg_length,
    GROUP_CONCAT(DISTINCT gene_primary) as genes
FROM proteins_focused 
WHERE tim_barrel_location IS NOT NULL 
  AND ec_number IS NOT NULL
GROUP BY ec_number
ORDER BY protein_count DESC;
```

### 3. Structure-Function Relationships
```sql
-- Compare proteins with experimental vs predicted structures
SELECT 
    CASE 
        WHEN xref_pdb IS NOT NULL THEN 'Experimental (PDB)'
        WHEN xref_alphafolddb IS NOT NULL THEN 'Predicted (AlphaFold)'
        ELSE 'No Structure'
    END as structure_type,
    COUNT(*) as protein_count,
    COUNT(CASE WHEN ec_number IS NOT NULL THEN 1 END) as with_ec_number,
    AVG(annotation_score) as avg_annotation_score
FROM proteins_focused 
WHERE tim_barrel_location IS NOT NULL
GROUP BY structure_type;
```

### 4. Quality Assessment
```sql
-- Filter by data quality
SELECT 
    protein_existence,
    reviewed,
    COUNT(*) as protein_count,
    AVG(annotation_score) as avg_score
FROM proteins_focused 
WHERE tim_barrel_location IS NOT NULL
GROUP BY protein_existence, reviewed
ORDER BY protein_count DESC;
```

## 🎯 Advantages of Focused Schema

### **Manageable Size**
- **67 fields** vs 282 in full schema
- Faster queries and reduced storage
- Easier to understand and maintain

### **Targeted for TIM Barrels**
- All essential structural fields included
- Functional classification capabilities
- Quality assessment built-in

### **Category-Based Organization**
- Easy filtering by data type
- Logical grouping of related fields
- Flexible collection strategies

### **Research-Ready**
- Direct support for structure-function analysis
- Quality indicators for data filtering
- Cross-reference integration

## 📊 Field Population Expectations

| Category | Expected Population | Notes |
|----------|-------------------|-------|
| **Names & Taxonomy** | 95-100% | Core identifiers always present |
| **Sequences** | 100% | Required for all proteins |
| **Function** | 60-80% | Higher for reviewed entries |
| **Structure** | 30-50% | Depends on experimental data |
| **Gene Ontology** | 70-90% | Well-annotated for human proteins |
| **3D Structure Databases** | 20-30% PDB, 80-90% AlphaFold | Growing coverage |
| **Family & Domains** | 80-95% | Good coverage for domain databases |

## 🚀 Next Steps

1. **Run Migration**: Execute `scripts/migrate_to_focused_schema.py`
2. **Update Collection**: Modify scripts to use `scripts/focused_uniprot_fields.py`
3. **Test Categories**: Validate data collection for each category
4. **Analyze Results**: Use category-based queries for research
5. **Expand Gradually**: Add more categories if needed

## 📁 Related Files

- **Schema**: `scripts/temp/focused_schema.sql`
- **Migration**: `scripts/migrate_to_focused_schema.py`
- **Field Management**: `scripts/focused_uniprot_fields.py`
- **Documentation**: `docs/focused-schema-guide.md`

This focused approach provides a perfect balance between comprehensive data collection and practical usability for TIM barrel protein research.