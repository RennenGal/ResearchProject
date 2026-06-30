# TIM Barrel Alternative Splicing — Data Collection & Analysis

Research tool for studying the effect of alternative splicing (AS) on TIM barrel domain structure.
Collects domain family entries, proteins, and isoforms from InterPro and UniProt into a local
SQLite database; annotates the (βα)₈ repeat motifs; maps splice events onto the motif structure.

**Institution**: Ben-Gurion University of the Negev  
**Supervisors**: Prof. Tal Shay, Prof. Chen Keasar

---

## Quick start

```bash
git clone https://github.com/RennenGal/ResearchProject.git
cd ResearchProject
python -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate
pip install -r requirements.txt

# Run the full pipeline (Phases 1–13)
python scripts/collect.py
```

---

## Data collection pipeline

A single command runs all 13 pipeline phases end-to-end:

```bash
python scripts/collect.py [--domain tim_barrel] [--organism homo_sapiens] [--db PATH]
```

| Phase | Step |
|-------|------|
| 1 | Fetch InterPro / UniProt domain family entries |
| 2 | Collect protein records from UniProt |
| 3 | Collect canonical and alternative isoforms from UniProt |
| 4 | Backfill protein metadata (name, reviewed status, annotation score) |
| 5 | Fetch and propagate gene names |
| 6 | Backfill domain locations from InterPro |
| 7 | Gene-level deduplication (keep best representative per gene) |
| 8 | Build `affected_isoforms` table (VSP-based AS detection) + fragment isoforms |
| 9 | Build `canonical_analysis` table |
| 10 | Annotate (βα)₈ motifs via AlphaFold + DSSP |
| 11 | Collect Ensembl transcripts + alignment analysis + exon boundary data |
| 12 | Backfill isoform exon junction data (UniProt isoforms) |
| 13 | Build final `analysis_proteins` table and views |

### Partial / maintenance modes

```bash
# Resume isoform collection for new proteins only
python scripts/collect.py --resume

# Update entries and proteins only (no isoform changes)
python scripts/collect.py --collect-proteins

# Re-fetch all isoforms from scratch
python scripts/collect.py --recollect-isoforms

# Re-query domain locations for canonical isoforms
python scripts/collect.py --backfill-domains
```

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--domain` | `tim_barrel` | Domain type to collect |
| `--organism` | `homo_sapiens` | Organism to collect |
| `--db` | `db/protein_data.db` | Database path |
| `--log-file` | — | Also write logs to file |
| `--log-level` | `INFO` | Logging verbosity |

### Adding a new domain

Add an entry to `DOMAINS` in `protein_data_collector/config.py`:

```python
"my_domain": DomainConfig(
    display_name="My Domain",
    interpro_annotation="...",
    entries_table="entries",
    table_prefix="",
    accession_col="my_domain_accession",
    location_col="my_domain_location",
    sequence_col="my_domain_sequence",
    cathgene3d_search="...",
)
```

Then run `python scripts/collect.py --domain my_domain`.

---

## Analysis

After collection, generate figures and results:

```bash
python scripts/run_alternative_analysis.py   # all figures → figures/alt/
python scripts/write_results_md.py           # regenerate Results.md
```

See `Results.md` for findings and `Statistical-Framework.md` for the formal statistical definitions.

---

## Database schema

### `entries`
```
accession (PK), entry_type  -- 'pfam' | 'interpro' | 'cathgene3d'
name, description, domain_annotation
```

### `proteins`
```
uniprot_id (PK), {domain}_accession (FK → entries)
protein_name, gene_name, organism, reviewed
protein_existence, annotation_score
canonical_uniprot_id   -- NULL = canonical; non-NULL = redundant, points to canonical
```

### `isoforms`
```
isoform_id (PK), uniprot_id (FK), is_canonical, is_fragment
sequence, sequence_length
exon_count, exon_annotations   -- JSON list of 1-based exon boundary positions
splice_variants                -- JSON UniProt VSP features
{domain}_location              -- JSON {domain_id, start, end, length, source}
{domain}_sequence              -- domain subsequence
ensembl_transcript_id, alphafold_id
```

*For `tim_barrel` the columns are `tim_barrel_location` and `tim_barrel_sequence`.*

### `affected_isoforms`
```
isoform_id (PK), uniprot_id (FK)
sequence, sequence_length, exon_annotations, splice_variants
domain_location, domain_sequence
canonical_domain_location, canonical_domain_sequence
identity_percentage, alignment_score
exon_boundary_in_domain, exon_boundaries_in_domain_count
vsp_domain_events     -- JSON VSPs overlapping the domain
detection_method      -- 'vsp_overlap' | 'sliding_window_fallback'
```

### `canonical_analysis`
```
uniprot_id (PK), domain_index, gene_name
sequence, domain_start, domain_end
exon_annotations    -- [{exon, start, end}] 1-based inclusive
motif_annotations   -- [{motif, beta_start, beta_end, alpha_start, alpha_end}]
dssp_source, pdb_source
```

### `analysis_proteins` + views
```
analysis_proteins   -- canonical and affected isoforms used in the analysis
view_canonical      -- canonical proteins with motif annotation
view_noncanonical   -- AS-affected isoforms with VSP domain events
```

### Ensembl tables
```
ensembl_transcripts  -- transcript sequences and exon annotations per gene
ensembl_affected     -- AS-affected Ensembl transcripts with domain alignment
```

---

## Querying

```python
from protein_data_collector.query.engine import QueryEngine

q = QueryEngine()
print(q.summary())

isoforms     = q.get_isoforms_for_protein("Q12794")   # HYAL1
alt_proteins = q.get_proteins_with_alternative_isoforms()
```

```python
from protein_data_collector.query.export import to_fasta, to_csv

fasta = to_fasta(q.get_all_isoforms())
csv   = to_csv(q.get_all_isoforms())
```

---

## Tests

```bash
pytest tests/ -v
```

| File | Coverage |
|------|----------|
| `test_models.py` | Pydantic model validation — sequence, length, bounds, fragment flag, domain sequence slicing |
| `test_database.py` | CRUD round-trips, FK enforcement, resume detection |
| `test_uniprot_parser.py` | Parser functions against a real UniProt response |
| `test_query.py` | QueryEngine queries, FASTA/CSV/JSON export |
| `test_motif_annotator.py` | `identify_ba_motifs()` and `_reformat_exon_annotations()` |

---

## Architecture

```
protein_data_collector/
  api/
    interpro_client.py    InterPro REST API
    uniprot_client.py     UniProt REST API
    ensembl_client.py     Ensembl REST API
  collector/
    interpro_collector.py   Phases 1–2: domain families and proteins
    uniprot_collector.py    Phase 3: isoforms and splice variants
    data_collector.py       Pipeline orchestrator
  database/
    schema.py     DDL — domain-parameterised table definitions
    connection.py get_connection() context manager
    storage.py    upsert_* / get_* CRUD functions
  models/
    entities.py   Pydantic models: DomainEntry, Protein, Isoform
  query/
    engine.py     QueryEngine
    export.py     to_fasta(), to_csv(), to_json()
  analysis/
    motif_annotator.py        identify_ba_motifs(): βα repeat detection from DSSP
    tim_barrel_alignment.py   populate_tim_barrel_isoforms(): VSP-based AS detection
  config.py   DomainConfig + OrganismConfig; DOMAINS and ORGANISMS registries

scripts/
  collect.py                  Pipeline entry point (Phases 1–13)

  Pipeline modules (called by collect.py — not standalone):
    backfill_protein_metadata.py
    fetch_gene_names.py
    backfill_domain_locations.py
    dedup_by_gene.py
    build_affected_isoforms.py
    build_canonical_analysis.py
    annotate_motifs.py
    collect_ensembl.py
    backfill_isoform_exons.py
    create_analysis_table.py

  Analysis scripts (run after collection):
    run_alternative_analysis.py   Generate all figures → figures/alt/
    write_results_md.py           Regenerate Results.md

tests/
  conftest.py
  test_models.py
  test_database.py
  test_uniprot_parser.py
  test_query.py
  test_motif_annotator.py
```
