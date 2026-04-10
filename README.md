# TIM Barrel Alternative Splicing — Data Collector

Research tool for studying the effect of alternative splicing on TIM barrel domain structure.
Collects protein families, human proteins, and isoforms from InterPro and UniProt into a local SQLite database.

**Institution**: Ben-Gurion University of the Negev  
**Supervisors**: Prof. Tal Shay, Prof. Chen Keasar

---

## Database contents (current)

| Table | Count |
|-------|-------|
| TIM barrel families (PFAM + InterPro) | 49 |
| Human proteins | 407 |
| Isoforms — canonical | 407 |
| Isoforms — alternative | 90 |
| **Isoforms — total** | **497** |

---

## Schema

Three tables with cascade foreign keys:

```
tim_barrel_entries
  accession (PK), entry_type, name, description, tim_barrel_annotation

proteins
  uniprot_id (PK), tim_barrel_accession (FK), protein_name, gene_name,
  organism, reviewed, protein_existence, annotation_score

isoforms
  isoform_id (PK), uniprot_id (FK), is_canonical,
  sequence, sequence_length, exon_count,
  exon_annotations,      -- JSON; null (future: Ensembl coordinate mapping)
  splice_variants,       -- JSON; UniProt VSP features per isoform
  tim_barrel_location,   -- JSON; null (known gap — see below)
  ensembl_gene_id, alphafold_id
```

---

## Setup

```bash
git clone https://github.com/RennenGal/ResearchProject.git
cd ResearchProject
python -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

---

## Data collection

**Full collection from scratch** (all three phases):
```bash
python scripts/collect.py
```

**Resume** — collect isoforms only for proteins not yet in the database:
```bash
python scripts/collect.py --resume
```

**Re-collect all isoforms** — wipe isoform table and re-fetch from UniProt:
```bash
python scripts/collect.py --recollect-isoforms
```

**Options**:
```
--db PATH         Override database path (default: db/protein_data.db)
--log-file PATH   Write logs to file in addition to stdout
--log-level LVL   Logging verbosity (default: INFO)
```

---

## Migrate from old schema

If you have a database in the pre-2025 schema (single `proteins` table with JSON blobs):
```bash
python scripts/migrate.py --old db/protein_data.db --new db/protein_data_v2.db
```

---

## Querying the database

```python
from protein_data_collector.query.engine import QueryEngine

q = QueryEngine()
print(q.summary())

isoforms = q.get_isoforms_for_protein("P04637")   # TP53
alt_proteins = q.get_proteins_with_alternative_isoforms()
```

**Export**:
```python
from protein_data_collector.query.export import to_fasta, to_csv

fasta = to_fasta(q.get_all_isoforms())
csv   = to_csv(q.get_all_isoforms())
```

---

## Tests

```bash
pytest tests/ -v          # 51 tests
pytest tests/ -v --tb=short
```

Test files:
- `tests/test_models.py` — Pydantic model validation (sequence characters, length, bounds)
- `tests/test_database.py` — CRUD round-trips, FK enforcement, resume detection
- `tests/test_uniprot_parser.py` — Pure parser functions against a real UniProt response
- `tests/test_query.py` — QueryEngine queries, FASTA/CSV export format

---

## Known gaps

1. **`tim_barrel_location` is null** for all isoforms.  
   Cause: `InterProClient.get_domain_boundaries()` filters by name containing "TIM barrel", but
   InterPro domain names are like "Aldolase_TIM".  
   Fix needed: match by accession from `tim_barrel_entries` instead of name string.

2. **`exon_annotations` is null** for all isoforms.  
   Cause: mapping genomic exon coordinates to protein coordinates requires the Ensembl REST API.  
   Planned as a future collection step.

---

## Architecture

```
protein_data_collector/
  api/
    interpro_client.py   InterPro REST API (pagination, domain boundaries)
    uniprot_client.py    UniProt REST API (protein JSON, isoform FASTA)
  collector/
    interpro_collector.py   Phase 1+2: families and human proteins
    uniprot_collector.py    Phase 3: isoforms and splice variant extraction
    data_collector.py       Pipeline orchestrator (CollectionReport)
  database/
    schema.py     SQL DDL constants, init_db()
    connection.py get_connection() context manager, ensure_db()
    storage.py    upsert_*/get_* CRUD functions
  models/
    entities.py   Pydantic models: TIMBarrelEntry, Protein, Isoform
  query/
    engine.py     QueryEngine — SQL-backed queries
    export.py     to_fasta(), to_csv(), to_json()
  config.py       Config dataclass (db_path, API URLs, retry settings)
  errors.py       Exception hierarchy
  retry.py        tenacity decorator

scripts/
  collect.py    Data collection entry point
  migrate.py    One-time schema migration from old database
```
