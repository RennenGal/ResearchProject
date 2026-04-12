# TIM Barrel Alternative Splicing — Data Collector

Research tool for studying the effect of alternative splicing on TIM barrel domain structure.
Collects protein families, human proteins, and isoforms from InterPro and UniProt into a local SQLite database.

**Institution**: Ben-Gurion University of the Negev  
**Supervisors**: Prof. Tal Shay, Prof. Chen Keasar

---

## Database contents (current)

| | Count |
|---|---|
| TIM barrel families (18 PFAM + 34 InterPro) | 52 |
| Human proteins | 572 |
| Isoforms — canonical | 572 |
| Isoforms — alternative | 131 |
| **Isoforms — total** | **703** |
| Fragments (sequence < 200 aa) | 185 |
| With `tim_barrel_location` | 572 (all canonical) |
| With `tim_barrel_sequence` | 404 (non-fragment canonical) |

Families were identified by name-based search (`TIM barrel`, `TIM-barrel`) combined with
structural database cross-checking (CATH 3.20.20 superfamily) to avoid missing entries
whose names don't contain "TIM".

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
  sequence, sequence_length,
  is_fragment,          -- 1 if sequence_length < 200 (UniProt fragment entry, not a full protein)
  exon_count,
  exon_annotations,     -- JSON; null (future: Ensembl coordinate mapping)
  splice_variants,      -- JSON; UniProt VSP features per isoform
  tim_barrel_location,  -- JSON; {domain_id, start, end, length, source}
                        --   populated for all canonical isoforms via InterPro
                        --   null for alternative isoforms (needs splice coordinate mapping)
  tim_barrel_sequence,  -- subsequence sequence[start-1:end] from tim_barrel_location
                        --   null for fragments and alternative isoforms
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

**Backfill domain locations** — populate `tim_barrel_location` for canonical isoforms where it is NULL:
```bash
python scripts/collect.py --backfill-domains
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
pytest tests/ -v          # 58 tests
pytest tests/ -v --tb=short
```

Test files:
- `tests/test_models.py` — Pydantic model validation (sequence, length, bounds, fragment flag, TIM barrel sequence slicing)
- `tests/test_database.py` — CRUD round-trips, FK enforcement, new fields, resume detection
- `tests/test_uniprot_parser.py` — Pure parser functions against a real UniProt response
- `tests/test_query.py` — QueryEngine queries, FASTA/CSV export format

---

## Known gaps

1. **`tim_barrel_location` for alternative isoforms is null.**  
   Alternative isoforms are structurally distinct from the canonical sequence due to AS events.
   Mapping their domain boundaries requires projecting splice variant coordinates onto the
   altered sequence — planned as a future analysis step.

2. **`tim_barrel_sequence` for alternative isoforms is null.**  
   Follows from gap 1 — once domain location is determined, slicing is trivial.

3. **`exon_annotations` is null for all isoforms.**  
   Mapping genomic exon coordinates to protein coordinates requires the Ensembl REST API.
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
    data_collector.py       Pipeline orchestrator (CollectionReport, backfill)
  database/
    schema.py     SQL DDL constants, init_db()
    connection.py get_connection() context manager, ensure_db()
    storage.py    upsert_*/get_* CRUD functions
  models/
    entities.py   Pydantic models: TIMBarrelEntry, Protein, Isoform
                  Isoform auto-computes is_fragment and tim_barrel_sequence
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
