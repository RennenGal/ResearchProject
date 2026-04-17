# Protein Domain Alternative Splicing — Data Collector

Research tool for studying the effect of alternative splicing on protein domain structure.
Collects domain family entries, proteins, and isoforms from InterPro and UniProt into a
local SQLite database. Supports multiple domain families and organisms.

**Institution**: Ben-Gurion University of the Negev  
**Supervisors**: Prof. Tal Shay, Prof. Chen Keasar

---

## Database contents (current)

### TIM barrel (Homo sapiens)

| | Count |
|---|---|
| TIM barrel families (38 Pfam/InterPro + 35 Gene3D CATH) | 73 |
| Proteins — canonical | 1,174 |
| Isoforms — canonical | 1,174 |
| Isoforms — alternative | 249 |
| **Isoforms — total** | **1,423** |
| AS-affected isoforms (domain disrupted) | 37 |

### Beta propeller (Homo sapiens)

| | Count |
|---|---|
| Beta propeller families (148 Pfam/InterPro + 9 Gene3D CATH) | 157 |
| Proteins — canonical | 2,853 |
| Isoforms — canonical | 2,853 |
| Isoforms — alternative | 712 |
| **Isoforms — total** | **3,565** |
| AS-affected isoforms (domain disrupted) | 143 |

### Additional organisms (TIM barrel only)

| Organism | Proteins | Isoforms | AS-affected |
|---|---|---|---|
| Mus musculus | 325 | 360 | 4 |
| Rattus norvegicus | 507 | 515 | 3 |

Proteins are deduplicated by `(protein_name, organism)` group: entries that share a name
and organism with a better-annotated representative are marked redundant via
`proteins.canonical_uniprot_id`. Isoforms are collected only for canonical proteins.

---

## Schema

Tables follow a `{domain_prefix}_{table}{organism_suffix}` naming convention:
- Domain prefixes: `tb_` (TIM barrel), `bp_` (beta propeller)
- Organism suffixes: `` (Homo sapiens), `_mus_musculus`, `_rattus_norvegicus`

### Domain entries

```
tb_entries / bp_entries
  accession (PK), entry_type  -- 'pfam' | 'interpro' | 'cathgene3d'
  name, description, domain_annotation
```

### Proteins

```
tb_proteins / bp_proteins [/ tb_proteins_mus_musculus / ...]
  uniprot_id (PK), tim_barrel_accession (FK → entries),
  protein_name, gene_name, organism, reviewed,
  protein_existence, annotation_score,
  canonical_uniprot_id  -- NULL = canonical representative
                        -- non-NULL = redundant, points to canonical
```

### Isoforms

```
tb_isoforms / bp_isoforms [/ tb_isoforms_mus_musculus / ...]
  isoform_id (PK), uniprot_id (FK → proteins), is_canonical,
  sequence, sequence_length,
  is_fragment,          -- 1 if sequence_length < 200 aa
  exon_count, exon_annotations,   -- exon_annotations: null (future: Ensembl mapping)
  splice_variants,      -- JSON; UniProt VSP features per isoform
  tim_barrel_location,  -- JSON {domain_id, start, end, length, source}
                        --   populated for canonical isoforms via InterPro
  tim_barrel_sequence,  -- subsequence sequence[start-1:end]
  ensembl_gene_id, alphafold_id
```

### AS-affected isoforms

```
tb_affected_isoforms / bp_affected_isoforms [/ tb_affected_isoforms_mus_musculus / ...]
  id (PK), isoform_id, uniprot_id,
  domain_location,           -- JSON {start, end} of domain in the alternative isoform
  domain_sequence,           -- domain subsequence in the alternative isoform
  canonical_domain_location, -- JSON {start, end} in the canonical isoform
  canonical_domain_sequence, -- domain subsequence in the canonical
  alignment_identity,        -- sequence identity between the two domain sequences
  insertion_detected         -- 1 if an insertion relative to canonical was detected
```

DB triggers (`trg_block_redundant_tb`, `trg_block_redundant_bp`, etc.) prevent inserting
isoforms for redundant proteins.

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

The `--domain` flag selects the protein domain (`tim_barrel` or `beta_propeller`).
The `--organism` flag selects the organism (`homo_sapiens`, `mus_musculus`, `rattus_norvegicus`).
Both default to `tim_barrel` / `homo_sapiens`.

**Full collection from scratch** (all three phases):
```bash
python scripts/collect.py --domain tim_barrel --organism homo_sapiens
python scripts/collect.py --domain beta_propeller --organism homo_sapiens
```

**Resume** — collect isoforms only for proteins not yet in the database:
```bash
python scripts/collect.py --resume --domain tim_barrel
```

**Re-collect all isoforms** — wipe isoform table and re-fetch from UniProt:
```bash
python scripts/collect.py --recollect-isoforms --domain beta_propeller
```

**Backfill domain locations** — populate `tim_barrel_location` for canonical isoforms where it is NULL:
```bash
python scripts/collect.py --backfill-domains --domain tim_barrel
```

**Update entries and proteins only** (Phase 1+2, no isoform changes):
```bash
python scripts/collect.py --collect-proteins --domain tim_barrel
```

**Options**:
```
--domain DOMAIN     Domain to collect (default: tim_barrel)
--organism ORG      Organism to collect (default: homo_sapiens)
--db PATH           Override database path (default: db/protein_data.db)
--log-file PATH     Write logs to file in addition to stdout
--log-level LVL     Logging verbosity (default: INFO)
```

---

## AS-affected isoform analysis

After collecting isoforms, identify alternative isoforms where alternative splicing disrupts
the domain sequence using a sliding-window local alignment:

```bash
python scripts/build_affected_isoforms.py --domain tim_barrel --organism homo_sapiens
python scripts/build_affected_isoforms.py --domain beta_propeller --organism homo_sapiens
```

Results are stored in the corresponding `*_affected_isoforms` table.

---

## Querying the database

```python
from protein_data_collector.query.engine import QueryEngine

# Default: TIM barrel, Homo sapiens
q = QueryEngine()
print(q.summary())

# Or specify domain / organism explicitly
q = QueryEngine(domain="beta_propeller", organism="homo_sapiens")

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
pytest tests/ -v
pytest tests/ -v --tb=short
```

Test files:
- `tests/test_models.py` — Pydantic model validation (sequence, length, bounds, fragment flag, domain sequence slicing)
- `tests/test_database.py` — CRUD round-trips, FK enforcement, new fields, resume detection
- `tests/test_uniprot_parser.py` — Pure parser functions against a real UniProt response
- `tests/test_query.py` — QueryEngine queries, FASTA/CSV/JSON export format

---

## Known gaps

1. **`tim_barrel_location` for alternative isoforms is null.**  
   Alternative isoforms are structurally distinct from the canonical sequence due to AS events.
   Mapping their domain boundaries requires projecting splice variant coordinates onto the
   altered sequence — planned as a future analysis step.

2. **`exon_annotations` is null for all isoforms.**  
   Mapping genomic exon coordinates to protein coordinates requires the Ensembl REST API.
   Planned as a future collection step.

---

## Architecture

```
protein_data_collector/
  api/
    interpro_client.py   InterPro REST API (pagination, multi-strategy entry search, domain boundaries)
    uniprot_client.py    UniProt REST API (protein JSON, isoform FASTA)
  collector/
    interpro_collector.py   Phase 1+2: domain families and proteins
                            collect_domain_entries() uses four strategies:
                              1. annotation= (InterPro/Pfam exact match)
                              2. search= (text fallback)
                              3. cathgene3d search= (CATH structurally-classified entries)
                              4. extra_accessions (always-included explicit accessions)
    uniprot_collector.py    Phase 3: isoforms and splice variant extraction
    data_collector.py       Pipeline orchestrator (CollectionReport, backfill, deduplication)
  database/
    schema.py     SQL DDL — 20 tables for 2 domains x 3 organisms, init_db()
    connection.py get_connection() context manager, ensure_db()
    storage.py    upsert_*/get_* CRUD functions (parameterised table names)
  models/
    entities.py   Pydantic models: TIMBarrelEntry, Protein, Isoform
  query/
    engine.py     QueryEngine — SQL-backed queries (domain + organism aware)
    export.py     to_fasta(), to_csv(), to_json()
  config.py       DomainConfig + OrganismConfig registries; DOMAINS and ORGANISMS dicts
  errors.py       Exception hierarchy
  retry.py        tenacity decorator

scripts/
  collect.py                Data collection entry point (--domain, --organism)
  build_affected_isoforms.py  AS-affected isoform detection and storage
```

### Domain + organism parameterisation

`DomainConfig` holds the InterPro query strategy for a domain: annotation term, text search
fallback, `cathgene3d_search` for structurally-classified CATH entries (e.g. `"3.20.20"` for
TIM barrel, `"2.130"` for beta propeller), and explicit `extra_accessions` for families with
non-standard names (e.g. WD40 superfamilies for beta propeller).
`OrganismConfig` holds the NCBI taxon ID and derives table names via
`protein_table(domain)`, `isoform_table(domain)`, and `affected_isoforms_table(domain)`.
