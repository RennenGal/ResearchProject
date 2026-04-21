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

Beta propeller is only collected for Homo sapiens. Mouse and rat are supported for TIM barrel only.

### Ensembl transcript expansion (TIM barrel, Homo sapiens)

| | Count |
|---|---|
| Proteins with Ensembl mapping | 799 |
| Transcripts collected | 2,714 |
| Duplicates (sequence already in UniProt isoforms) | 1,490 |
| Fragments (< 200 aa) | 355 |
| Novel transcripts | 1,224 |
| **AS-affected novel transcripts** | **391** |

### Exon boundary analysis (AS-affected Ensembl transcripts)

| | Count |
|---|---|
| AS-affected transcripts with exon data | 391 |
| With ≥1 exon junction inside domain | 391 (100%) |
| Average exon junctions inside domain | 7.1 |
| Transcripts with exactly 1 junction in domain | 19 |

All 391 AS-affected transcripts have at least one exon junction inside the TIM barrel domain
(average 7.1 per transcript). This is expected given domain lengths of ~280 aa and 8–15 coding
exon boundaries per gene. Transcripts with fewer intra-domain junctions (e.g. exactly 1) are
candidates for highly targeted splice events.

Results stored in `tb_ensembl_transcripts` and `tb_ensembl_affected` (separate from UniProt isoform tables).

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
tb_proteins / bp_proteins [/ tb_proteins_mus_musculus / tb_proteins_rattus_norvegicus]
  uniprot_id (PK), tim_barrel_accession (FK → entries),
  protein_name, gene_name, organism, reviewed,
  protein_existence, annotation_score,
  canonical_uniprot_id  -- NULL = canonical representative
                        -- non-NULL = redundant, points to canonical
```

### Isoforms

```
tb_isoforms / bp_isoforms [/ tb_isoforms_mus_musculus / tb_isoforms_rattus_norvegicus]
  isoform_id (PK), uniprot_id (FK → proteins), is_canonical,
  sequence, sequence_length,
  is_fragment,              -- 1 if sequence_length < 200 aa
  exon_count, exon_annotations,
  splice_variants,          -- JSON; UniProt VSP features per isoform
  tim_barrel_location,      -- JSON {domain_id, start, end, length, source}
  tim_barrel_sequence,      -- subsequence sequence[start-1:end]
  ensembl_transcript_id,    -- ENST ID from UniProt cross-reference
  alphafold_id
```

### AS-affected isoforms (UniProt)

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

### Ensembl transcript expansion (TIM barrel, Homo sapiens)

```
tb_ensembl_transcripts
  enst_id (PK), ensg_id, ensp_id,
  uniprot_id (FK → tb_proteins), gene_name,
  sequence, sequence_length, is_fragment,
  is_mane_select,          -- 1 if canonical Ensembl transcript
  biotype,
  duplicate_isoform_id,    -- isoform_id if sequence matches an existing UniProt isoform
  exon_annotations         -- JSON list of 1-based protein positions of each exon's last AA
                           --   (all coding exons except the final one; ceiling(cumulative_cds/3))

tb_ensembl_affected
  id (PK), enst_id (FK), uniprot_id (FK),
  domain_location, domain_sequence,
  canonical_domain_location, canonical_domain_sequence,
  alignment_identity, alignment_score, insertion_detected,
  exon_boundary_in_domain,          -- 1 if any exon junction falls inside [domain_start, domain_end)
  exon_boundaries_in_domain_count   -- number of exon junctions inside the domain
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

## Ensembl transcript expansion and exon boundary analysis

Expand TIM barrel coverage with Ensembl protein-coding transcripts (human only):
```bash
python scripts/collect_ensembl.py
python scripts/collect_ensembl.py --rebuild    # drop and re-collect
python scripts/collect_ensembl.py --limit 50   # test run (first 50 proteins)
```

Backfill exon boundary data and flag domain-disrupting junctions:
```bash
python scripts/backfill_exons.py
python scripts/backfill_exons.py --phase1-only  # fetch exon data only
python scripts/backfill_exons.py --phase2-only  # flag domain boundaries only
```

Results are stored in `tb_ensembl_transcripts` (with `exon_annotations`) and
`tb_ensembl_affected` (with `exon_boundary_in_domain` and `exon_boundaries_in_domain_count`).

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
    ensembl_client.py    Ensembl REST API (ENST→ENSG, transcript listing, protein sequences)
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
    schema.py     SQL DDL — tables for 2 domains x organisms + Ensembl expansion, init_db()
    connection.py get_connection() context manager, ensure_db()
    storage.py    upsert_*/get_* CRUD functions (parameterised table names)
  models/
    entities.py   Pydantic models: TIMBarrelEntry, Protein, Isoform
  query/
    engine.py     QueryEngine — SQL-backed queries (domain + organism aware)
    export.py     to_fasta(), to_csv(), to_json()
  config.py       DomainConfig + OrganismConfig registries; DOMAINS and ORGANISMS dicts
                  Note: beta_propeller is Homo sapiens only; mouse/rat for TIM barrel only
  errors.py       Exception hierarchy
  retry.py        tenacity decorator

scripts/
  collect.py                  Data collection entry point (--domain, --organism)
  build_affected_isoforms.py  AS-affected isoform detection and storage
  collect_ensembl.py          Ensembl transcript expansion for TIM barrel (Homo sapiens)
  backfill_exons.py           Fetch exon boundary data (protein-space) and flag
                              AS-affected transcripts whose exon junctions fall inside the domain
  run_hmmer.py                HMMER3 domain boundary scan using pyhmmer
```

### Domain + organism parameterisation

`DomainConfig` holds the InterPro query strategy for a domain: annotation term, text search
fallback, `cathgene3d_search` for structurally-classified CATH entries (e.g. `"3.20.20"` for
TIM barrel, `"2.130"` for beta propeller), and explicit `extra_accessions` for families with
non-standard names (e.g. WD40 superfamilies for beta propeller).
`OrganismConfig` holds the NCBI taxon ID and derives table names via
`protein_table(domain)`, `isoform_table(domain)`, and `affected_isoforms_table(domain)`.
