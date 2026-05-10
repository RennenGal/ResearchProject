# TIM Barrel Alternative Splicing — Data Collector & Analysis

Research tool for studying the effect of alternative splicing (AS) on TIM barrel domain structure.
Collects domain family entries, proteins, and isoforms from InterPro and UniProt into a local
SQLite database; annotates the (βα)₈ repeat motifs; maps splice events onto the motif structure.

**Institution**: Ben-Gurion University of the Negev  
**Supervisors**: Prof. Tal Shay, Prof. Chen Keasar

---

## Database contents (current)

### TIM barrel — Homo sapiens

| | Count |
|---|---|
| TIM barrel families (38 Pfam/InterPro + 35 Gene3D CATH) | 73 |
| Proteins — canonical (`canonical_proteins` view) | 399 |
| &nbsp;&nbsp;of which Swiss-Prot reviewed | 193 |
| &nbsp;&nbsp;of which TrEMBL unreviewed | 206 |
| Isoforms — canonical | 399 |
| Isoforms — alternative | 249 |
| **AS-affected isoforms (domain disrupted)** | **132** |
| Genes with ≥1 AS-affected isoform | 75 |

The `proteins` table holds all 1,174 collected entries. The `canonical_proteins` view
filters to `canonical_uniprot_id IS NULL` — the best representative per
`(protein_name, organism)` group ranked by reviewed > annotation_score > isoform count.
775 proteins are marked redundant (e.g. TrEMBL fragments superseded by a Swiss-Prot entry).

AS-affected detection is VSP-based: any alternative isoform with at least one UniProt VSP feature
whose coordinates overlap `[domain_start, domain_end]` is included.

### Structural annotation (`canonical_analysis`, 810 proteins)

| Method | Tool | Coverage | Full 8-motif rate |
|---|---|---|---|
| AlphaFold DSSP | pydssp on AF2 structure | 764 / 810 | 386 / 764 (51%) |
| Per-family HMM | pyhmmer jackhmmer / phmmer | 767 / 810 | — |
| Experimental PDB | X-ray / EM ≤ 3.0 Å via PDBe | 86 / 810 | 60 / 71 AF=8 confirmed (85%) |

### Ensembl transcript expansion

| | Count |
|---|---|
| Proteins with Ensembl mapping | 799 |
| Novel unique transcripts | 1,097 |
| **AS-affected novel transcripts** | **359** |

---

## Schema

### Domain entries

```
entries
  accession (PK), entry_type  -- 'pfam' | 'interpro' | 'cathgene3d'
  name, description, domain_annotation
```

### Proteins

```
proteins
  uniprot_id (PK), tim_barrel_accession (FK → entries),
  protein_name, gene_name, organism, reviewed,
  protein_existence, annotation_score,
  canonical_uniprot_id  -- NULL = canonical representative
                        -- non-NULL = redundant, points to canonical

canonical_proteins      -- VIEW: SELECT * FROM proteins WHERE canonical_uniprot_id IS NULL
```

### Isoforms

```
isoforms
  isoform_id (PK), uniprot_id (FK), is_canonical,
  sequence, sequence_length, is_fragment,
  exon_count, exon_annotations,   -- JSON list of 1-based exon boundary positions
  splice_variants,                 -- JSON; UniProt VSP features
  tim_barrel_location,             -- JSON {domain_id, start, end, length, source}
  tim_barrel_sequence,
  ensembl_transcript_id, alphafold_id
```

### AS-affected isoforms (UniProt)

```
affected_isoforms
  isoform_id (PK), uniprot_id (FK),
  sequence, sequence_length, exon_count, exon_annotations, splice_variants,
  domain_location,                    -- JSON {start, end} in the alternative isoform
  domain_sequence,                    -- domain subsequence in the alternative isoform
  canonical_domain_location,          -- JSON {start, end} in the canonical
  canonical_domain_sequence,
  identity_percentage,
  alignment_score,
  exon_boundary_in_domain,
  exon_boundaries_in_domain_count,
  vsp_domain_events,                  -- JSON list of VSPs overlapping the domain
  detection_method,                   -- 'vsp_overlap' | 'sliding_window_fallback'
  disrupted_motifs                    -- JSON; which of the 8 motifs each VSP disrupts
```

### Ensembl transcript expansion

```
ensembl_transcripts
  enst_id (PK), ensg_id, ensp_id,
  uniprot_id (FK), gene_name,
  sequence, sequence_length, is_fragment,
  is_mane_select, biotype,
  duplicate_isoform_id, duplicate_enst_id,
  exon_annotations

ensembl_affected
  id (PK), enst_id (FK), uniprot_id (FK),
  domain_location, domain_sequence,
  canonical_domain_location, canonical_domain_sequence,
  alignment_identity, alignment_score, insertion_detected,
  exon_boundary_in_domain,
  exon_boundaries_in_domain_count,
  disrupted_motifs   -- JSON; which motifs are affected by each changed region
```

### Canonical analysis (structural annotation)

```
canonical_analysis
  uniprot_id (PK), gene_name,
  sequence, domain_start, domain_end, domain_sequence,
  exon_annotations,           -- [{exon, start, end}] 1-based inclusive
  motif_annotations,          -- [{motif, start, end, beta_start, beta_end,
                              --   alpha_start, alpha_end}] from AlphaFold DSSP
  dssp_source,
  hmmer_annotations,          -- JSON hit details from per-family HMM
  hmmer_source,
  pdb_motif_annotations,      -- [{motif, ...}] from experimental PDB
  pdb_source
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

**Full collection from scratch**:
```bash
python scripts/collect.py
```

**Resume** — collect isoforms only for proteins not yet in the database:
```bash
python scripts/collect.py --resume
```

**Re-collect all isoforms**:
```bash
python scripts/collect.py --recollect-isoforms
```

**Update entries and proteins only** (no isoform changes):
```bash
python scripts/collect.py --collect-proteins
```

**Backfill domain locations** — populate `tim_barrel_location` where it is NULL:
```bash
python scripts/collect.py --backfill-domains
```

**Options**:
```
--db PATH        Override database path (default: db/protein_data.db)
--log-file PATH  Write logs to file in addition to stdout
--log-level LVL  Logging verbosity (default: INFO)
```

---

## AS-affected isoform analysis

Identify alternative isoforms where splicing disrupts the TIM barrel domain:

```bash
python scripts/build_affected_isoforms.py
```

Results go into `affected_isoforms`. Then annotate which of the 8 motifs each splice event
disrupts:

```bash
python scripts/annotate_disrupted_motifs.py
```

---

## Ensembl transcript expansion

Expand coverage with Ensembl protein-coding transcripts:
```bash
python scripts/collect_ensembl.py
python scripts/collect_ensembl.py --rebuild    # drop and re-collect
```

Backfill exon boundary data:
```bash
python scripts/backfill_exons.py
```

---

## Structural motif annotation pipeline

Five steps annotate the (βα)₈ repeat units and cross-validate the assignments.

### Step 1 — Populate gene names
```bash
python scripts/fetch_gene_names.py
```

### Step 2 — Build canonical analysis table
```bash
python scripts/build_canonical_analysis.py
python scripts/build_canonical_analysis.py --rebuild
```

### Step 3 — Annotate β-α motifs (AlphaFold + DSSP)
```bash
python scripts/annotate_motifs.py
python scripts/annotate_motifs.py --rerun
```

Motif format:
```json
[{"motif": 1, "start": 14, "end": 37,
  "beta_start": 14, "beta_end": 20,
  "alpha_start": 27, "alpha_end": 37}, ...]
```
All positions 1-based, inclusive, in full-protein coordinates.

### Step 4 — Cross-validate with per-family HMMs
```bash
python scripts/cross_validate_hmmer.py
```

### Step 5 — Validate with experimental PDB structures
```bash
python scripts/validate_pdb_experimental.py
python scripts/validate_pdb_experimental.py --resolution 2.5
```

---

## Exon / motif analysis

```bash
python scripts/analyze_exon_junctions.py
python scripts/analyze_exon_junctions.py --full-only   # restrict to 8-motif proteins
```

Answers two questions:
- **Q1**: How are constitutive exon junctions distributed across the β-α motif elements (gene architecture)?
- **Q2**: Which junctions are exploited by AS events (VSP spans cross-referenced against canonical junctions)?

Key finding: α-helices are the most enriched element in gene structure (1.49×), robust across all
domain-length bins. AS events show near-neutral distribution relative to that baseline.
See `Statistical-Analysis.md` for the full formal results.

### Domain-length subgroup analysis

```bash
python scripts/analyze_domain_length_subgroups.py
python scripts/analyze_domain_length_subgroups.py --min 200 --max 400 --step 50
```

Splits proteins into 50 aa domain-length bins and computes enrichment per structural element in
each bin, testing whether the Q1 signal is an artefact of domain-length heterogeneity.

### Junction alignment plots

```bash
# Four panels: 200-400 aa in 50 aa bins
python scripts/plot_junction_alignment.py

# Single panel: all proteins
python scripts/plot_junction_alignment.py --min 1 --max 9999 --no-subgroups \
    --out figures/junction_alignment_all.png
```

One row per protein (sorted by domain length), coloured dots mark exon junction positions
normalised to [0, 1] within the domain. Background stripes show structural element positions
from a representative protein. Output written to `figures/`.

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
pytest tests/ -v
pytest tests/ -v --tb=short
```

Test files:
- `tests/test_models.py` — Pydantic model validation (sequence, length, bounds, fragment flag, domain sequence slicing)
- `tests/test_database.py` — CRUD round-trips, FK enforcement, new fields, resume detection
- `tests/test_uniprot_parser.py` — Pure parser functions against a real UniProt response
- `tests/test_query.py` — QueryEngine queries, FASTA/CSV/JSON export format
- `tests/test_motif_annotator.py` — `identify_ba_motifs()` (16 tests: clean barrel, partial barrel, MERGE_GAP fusion, domain boundary edge cases, min-length filtering) and `_reformat_exon_annotations()` (8 tests)

---

## Architecture

```
protein_data_collector/
  api/
    interpro_client.py   InterPro REST API (pagination, multi-strategy entry search)
    uniprot_client.py    UniProt REST API (protein JSON, isoform FASTA)
    ensembl_client.py    Ensembl REST API (transcript listing, protein sequences)
  collector/
    interpro_collector.py   Phase 1+2: domain families and proteins
    uniprot_collector.py    Phase 3: isoforms and splice variant extraction
    data_collector.py       Pipeline orchestrator (CollectionReport, backfill, deduplication)
  database/
    schema.py     SQL DDL — TIM barrel / Homo sapiens tables, init_db()
    connection.py get_connection() context manager
    storage.py    upsert_*/get_* CRUD functions
  models/
    entities.py   Pydantic models: TIMBarrelEntry, Protein, Isoform
  query/
    engine.py     QueryEngine — SQL-backed queries
    export.py     to_fasta(), to_csv(), to_json()
  analysis/
    motif_annotator.py         identify_ba_motifs(): β-α repeat detection from SS array
    tim_barrel_alignment.py    populate_tim_barrel_isoforms(): VSP-based AS detection
  config.py       DomainConfig + OrganismConfig; DOMAINS and ORGANISMS registries
  errors.py       Exception hierarchy
  retry.py        tenacity decorator

scripts/
  collect.py                    Data collection entry point
  build_affected_isoforms.py    Build affected_isoforms (VSP-based detection)
  annotate_disrupted_motifs.py  Map each splice event to the motifs it disrupts
  collect_ensembl.py            Ensembl transcript expansion
  backfill_isoform_exons.py     Backfill exon junction data for UniProt isoforms
  backfill_exons.py             Backfill exon junction data for Ensembl transcripts
  run_hmmer.py                  HMMER3 domain boundary scan using pyhmmer
  fetch_gene_names.py           Batch-fetch gene names from UniProt
  build_canonical_analysis.py   Build canonical_analysis
  annotate_motifs.py            AlphaFold + pydssp motif annotation
  cross_validate_hmmer.py       Per-family HMM cross-validation
  validate_pdb_experimental.py  Experimental PDB validation via PDBe + RCSB
  analyze_exon_junctions.py          Q1/Q2 exon–motif junction analysis
  analyze_domain_length_subgroups.py Junction enrichment by domain-length bin
  plot_junction_alignment.py         Junction alignment dot plot (per-protein rows)

tests/
  test_models.py
  test_database.py
  test_uniprot_parser.py
  test_query.py
  test_motif_annotator.py
```
