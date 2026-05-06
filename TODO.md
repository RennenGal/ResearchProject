# Project TODO

## Current state (2025-05-06)

### Data collected
| Table | Count |
|---|---|
| TIM barrel families (Pfam/InterPro + CATH Gene3D) | 73 |
| Canonical human proteins (`tb_proteins`) | 1,174 |
| Canonical isoforms (`tb_isoforms`) | 1,174 |
| Alternative isoforms | 249 |
| AS-affected isoforms (`tb_affected_isoforms`) | 37 |
| Ensembl novel transcripts (`tb_ensembl_transcripts`) | 1,097 unique |
| AS-affected Ensembl transcripts (`tb_ensembl_affected`) | 359 |

### Structural annotation (`tb_canonical_analysis`, 810 proteins)
| Step | Method | Result |
|---|---|---|
| Motif annotation | AlphaFold DSSP | 386 with full 8 motifs, 307 partial, 46 no structure |
| Cross-validation | Per-family HMM (jackhmmer / phmmer) | 85% of 7-motif proteins confirmed as TIM barrels |
| Independent validation | Experimental X-ray/EM structures (≤3.0 Å) | 86 proteins; 60/71 AF=8 confirmed (85%) |

---

## TODO

### High priority

- [ ] **Annotate motif locations for AS-affected isoforms**
  For each of the 37 UniProt AS-affected isoforms and 359 Ensembl AS-affected transcripts,
  map each AS event (exon skipping, insertion, substitution) onto the canonical motif
  coordinates to determine which of the 8 β-α repeat units is disrupted. Store as
  `disrupted_motifs` JSON in `tb_affected_isoforms` / `tb_ensembl_affected`.

- [ ] **Collect missing Gene3D TIM barrel proteins**
  ~292 proteins classified under CATH Gene3D TIM barrel entries were not collected
  in the initial run. Identify the missing entries and re-run collection.

- [ ] **Restore missing WD40 / beta-propeller proteins**
  ~1,500 beta-propeller proteins were lost when `bp_entries` was cleaned up.
  Restore by re-running `collect.py --domain beta_propeller`.

### Analysis

- [ ] **Motif-level disruption statistics**
  Count how often each of the 8 motif positions is disrupted across all AS events.
  Determine whether disruptions are uniform or biased toward specific barrel positions
  (e.g. motifs 1–2 vs. 7–8).

- [ ] **Exon–motif boundary overlap**
  Cross-reference `exon_annotations` positions with `motif_annotations` positions
  to identify exon boundaries that fall precisely at β→loop or loop→α transitions.

- [ ] **Extend structural annotation to mouse and rat**
  Run `build_canonical_analysis.py`, `annotate_motifs.py`, and validation scripts
  for `tb_proteins_mus_musculus` and `tb_proteins_rattus_norvegicus`.

- [ ] **Write motif-disruption mapping script**
  Join `tb_canonical_analysis.motif_annotations` with `tb_affected_isoforms.domain_location`
  and the VSP splice variant coordinates to determine which of the 8 β-α units each AS event
  overlaps. `domain_location` already exists in `tb_affected_isoforms`; the missing piece is
  the overlap logic and `disrupted_motifs` storage. Same for `tb_ensembl_affected`.

- [ ] **Investigate partial barrels (7-motif proteins)**
  116 proteins show 7 DSSP motifs; HMMs confirm ~85% are genuine TIM barrels. The missing
  motif is likely a domain boundary artefact (InterPro clips the N- or C-terminal strand).
  Try expanding `domain_start`/`domain_end` by up to 15 residues and re-running
  `identify_ba_motifs()` to see how many recover a full 8th motif.

- [ ] **Export annotated dataset to CSV**
  Write an export script (or QueryEngine method) that joins `tb_canonical_analysis`,
  `tb_proteins`, `tb_affected_isoforms`, and `tb_ensembl_affected` into a tidy flat CSV
  suitable for R / pandas analysis and sharing with supervisors.

### Infrastructure

- [ ] **Fix schema.py to include all runtime-added columns**
  `hmmer_annotations`, `hmmer_source`, `pdb_motif_annotations`, and `pdb_source` are
  added via `ALTER TABLE` inside script `ensure_columns()` calls but are absent from the
  `CREATE TABLE` DDL in `schema.py`. A fresh `init_db()` produces a schema that diverges
  from the live DB, breaking reruns from scratch.

- [ ] **Add sanity-check for PDB coordinate mapping**
  The `validate_pdb_experimental.py` offset logic (`min_uni + shift`) was verified on
  P00813 but not systematically. Write a check that, for all 86 annotated proteins,
  confirms that the first and last motif positions fall within `[domain_start, domain_end]`
  and flag any that don't (possible insertion-code or non-standard chain numbering issues).

- [ ] **Add tests for motif annotator**
  Unit tests for `identify_ba_motifs()` covering: clean 8-motif barrel, partial
  barrel, merged strands/helices (MERGE_GAP), edge cases (domain at protein ends).

- [ ] **Add tests for `build_canonical_analysis._reformat_exon_annotations()`**
  Cover: single exon, multi-exon, unsorted boundaries, seq_len at boundary.
