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

### Infrastructure

- [ ] **Add tests for motif annotator**
  Unit tests for `identify_ba_motifs()` covering: clean 8-motif barrel, partial
  barrel, merged strands/helices (MERGE_GAP), edge cases (domain at protein ends).

- [ ] **Add tests for `build_canonical_analysis._reformat_exon_annotations()`**
  Cover: single exon, multi-exon, unsorted boundaries, seq_len at boundary.
