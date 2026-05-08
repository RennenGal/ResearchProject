# Project TODO

## Current state (2026-05-08)

### Data collected
| Table | Count |
|---|---|
| TIM barrel families (Pfam/InterPro + CATH Gene3D) | 73 |
| Proteins total collected | 1,174 |
| Canonical proteins (deduplicated) | 399 (193 reviewed + 206 TrEMBL) |
| Canonical isoforms (`isoforms`) | 399 |
| Alternative isoforms (canonical proteins) | 249 |
| AS-affected isoforms (`affected_isoforms`) | 132 (VSP-based detection) |
| Ensembl novel transcripts (`ensembl_transcripts`) | 1,097 unique |
| AS-affected Ensembl transcripts (`ensembl_affected`) | 359 |

### Structural annotation (`canonical_analysis`, 810 proteins)
| Step | Method | Result |
|---|---|---|
| Motif annotation | AlphaFold DSSP | 386 with full 8 motifs, 307 partial, 46 no structure |
| Cross-validation | Per-family HMM (jackhmmer / phmmer) | 85% of 7-motif proteins confirmed as TIM barrels |
| Independent validation | Experimental X-ray/EM structures (≤3.0 Å) | 86 proteins; 60/71 AF=8 confirmed (85%) |

---

## TODO

### High priority

- [x] **Annotate motif locations for AS-affected isoforms**
  Done in `scripts/annotate_disrupted_motifs.py`. UniProt: VSP overlap_start/end mapped to
  motifs (132 annotated). Ensembl: position-by-position sequence comparison of canonical vs
  alt domain sequence, changed regions mapped to motifs (330 annotated, 11 skipped — no
  motif data). `disrupted_motifs` column added to both tables.

- [x] **Collect missing Gene3D TIM barrel proteins**
  Investigated: all proteins reachable via Gene3D entries are already in proteins
  under their Pfam/InterPro accessions (verified for the two largest gaps: G3DSA:3.20.20.70
  has 293 proteins in the API — all 293 already in DB; G3DSA:3.20.20.190 has 118 — all 118
  in DB). The collection pipeline deduplicates by uniprot_id so they were collected
  via the richer Pfam/InterPro entries. No action needed.

### Low priority

- [ ] **Restore missing WD40 / beta-propeller proteins**
  ~1,500 beta-propeller proteins were lost when `bp_entries` was cleaned up.
  Restore by re-running `collect.py --domain beta_propeller`.

### Analysis

- [ ] **Motif-level disruption statistics**
  Count how often each of the 8 motif positions is disrupted across all AS events.
  Determine whether disruptions are uniform or biased toward specific barrel positions
  (e.g. motifs 1–2 vs. 7–8).

- [x] **Exon–motif boundary overlap**
  Done in `scripts/analyze_exon_junctions.py`. Q1 (all domain junctions) and Q2
  (AS-exploited junctions) both classify every junction against the 8-motif structure.
  Key finding: alpha helices are the most enriched element in gene structure (1.56x);
  AS events show near-neutral distribution relative to that baseline (see `results.md`).

- [x] **Write motif-disruption mapping script**
  Done — see `scripts/annotate_disrupted_motifs.py` and the high-priority item above.

- [ ] **Investigate partial barrels (7-motif proteins)**
  116 proteins show 7 DSSP motifs; HMMs confirm ~85% are genuine TIM barrels. The missing
  motif is likely a domain boundary artefact (InterPro clips the N- or C-terminal strand).
  Try expanding `domain_start`/`domain_end` by up to 15 residues and re-running
  `identify_ba_motifs()` to see how many recover a full 8th motif.

- [ ] **Export annotated dataset to CSV**
  Write an export script (or QueryEngine method) that joins `canonical_analysis`,
  `proteins`, `affected_isoforms`, and `ensembl_affected` into a tidy flat CSV
  suitable for R / pandas analysis and sharing with supervisors.

### Infrastructure

- [x] **Backfill protein metadata and fix deduplication**
  Root cause: `interpro_collector.py` never fetched `protein_name`, `reviewed`, or
  `annotation_score` from UniProt, so `deduplicate_proteins()` (which groups by
  `protein_name`) had never functioned. Fixed by:
  - Adding `batch_protein_metadata()` to `UniProtClient` (fields: protein_name, reviewed,
    annotation_score; handles both `recommendedName` and `submissionNames` for TrEMBL).
  - Writing `scripts/backfill_protein_metadata.py` which fetches metadata for all proteins
    with NULL fields and then re-runs deduplication.
  Result: 1,174 total → 399 canonical (193 Swiss-Prot reviewed + 206 TrEMBL); 775 redundant.
  B4DHM2 (the triggering case) correctly marked redundant to Q9NZK5.

- [x] **Fix schema.py to include all runtime-added columns**
  Added `hmmer_annotations`, `pdb_motif_annotations`, `pdb_source` to `canonical_analysis`;
  `vsp_domain_events`, `detection_method`, `disrupted_motifs` to `affected_isoforms`;
  `disrupted_motifs` to `ensembl_affected`. DDL now matches live DB.

- [x] **Add sanity-check for PDB coordinate mapping**
  `check_pdb_coordinates()` added to `scripts/validate_pdb_experimental.py`. Runs
  automatically at the end of the script; flags any protein whose motif span falls
  outside `[domain_start, domain_end]`. (PDB annotations currently absent from DB —
  re-run `validate_pdb_experimental.py` to regenerate them.)

- [x] **Add tests for motif annotator**
  16 tests in `tests/test_motif_annotator.py`: clean 8-motif barrel, partial barrel,
  MERGE_GAP fusion, gap-too-large separation, domain offset, protein-boundary edge
  cases, min-length filtering for strands and helices.

- [x] **Add tests for `build_canonical_analysis._reformat_exon_annotations()`**
  8 tests in same file: single exon, multi-exon, unsorted boundaries, boundary at
  seq_len, sequential numbering, first starts at 1, last ends at seq_len, contiguous.
