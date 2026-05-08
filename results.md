# TIM Barrel Alternative Splicing — Results

Data source: UniProt (Homo sapiens, Swiss-Prot reviewed + TrEMBL).
Structural annotations: AlphaFold v4 + DSSP.
Analysis includes all proteins with motif annotations (full and partial) unless noted.

---

## 1. Dataset

| | Count |
|---|---|
| Proteins collected (total) | 1,174 |
| Canonical proteins after deduplication | 399 |
| — Swiss-Prot reviewed | 193 |
| — TrEMBL unreviewed | 206 |
| Canonical proteins in `canonical_analysis` (non-fragment) | 311 |
| — with full 8-motif annotation | 188 |
| — with partial motif annotation (1–7) | 110 |
| — no AlphaFold structure | 8 |
| — fragment (excluded from structural analysis) | 88 |
| Distinct gene names (Swiss-Prot reviewed) | 199 |
| Proteins without gene name (TrEMBL) | 47 |
| Alternative isoforms (from canonical proteins) | 249 |
| **AS-affected isoforms (domain-disrupting)** | **132** |
| Genes with at least one AS-affecting isoform | 75 |

Deduplication groups proteins by `(protein_name, organism)`; the best
representative (reviewed > annotation_score > isoform count) is kept canonical.
775 proteins were marked redundant (e.g. TrEMBL fragments with a Swiss-Prot counterpart).

AS detection method: VSP overlap. Any UniProt VSP feature whose canonical
coordinates overlap [domain_start, domain_end] by at least 1 residue is
classified as domain-affecting.

### Distribution of AS variants per gene

| Variants per gene | Genes |
|---|---|
| 1 | 42 |
| 2 | 21 |
| 3 | 5 |
| 4 | 5 |
| 6 | 1 |
| 7 | 1 |

Top genes by variant count: MOCS1 (7), SLC3A1 (6), TATDN3 / KCNAB2 / HYAL1 (4 each).

---

## 2. Domain length distribution

Before performing junction analysis, the distribution of TIM barrel domain
lengths was examined across all 311 canonical non-fragment proteins.

| Statistic | Value |
|---|---|
| Min | 49 aa |
| Max | 763 aa |
| Mean | 290 aa |
| Median | 302 aa |
| Std dev | 116 aa |
| CV | 40% |

**Length distribution:**

| Range (aa) | Count | % |
|---|---|---|
| < 200 | 62 | 19.9% |
| 200–250 | 33 | 10.6% |
| 250–300 | 54 | 17.4% |
| 300–350 | 70 | 22.5% |
| 350–400 | 54 | 17.4% |
| 400–500 | 29 | 9.3% |
| > 500 | 9 | 2.9% |

The bulk (68%) fall in the 200–400 aa range, consistent with expected TIM barrel
domain size. Two notable outlier classes exist:

- **Very short (<200 aa, 20%):** Entries such as GLA (49 aa, residues 32–80)
  likely reflect mis-annotated or clipped domain boundaries from InterPro;
  a TIM barrel cannot fold in fewer than ~200 residues.
- **Very long (>500 aa, 3%):** PLCG1/PLCG2 (742–763 aa) are multi-domain
  phospholipase C proteins where the TIM barrel catalytic core spans most of
  the protein length.

The 40% CV confirms that domain lengths are not uniform, making length
normalization essential for interpreting junction placement statistics.

---

## 3. Q1 — Exon structure of the TIM barrel domain

**Question:** Where do constitutive exon-intron boundaries fall within the
(βα)₈ motif structure? This reflects gene architecture, independent of AS.

Dataset: 231 proteins with motif annotation and exon data (188 full 8-motif +
43 partial). Total exon junctions inside domain: 1,455. Per-protein mean 6.3,
median 7.

### Broad classification

| Category | Count | Raw % | Total residues | Junctions/residue | Enrichment |
|---|---|---|---|---|---|
| alpha helices | 459 | 31.5% | 19,694 | 0.0233 | **1.49x** |
| inter-motif gaps | 200 | 13.7% | 9,653 | 0.0207 | **1.32x** |
| flanking (outside domain) | 285 | 19.6% | 15,055 | 0.0189 | 1.21x |
| beta strands | 145 | 10.0% | 9,316 | 0.0156 | 0.99x |
| loops (beta→alpha) | 366 | 25.2% | 39,277 | 0.0093 | **0.60x** |

Enrichment = (junctions/residue) / (total junctions / total domain residues).

**Key finding:** Raw counts overstate the contribution of loops because
β→α loops are very long (~39k total residues vs ~20k for helices). When
normalized by length, alpha helices are the most enriched element (1.49×),
inter-motif gaps are modestly enriched (1.32×), and loops are strongly
depleted (0.60×). The gene structure of human TIM barrel proteins
preferentially places exon boundaries inside helices.

### Inter-motif gap detail

| Gap | Count | Raw % | Residues | Enrichment |
|---|---|---|---|---|
| motif 1 > 2 | 39 | 2.7% | 2,494 | 1.00x |
| motif 2 > 3 | 28 | 1.9% | 1,229 | 1.46x |
| motif 3 > 4 | 35 | 2.4% | 1,530 | 1.46x |
| motif 4 > 5 | 43 | 3.0% | 1,798 | 1.53x |
| motif 5 > 6 | 24 | 1.6% | 945 | **1.62x** |
| motif 6 > 7 | 12 | 0.8% | 841 | 0.91x |
| motif 7 > 8 | 19 | 1.3% | 816 | 1.49x |

Gaps between motifs 5–6 and 4–5 are the most enriched inter-motif positions.

---

## 4. Domain-length subgroup analysis (Q1)

To test whether the Q1 enrichment pattern is an artefact of domain length
heterogeneity, the 200–400 aa proteins (165 proteins, 1,083 junctions) were
split into four 50 aa bins and the enrichment computed independently for each.

### Enrichment by subgroup

| Category | 200–250 (22p) | 250–300 (44p) | 300–350 (60p) | 350–400 (39p) |
|---|---|---|---|---|
| alpha helices | 1.21x | 1.38x | **1.61x** | **1.64x** |
| inter-motif gaps | 0.95x | 1.31x | 1.53x | 1.28x |
| flanking | 1.02x | 0.94x | 1.13x | 1.24x |
| beta strands | **1.68x** | 0.89x | 0.64x | 0.98x |
| loops (beta→alpha) | 0.43x | 0.65x | 0.69x | 0.47x |

### Key observations

- **α-helix enrichment is the most stable signal** — present in all four
  bins (1.21×–1.64×), confirming it is not a domain-length artefact.
- **Loops are depleted in all bins** (0.43×–0.69×), also robust to length.
- **β-strand enrichment is length-dependent**: strongly enriched at 200–250 aa
  (1.68×) but drops to neutral at longer lengths. Shorter annotated domains
  may clip the N- or C-terminal strands, inflating the apparent strand count
  relative to the compressed region.
- **Flanking enrichment rises with length** (1.02× → 1.24× at 350–400 aa),
  consistent with longer domains leaving more sequence outside the annotated
  barrel core.
- **Between-motif gaps** peak at 300–350 aa (1.53×) and are less enriched at
  the extremes, suggesting the gap regions are best resolved at intermediate
  domain lengths.

**Conclusion:** The α-helix enrichment (Q1 key finding) is robust across all
domain length bins and is not explained by the 40% variation in domain length
across the dataset.

---

## 5. Q2 — Where do AS events hit the domain?

**Question:** Among the exon junctions that exist in the canonical gene, which
ones are actually exploited by AS events?

Method: For each AS-affected isoform, the VSP canonical span [can_start,
can_end] defines the region changed by the splice event. Canonical exon
junctions that fall within this span are the actual splice boundaries exploited.

Dataset: 131 AS-affected isoforms with exon and motif data.
94 isoforms have at least 1 canonical junction inside their VSP span.
Total AS-exploited junction instances: 218.

### Q2 vs. Q1 comparison

| Category | Q2 count | Q2 % | Q1 % | AS enrichment |
|---|---|---|---|---|
| beta strands | 31 | 14.2% | 10.0% | **1.43x** |
| loops (beta→alpha) | 57 | 26.1% | 25.2% | ~1.04x (neutral) |
| inter-motif gaps | 24 | 11.0% | 13.7% | 0.80x |
| alpha helices | 54 | 24.8% | 31.5% | **0.79x** (avoided) |
| flanking | 52 | 23.9% | 19.6% | ~1.22x |

AS enrichment = Q2% / Q1%. A value >1 means that element type is
over-represented in AS events relative to its frequency in the gene structure.

**Key finding:** When controlling for the underlying exon structure (Q1), AS
events show no preference for inter-motif gaps or loops. The distribution is
close to neutral across all categories. The mild exceptions are a slight
enrichment of beta-strand junctions (1.43×) and a slight avoidance of
alpha-helix junctions (0.79×).

The largest single Q2 category is `after_barrel` (41/218 = 18.8%): many AS
events span from inside the domain to beyond the C-terminal domain boundary,
removing the tail of the barrel.

### Top Q2 elements

| Element | Count | % |
|---|---|---|
| after_barrel | 41 | 18.8% |
| loop_3 | 11 | 5.0% |
| loop_1 | 11 | 5.0% |
| before_barrel | 11 | 5.0% |
| loop_2 | 10 | 4.6% |
| loop_5 | 9 | 4.1% |
| alpha_7 | 9 | 4.1% |
| inter_2_3 | 8 | 3.7% |

---

## 6. VSP boundary placement

A separate view: where do the VSP start/end coordinates themselves fall in the
motif structure? These are protein-level change boundaries, not exon-level
splice sites.

Dataset: 132 isoforms, 171 VSP domain events.

| Category | Start % | End % | Combined % |
|---|---|---|---|
| flanking (before/after domain) | 51.5% | 45.0% | **48.2%** |
| alpha helices | 24.0% | 18.7% | 21.3% |
| loops (beta→alpha) | 13.5% | 21.6% | 17.5% |
| beta strands | 8.2% | 9.9% | 9.1% |
| inter-motif gaps | 2.9% | 4.7% | 3.8% |

~48% of VSP boundaries fall outside the domain. This is expected: a VSP that
spans e.g. residues 56–330 can affect the domain even if the splice site itself
is upstream/downstream of the domain boundary.

Top start-boundary elements: `before_barrel` (31.0%), `after_barrel` (20.5%),
`alpha_4` (5.3%), `alpha_7` (5.3%).

---

## 7. Relation to Ochoa-Leyva et al. 2013

The paper reported that most AS events in TIM barrel proteins fall between
the βα motifs (inter-motif loops and linker regions). Our data does not
confirm this for the full set of human TIM barrel proteins.

The discrepancy likely has two sources:

1. **Raw vs. normalized counts.** Without normalizing by element length,
   loops dominate (25.2% raw) because they are long. The paper may not
   have applied length normalization.

2. **Gene structure vs. AS preference (Q1 vs. Q2).** Even after normalization,
   our Q2 analysis shows that AS events do not specifically target inter-motif
   regions beyond what the baseline gene structure predicts (0.80×). The
   paper's dataset (specific enzyme families) may have a different baseline
   gene structure.

The Q1 finding is the stronger result: human TIM barrel genes are structured
such that alpha helices are the most exon-boundary-enriched element (1.49×,
stable across all domain length bins). Whether AS preferentially targets any
particular element relative to this baseline is not strongly supported by our
data.

---

## 8. Open questions

- 37 AS-affected isoforms have VSP spans containing no canonical exon junction.
  These are likely large deletions spanning whole exons, or isoforms where
  exon annotation is incomplete.
- The Q2 enrichment of beta-strand junctions (1.43×) warrants closer
  inspection: are these concentrated in specific motif positions or proteins?
- The `after_barrel` dominance in Q2 (18.8%) suggests a common mode of
  C-terminal domain truncation via AS; these isoforms may produce truncated
  or structurally distinct domain variants.
- The very short domain annotations (<200 aa, 20% of proteins) likely reflect
  InterPro boundary clipping. Expanding domain boundaries by 10–15 residues
  and re-running motif detection may recover a full 8th motif for some of
  these proteins.
