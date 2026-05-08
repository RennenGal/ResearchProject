# TIM Barrel Alternative Splicing — Results

Data source: UniProt (Homo sapiens, Swiss-Prot reviewed + TrEMBL).
Structural annotations: AlphaFold v4 + DSSP.
Analysis restricted to 8-motif proteins unless noted.

---

## 1. Dataset

| | Count |
|---|---|
| Canonical proteins with TIM barrel annotation | 810 |
| Distinct gene names (Swiss-Prot reviewed) | 222 |
| Proteins without gene name (TrEMBL) | 183 |
| Alternative isoforms (non-fragment) | 224 |
| **AS-affected isoforms (domain-disrupting)** | **132** |
| Genes with at least one AS-affecting isoform | 75 |

AS detection method: VSP overlap. Any UniProt VSP feature whose canonical
coordinates overlap [domain_start, domain_end] by at least 1 residue is
classified as domain-affecting. Sliding-window identity (12.5%–95%) is used
only as a fallback for the 6 isoforms with no VSP annotations.

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

## 2. Q1 — Exon structure of the TIM barrel domain

**Question:** Where do constitutive exon-intron boundaries fall within the
(ba)8 motif structure? This reflects gene architecture, independent of AS.

Dataset: 241 proteins with full 8-motif annotation and exon data.
Total exon junctions inside domain: 1842. Per-protein mean 7.6, median 7.

### Broad classification

| Category | Count | Raw % | Total residues | Junctions/residue | Enrichment |
|---|---|---|---|---|---|
| alpha helices | 586 | 31.8% | 23,768 | 0.0247 | **1.56x** |
| inter-motif gaps | 233 | 12.6% | 11,132 | 0.0209 | **1.32x** |
| flanking (outside domain) | 393 | 21.3% | 19,308 | 0.0204 | 1.29x |
| beta strands | 183 | 9.9% | 11,205 | 0.0163 | 1.03x |
| loops (beta->alpha) | 447 | 24.3% | 51,124 | 0.0087 | **0.55x** |

Enrichment = (junctions/residue) / (total junctions / total domain residues).

**Key finding:** Raw counts overstate the contribution of loops because
beta->alpha loops are very long (~51k total residues across 241 proteins vs
~24k for helices). When normalized by length, alpha helices are the most
enriched element (1.56x), inter-motif gaps are modestly enriched (1.32x),
and loops are strongly depleted (0.55x). The gene structure of human TIM
barrel proteins preferentially places exon boundaries inside helices, not
between motifs.

### Inter-motif gap detail

| Gap | Count | Raw % | Residues | Enrichment |
|---|---|---|---|---|
| motif 1 > 2 | 49 | 2.7% | 2,992 | 1.04x |
| motif 2 > 3 | 27 | 1.5% | 1,255 | 1.36x |
| motif 3 > 4 | 38 | 2.1% | 1,755 | 1.37x |
| motif 4 > 5 | 35 | 1.9% | 1,559 | 1.42x |
| motif 5 > 6 | 33 | 1.8% | 1,154 | **1.81x** |
| motif 6 > 7 | 15 | 0.8% | 1,076 | 0.88x |
| motif 7 > 8 | 36 | 2.0% | 1,341 | **1.70x** |

The gap between motifs 5 and 6 and between 7 and 8 are the most enriched
inter-motif positions in the gene structure.

---

## 3. Q2 — Where do AS events hit the domain?

**Question:** Among the exon junctions that exist in the canonical gene, which
ones are actually exploited by AS events?

Method: For each AS-affected isoform, the VSP canonical span [can_start,
can_end] defines the region changed by the splice event. Canonical exon
junctions that fall within this span are the actual splice boundaries exploited.

Dataset: 108 AS-affected isoforms with 8-motif annotation.
78 isoforms have at least 1 canonical junction inside their VSP span.
Total AS-exploited junction instances: 193.

### Q2 vs. Q1 comparison

| Category | Q2 count | Q2 % | Q1 % | AS enrichment |
|---|---|---|---|---|
| beta strands | 27 | 14.0% | 9.9% | **1.41x** |
| loops (beta->alpha) | 48 | 24.9% | 24.3% | ~1.0x (neutral) |
| inter-motif gaps | 23 | 11.9% | 12.6% | ~0.94x (neutral) |
| alpha helices | 49 | 25.4% | 31.8% | **0.80x** (avoided) |
| flanking | 46 | 23.8% | 21.3% | ~1.1x |

AS enrichment = Q2% / Q1%. A value >1 means that element type is
over-represented in AS events relative to its frequency in the gene structure.

**Key finding:** When controlling for the underlying exon structure (Q1), AS
events show no preference for inter-motif gaps or loops. The distribution is
close to neutral across all categories. The mild exceptions are a slight
enrichment of beta-strand junctions (1.41x) and a slight avoidance of
alpha-helix junctions (0.80x).

The largest single Q2 category is `after_barrel` (40/193 = 20.7%): many AS
events span from inside the domain to beyond the C-terminal domain boundary,
removing the tail of the barrel.

### Top Q2 elements

| Element | Count | % |
|---|---|---|
| after_barrel | 40 | 20.7% |
| loop_3 | 10 | 5.2% |
| loop_5 | 9 | 4.7% |
| loop_1 | 9 | 4.7% |
| inter_2_3 | 8 | 4.1% |
| alpha_7 | 8 | 4.1% |
| alpha_5 | 7 | 3.6% |
| alpha_4 | 7 | 3.6% |

---

## 4. VSP protein boundary placement

A separate (cruder) view: where do the VSP start/end coordinates themselves
fall in the motif structure? These are protein-level change boundaries, not
exon-level splice sites.

Dataset: 109 isoforms, 144 VSP domain events.

| Category | Start % | End % | Combined % |
|---|---|---|---|
| flanking (before/after domain) | 54.9% | 46.5% | **50.7%** |
| alpha helices | 21.5% | 18.8% | 20.1% |
| loops (beta->alpha) | 13.2% | 20.1% | 16.7% |
| beta strands | 8.3% | 11.1% | 9.7% |
| inter-motif gaps | 2.1% | 3.5% | 2.8% |

50.7% of VSP boundaries fall outside the domain. This is expected: a VSP
that spans e.g. residues 56–330 can affect the domain even if the splice
site itself is upstream/downstream of the domain. This analysis classifies
the protein-level change extent, not the splice site position.

---

## 5. Relation to Ochoa-Leyva et al. 2013

The paper reported that most AS events in TIM barrel proteins fall between
the ba motifs (inter-motif loops and linker regions). Our data does not
confirm this for the full set of human TIM barrel proteins.

The discrepancy likely has two sources:

1. **Raw vs. normalized counts.** Without normalizing by element length,
   loops dominate (24.3% raw) because they are long. The paper may not
   have applied length normalization.

2. **Gene structure vs. AS preference (Q1 vs. Q2).** Even after normalization,
   our Q2 analysis shows that AS events do not specifically target inter-motif
   regions beyond what the baseline gene structure predicts (0.94x). The
   paper's dataset (specific enzyme families) may have a different baseline
   gene structure.

The Q1 finding is the stronger result: human TIM barrel genes are structured
such that alpha helices are the most exon-boundary-enriched element. Whether
AS preferentially targets any particular element relative to this baseline
is not strongly supported by our data.

---

## 6. Open questions

- 30 AS-affected isoforms (of 108 in Q2) have VSP spans containing no
  canonical exon junction. These are likely large deletions or substitutions
  that span whole exons without a junction inside the affected region, or
  isoforms where the exon annotation is incomplete.
- The Q2 enrichment of beta-strand junctions (1.41x) warrants closer
  inspection: are these concentrated in specific motif positions or proteins?
- The `after_barrel` dominance in Q2 (20.7%) suggests a common mode of
  C-terminal domain truncation via AS; these isoforms may produce truncated
  or structurally distinct domain variants.
