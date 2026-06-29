# Statistical Framework for TIM-barrel Junction Analysis

## 1. Biological background and data definitions

A **TIM barrel** (triose-phosphate isomerase barrel) is a protein fold consisting of eight
$\beta\alpha$ repeat units arranged in a closed barrel topology, conventionally numbered $k = 1, \ldots, 8$.
Each repeat unit contains a $\beta$-strand, a connecting $\beta \to \alpha$ loop, and an $\alpha$-helix;
adjacent repeat units are separated by an α→β loop. We use "inter" as the formal category
label for α→β loop positions throughout.

**Proteins.** The analysis uses UniProt (Homo sapiens) data only. A **canonical protein** is the
single best representative of each $(protein\_name, organism)$ group after deduplication: Swiss-Prot
reviewed entries are preferred, then ranked by annotation score and isoform count. Proteins in
$\mathcal{P}$ are canonical proteins with a TIM barrel domain annotation from InterPro
(Pfam or Gene3D/CATH) that have both motif annotation and exon junction data; fragment proteins
are included if they satisfy these two criteria.

**Structural annotation.** Motif boundaries ($\beta$-strand, loop, $\alpha$-helix) are derived from
the AlphaFold v4 predicted structure of each protein using DSSP secondary-structure assignment,
followed by a $\beta\alpha$ repeat detection algorithm. Proteins may have a partial annotation
($K_p < 8$ motifs) when the structure is incomplete or the algorithm cannot identify all eight repeats.

**Exon annotations.** UniProt records the last residue position of each exon in the canonical
sequence. An **exon junction** at position $j$ means that residue $j$ is the last residue of one
exon and residue $j+1$ is the first residue of the next exon. All positions are 1-based and
inclusive, in full-protein coordinates. The junction between the last exon and the 3′ UTR is
excluded (it carries no intron information within the coding sequence).

## 2. Setup and notation

Throughout, superscript $\cdot^s$ denotes the **start** (first residue) and $\cdot^e$ denotes
the **end** (last residue) of an interval; all coordinates are 1-based and inclusive in
full-protein space.

Let $\mathcal{P}$ be the set of canonical TIM barrel proteins with both motif annotation and
exon data. The enrichment analysis is restricted to the subset
$\mathcal{P}_{\text{analyzed}} = \{p \in \mathcal{P} : n_p > 0\}$ — proteins that contain
at least one eligible domain-internal exon junction.

For each protein $p \in \mathcal{P}$:

- **Domain** $D_p = [d_p^s,\, d_p^e] \subset \mathbb{Z}$, where $d_p^s$ is the first domain residue
  and $d_p^e$ is the last domain residue; domain length $L_p = d_p^e - d_p^s + 1$.
- **Motif annotation** $\mathcal{M}_p = \{(\beta_k^p,\, \lambda_k^p,\, \alpha_k^p)\}_{k=1}^{K_p}$,
  with $1 \le K_p \le 8$, where
  - $\beta_k^p = [b_k^s, b_k^e]$ is the $k$-th $\beta$-strand (start $b_k^s$, end $b_k^e$),
  - $\lambda_k^p = [b_k^e+1,\, a_k^s-1]$ is the $\beta \to \alpha$ loop,
  - $\alpha_k^p = [a_k^s, a_k^e]$ is the $k$-th $\alpha$-helix (start $a_k^s$, end $a_k^e$).
- **α→β loop linker** $\gamma_k^p = [a_k^e+1,\, b_{k+1}^s-1]$ for $k = 1,\ldots,K_p-1$.
  A linker exists only when $b_{k+1}^s > a_k^e + 1$.
- **Eligible junction positions** $E_p = \{d_p^s,\, d_p^s+1,\, \ldots,\, d_p^e - 1\}$, with
  $|E_p| = L_p - 1$. This is the set of domain positions at which an internal exon junction can
  occur; the terminal position $d_p^e$ is excluded because a junction after this residue would
  lie at or beyond the domain boundary rather than within the domain.
- **Exon junction set** $\mathcal{J}_p \subseteq E_p$, the subset of eligible positions that coincide
  with a UniProt exon boundary; $n_p = |\mathcal{J}_p|$. Positions are drawn without replacement
  in the permutation procedure because each coding boundary can contain at most one exon junction.

## 3. Element partition

For each position $x \in D_p$ define its structural element type $\tau(x, p)$:

$$\tau(x,p) = \begin{cases}
\beta & \text{if } x \in \bigcup_k \beta_k^p \\
\alpha & \text{if } x \in \bigcup_k \alpha_k^p \\
\text{inter} & \text{if } x \in \bigcup_k \gamma_k^p \\
\text{loop} & \text{if } x \in \bigcup_k \lambda_k^p \\
\text{flanking} & \text{otherwise}
\end{cases}$$

The five element types are mutually exclusive and together cover every position in $D_p$
(the partition is exhaustive by construction). Here, "flanking" refers only to unassigned positions
**within** the annotated domain $D_p$ — specifically positions before the first annotated $\beta$-strand
or after the last annotated $\alpha$-helix — and not to residues outside the domain boundaries.
Within each type all positions are pooled: motif index $k$ is not distinguished.

_Note on partial annotations:_ In proteins with $K_p < 8$, residues belonging to unannotated
motifs may be absorbed into the flanking category. Enrichment of the flanking category in partially
annotated proteins should therefore be interpreted cautiously.

Two coarsening levels are considered depending on the question:

**Simplified model** ($m = 3$) — three categories:

$$\tau_3(x,p) = \begin{cases}
\beta & \text{if } \tau(x,p) = \beta \\
\alpha & \text{if } \tau(x,p) = \alpha \\
\text{other} & \text{if } \tau(x,p) \in \{\text{inter},\, \text{loop},\, \text{flanking}\}
\end{cases}$$

$\mathcal{T}_3 = \{\beta,\, \alpha,\, \text{other}\}$.

**Full model** ($m = 5$) — five categories, $\tau_5 \equiv \tau$ as defined above.
$\mathcal{T}_5 = \{\beta,\, \alpha,\, \text{inter},\, \text{loop},\, \text{flanking}\}$.

## 4. Null model

**$H_0$ (uniform placement):** Exon junctions are placed uniformly at random over the eligible
positions $E_p$ within each protein domain.

Let $t \in \mathcal{T}_m$ denote a fixed element type under model $m$. For protein $p$,
define the per-protein probability of landing in $t$:

$$q_{pt} = \frac{\Lambda_t^p}{|E_p|} = \frac{\Lambda_t^p}{L_p - 1},
\qquad \Lambda_t^p = |\{x \in E_p : \tau_m(x,p) = t\}|$$

where $\Lambda_t^p$ counts eligible positions in $E_p$ that belong to element $t$.

The expected fraction of all junctions falling in $t$, aggregated across proteins, is:

$$\pi_t^0 = \frac{1}{N} \sum_{p} n_p\, q_{pt}
= \frac{\sum_p n_p\, \Lambda_t^p\,/\,(L_p - 1)}{N}$$

where $N = \sum_p n_p$. This weighting matches the conditional null in which the observed number
of junctions $n_p$ is fixed for each protein, making the analytical null consistent with the
permutation procedure.

_Note:_ This pooled model treats junctions as independent observations and uses a single global
expected category distribution. It therefore does not fully account for within-protein dependence,
protein-specific element composition, or variation in the number of junctions per protein. The
within-protein permutation test (§5) is the preferred inferential procedure.

**Numerical example.** For a domain of length $L = 200$ aa with eight β-strands of 5 residues
each: $\Lambda_\beta = 8 \times 5 = 40$, so $\pi_\beta^0 = 40/199 \approx 0.201$. With $N = 6$
junctions the expected number in β-strands is $E_\beta = 6 \times 0.201 \approx 1.2$.

## 5. Element-level enrichment

Let $N_t = \sum_p |\{j \in \mathcal{J}_p : \tau_m(j,p) = t\}|$ be the observed junction count in
element $t$ under model $m$, and $N = \sum_p n_p$ the total junction count.

The **observed fraction** is $f_t = N_t / N$.

The **enrichment ratio** is:

$$\rho_t = \frac{f_t}{\pi_t^0}$$

$\rho_t > 1$ denotes enrichment (junctions more frequent than expected); $\rho_t < 1$ denotes depletion.

**Analytical approximation (chi-square).** As a first-pass test, compare the observed count vector
against the pooled expected counts using a chi-square goodness-of-fit statistic. Because proteins
differ in domain length, element composition, and number of junctions, and because junction
positions within a protein are sampled without replacement, the pooled count vector is not exactly
a single multinomial sample; the chi-square test is therefore an approximation.

**Preferred: within-protein permutation test.**
Let $B$ be the number of permutation replicates (default $B = 2000$).
For each replicate $b = 1, \ldots, B$:

1. For each protein $p$, draw $n_p$ positions uniformly at random from $E_p$ without replacement.
2. Classify each drawn position using $\tau_m(\cdot, p)$ to obtain permuted counts
   $N_t^{(b)}$, observed fraction $f_t^{(b)} = N_t^{(b)} / N$, and enrichment ratio
   $\rho_t^{(b)} = f_t^{(b)} / \pi_t^0$. Note that $\pi_t^0$ is computed once from the observed
   protein set and fixed $n_p$ values, and is kept fixed across all replicates.

This procedure preserves the number of junctions $n_p$, the domain length $L_p$, and the element
structure of each protein, but does not preserve exon-length distribution or spacing between
adjacent junctions within a protein.

The two-sided permutation p-value for element $t$ is:

$$\hat{p}_t = \frac{1 + \sum_{b=1}^{B} \mathbf{1}\!\left[|\rho_t^{(b)} - 1| \ge |\rho_t - 1|\right]}{B + 1}$$

The $|\rho - 1|$ distance centres the test on the null ($\rho_t = 1$), capturing both enrichment
and depletion. The $+1$ correction avoids reporting $\hat{p}_t = 0$ for any finite $B$.

**Multiple testing correction.** Per-element permutation p-values are adjusted across categories
within each model using Benjamini–Hochberg FDR correction.

**Confidence intervals.** BH-adjusted 95% confidence intervals for $\rho_t$ are derived from the
permutation null distribution (2.5th and 97.5th percentiles of $\rho_t^{(b)}$).

## 6. Motif-specific enrichment (§8A)

For each junction $j \in \mathcal{J}_p$, assign the **motif-element label**

$$\tau_\text{motif}(j, p) = (t,\, k)$$

where $t \in \{\beta, \text{loop}, \alpha\}$ is the structural element type and $k \in \{1, \ldots, K_p\}$
is the motif number. The 31 **primary categories** are $(t, k)$ for $t \in \{\beta, \text{loop}, \alpha\}$
and $k = 1, \ldots, 8$ (24 categories), plus $(\text{inter}, k)$ for $k = 1, \ldots, 7$
(7 inter-motif linker positions). Flanking positions are tracked but excluded from the primary test.

The null expectation follows §4:

$$\pi_{(t,k)}^0 = \frac{1}{N} \sum_{p} n_p\, q_{p,(t,k)}$$

where $q_{p,(t,k)} = |E_{p,(t,k)}| / L_p$ is the fraction of domain positions in element $(t, k)$
for protein $p$, summing only over proteins that have motif $k$. The enrichment ratio and
BH-corrected chi-square test are applied across all 31 primary categories.

## 7. Within-element phase consistency (§8B)

For element types $t \in \{\alpha, \beta, \text{loop}\}$, let $[s, e]$ denote the boundaries of
an element instance. For each junction $j$ with $s \le j \le e$ and $e > s$, define the
**within-element phase**

$$\phi_j = \frac{j - s}{e - s} \in [0,1]$$

$\phi = 0$ is the first residue of the element; $\phi = 1$ is the last. Element instances of
length one ($e = s$) are excluded. Phases are pooled across all proteins and all motif instances
of type $t$.

Under $H_0$, given that a junction falls in element $[s,e]$, its position is discrete
$\text{Uniform}\{s, \ldots, e\}$, approximating $\text{Uniform}[0,1]$ after normalisation.

The KS statistic $D_t$ of $\{\phi_j\}$ against $\text{Uniform}[0,1]$ is computed. The permutation
null re-draws, for each observed junction-in-element pair $(j, [s,e])$, a position uniformly from
$\{s, \ldots, e\}$:

$$\hat{p}_{B,t} = \frac{1 + \sum_{b=1}^{B} \mathbf{1}\!\left[D_t^{(b)} \ge D_t\right]}{B+1}$$

The three permutation p-values are adjusted using Benjamini–Hochberg FDR correction.

## 8. AS junction placement (§9)

### 8.1 Dataset

Let $\mathcal{A}$ be the set of AS-affected isoforms satisfying:

1. the canonical protein $p(a) \in \mathcal{P}_{\text{analyzed}}$ (motif annotation and exon data present),
2. at least one canonical junction falls within the VSP span.

For each $a \in \mathcal{A}$, let $\mathcal{V}_a = \{[v_{ar}^s,\, v_{ar}^e]\}_r$ be the set of
**VSP canonical spans** — the intervals of the canonical sequence altered by the splice event.
The **AS-affected junction set** is

$$\mathcal{J}_a^{AS} = \bigl\{j \in \mathcal{J}_{p(a)} : \exists\, r,\ v_{ar}^s \le j \le v_{ar}^e\bigr\}$$

Total AS-affected junction-instance count: $N^{AS} = \sum_{a \in \mathcal{A}} |\mathcal{J}_a^{AS}|$.

An "AS-affected canonical junction" means a canonical exon junction lying within the canonical
span altered by the AS event — not necessarily that the alternative isoform creates a novel
boundary at that exact coordinate.

### 8.2 Element-level distribution of AS junctions (9A)

The **AS enrichment ratio relative to the canonical junction baseline** is

$$\rho_t^{AS} = \frac{f_t^{AS}}{f_t}$$

where $f_t^{AS} = N_t^{AS} / N^{AS}$ and $f_t = N_t / N$ is the canonical observed fraction from §5.

**Null** $H_0^{AS}$: for each isoform $a$, the set $\mathcal{J}_a^{AS}$ is a uniformly random
subset of $\mathcal{J}_{p(a)}$ of size $|\mathcal{J}_a^{AS}|$, drawn without replacement.

Permutation p-value:

$$\hat{p}_t^{AS} = \frac{1 + \sum_{b=1}^{B} \mathbf{1}\!\left[|\rho_t^{AS,(b)} - 1| \ge |\rho_t^{AS} - 1|\right]}{B+1}$$

Per-element p-values are adjusted using Benjamini–Hochberg FDR correction.

### 8.3 Positional distribution of AS junctions (9B)

For each junction $j \in \mathcal{J}_a^{AS}$, define the normalised domain position

$$x_{aj} = \frac{j - d_{p(a)}^s}{d_{p(a)}^e - d_{p(a)}^s} \in [0,1)$$

The test statistic

$$D_N^{AS} = \sup_{x \in [0,1]} \bigl|F_N^{AS}(x) - F_N^{const}(x)\bigr|$$

measures deviation of the AS-affected junction distribution from the canonical baseline CDF
$F_N^{const}$. The permutation null draws $|\mathcal{J}_a^{AS}|$ canonical junctions uniformly
per isoform and normalises; p-value computed as in §5.

### 8.4 Within-protein hotspot reuse (9C)

For each canonical protein $p$ with $k_p \ge 2$ AS isoforms, define the
**junction usage count** $u_p(j) = |\{a \in \mathcal{A} : p(a) = p,\, j \in \mathcal{J}_a^{AS}\}|$.

A junction is a **hotspot** if $u_p(j) \ge 2$. The aggregate hotspot fraction across all
multi-isoform proteins $\mathcal{P}^{(2)} = \{p : k_p \ge 2\}$ is

$$\bar{H} = \frac{\sum_{p \in \mathcal{P}^{(2)}} |\{j : u_p(j) \ge 2\}|}{\sum_{p \in \mathcal{P}^{(2)}} |\bigcup_a \mathcal{J}_a^{AS}|}$$

Under $H_0^{AS}$, each isoform independently draws junctions uniformly from the canonical pool.
The permutation null is generated consistently with §8.2.

## 9. VSP boundary placement (§10)

For each VSP with domain overlap, the domain-clipped start $s_v = \max(v_s,\, d_s)$ and end
$e_v = \min(v_e,\, d_e - 1)$ are assigned to structural elements using the $\tau_5$ classifier.
The enrichment ratio

$$\rho_t = \frac{f_t}{\pi_t^0}$$

compares the observed fraction of VSP boundaries in element $t$ to the length-weighted null
computed from the same canonical protein set. Significance is assessed with a chi-square z-score
test: $z_t = (O_t - E_t)/\sqrt{E_t}$, BH-corrected across five element types.

Start and end positions are tested separately (10A, 10B). Start-end asymmetry (10C) tests whether
start and end positions are drawn from the same structural distribution, using a permutation null
that randomly swaps start/end labels within each VSP.

## References

Ochoa-Leyva A, Montero-Morán G, Saab-Rincón G, Brieba LG, Soberón X. 2013. Alternative Splice
Variants in TIM Barrel Proteins from Human Genome Correlate with the Structural and Evolutionary
Modularity of this Versatile Protein Fold. *PLoS ONE* 8(8): e70582.
