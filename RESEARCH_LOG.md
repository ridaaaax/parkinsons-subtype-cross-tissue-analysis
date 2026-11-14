# Research Log and Analytical Decision Record

## Project

**Evaluating Peripheral Transcriptomic Signatures Associated with Brain-Defined Parkinson's Disease Subtypes**

**Project type:** Self-directed exploratory bioinformatics investigation  
**Research approach:** Public transcriptomic data analysis, differential gene expression, gene-set enrichment and cross-tissue comparison.

## Purpose of this document

This document records the research rationale, analytical decisions, methodological considerations, challenges and interpretations developed during the project. It is a retrospective reconstruction of the research process, intended to document the reasoning behind the analysis rather than simply reproduce its final results.

## 1. Research motivation and question

The project originated from reading the study *A Hybrid Machine Learning and Network Analysis Approach Reveals Two Parkinson's Disease Subtypes from 115 RNA-Seq Post-Mortem Brain Samples*.

I was particularly interested in how an unsupervised machine-learning approach identified two molecularly distinct Parkinson's disease subtypes from transcriptomic data. This raised a further question: whether molecular signatures associated with these brain-defined subtypes could also be detected in peripheral tissues.

The investigation therefore focused on:

**To what extent are brain-derived Parkinson's disease subtype-associated transcriptional signatures represented in peripheral blood and skeletal muscle, and what tissue-specific patterns might influence their correspondence?**

The intention was to explore the feasibility of detecting peripheral representation, not to develop or validate a peripheral subtype classifier.

## 2. Dataset selection and reference material

Three sources formed the basis of the investigation:

- **Brain reference study:** Original subtype-associated gene sets and supplementary differential-expression results.
- **GSE165082:** Peripheral whole-blood RNA-sequencing data.
- **GSE128177:** Skeletal-muscle transcriptomic data.

The original brain study provided the molecular reference. The two peripheral datasets were analysed independently because they represented different tissues and experimental datasets.

An important design constraint was identified early: the peripheral datasets did not contain PDC1/PDC2 subtype labels or matched brain and peripheral samples. Consequently, they could not directly establish whether an individual peripheral sample belonged to either brain-defined subtype.

## 3. Analytical decisions and methodological reasoning

### 3.1 Using brain-derived signatures as reference sets

The original study's subtype-associated gene sets were used to examine their representation in peripheral tissues.

The 42 genes distinguishing PDC1 from PDC2 were considered separately from the six highlighted PDC2-associated inflammatory genes. These gene groups represent different comparisons in the original study and were not treated as interchangeable signatures.

### 3.2 Blood analysis: GSE165082

The dataset contained 26 samples: 12 Parkinson's disease samples and 14 controls.

The analysis involved expression filtering, TMM normalisation, voom transformation and differential-expression testing using limma.

A non-TMM analysis was also retained as a sensitivity comparison to examine whether the main findings depended strongly on the normalisation approach.

No individual genes reached statistical significance after FDR correction in either analysis.

The original study reported adjustment for age, sex and RNA integrity. However, individual-level age, sex and RIN information could not be reconstructed from the available GEO sample metadata. These covariates therefore could not be included in the present reanalysis.

### 3.3 Skeletal-muscle analysis: GSE128177

The dataset contained 36 samples:

- 12 Parkinson's disease samples
- 12 older controls
- 12 young controls

The primary comparison used the 12 Parkinson's disease samples and 12 older controls. Young controls were excluded from this comparison to avoid introducing an additional age-group contrast.

Although the older controls were described as age- and sex-matched at the cohort level, explicit individual pairing identifiers were not available in the GEO metadata. An unpaired design was therefore used rather than assuming sample pairs from their ordering.

The baseline analysis used CPM-based expression filtering, log-transformed expression values and limma empirical Bayes testing.

A TMM–voom sensitivity analysis was subsequently prepared but remains pending validation.

## 4. Findings and cross-tissue interpretation

### 4.1 Brain-derived subtype-associated genes

Of the 42 PDC1-versus-PDC2 genes from the original study, 11 passed the expression filter in each peripheral tissue. None reached FDR-adjusted statistical significance in either tissue.

For the six highlighted PDC2-associated inflammatory genes, none was individually FDR-significant in blood. In muscle, four passed the expression filter, and none was FDR-significant. P2RX7 and C3AR1 did not pass the muscle expression filter and were not tested in the differential-expression analysis.

These results did not provide statistically significant individual-gene evidence that the brain-derived subtype-associated genes were differentially expressed in the peripheral datasets.

### 4.2 Immune-related patterns

Immune-related pathway enrichment was observed in both tissues, but the patterns were heterogeneous rather than uniformly directional.

Exploratory examination of leading-edge genes identified overlaps with PDC2-unique genes from the original brain study. These overlaps were treated as descriptive observations, not evidence of a replicated or conserved PDC2 immune signature.

No formal statistical gene-set overlap test was performed.

### 4.3 Mitochondrial pathway comparison

Six selected mitochondrial pathways were significantly enriched in both tissues.

Their enrichment directions were opposing: the pathways had negative NES values in blood and positive NES values in skeletal muscle.

A comparison of 42 shared leading-edge genes showed a weak negative correlation between tissue-specific log-fold-change estimates (Pearson r = -0.240; Spearman r = -0.238). None of these genes was individually FDR-significant in either tissue.

These findings suggest tissue-dependent differences in the representation of mitochondrial transcriptional patterns. They do not establish a formal tissue interaction or demonstrate differences in mitochondrial function.

## 5. Challenges and limitations identified

Several limitations shaped the interpretation of the results:

- Peripheral samples lacked brain-defined PDC1/PDC2 labels.
- Brain and peripheral samples were not matched by individual.
- Individual-level blood age, sex and RIN covariates were unavailable.
- Individual pairing identifiers were unavailable for the muscle comparison.
- The cohorts were relatively small.
- The tissues originated from different datasets and experimental contexts.
- No independent peripheral cohort was available for replication.
- Non-significant findings could not be interpreted as proof that a biological signature was absent.

These constraints meant that the project could assess exploratory cross-tissue representation but could not validate peripheral biomarkers or establish subtype prediction.

## 6. Further work and current status

- Validate and document the skeletal-muscle TMM–voom sensitivity analysis.
- Continue improving script documentation and reproducibility.
- Review the interpretation of pathway-level findings and their limitations.
- Preserve the distinction between exploratory observations and statistically supported conclusions.

## 7. Research reflection

An important part of this investigation was learning to distinguish an interesting biological pattern from evidence that can actually support a scientific conclusion.

The absence of individually significant differential-expression results did not make the analysis uninformative. It required closer examination of expression filtering, metadata availability, pathway-level findings and the limitations of comparing independently collected tissues.

The project consequently developed beyond a straightforward search for overlapping genes into an exercise in evaluating what public transcriptomic datasets can—and cannot—establish about cross-tissue molecular correspondence.
