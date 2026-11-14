# Evaluating Peripheral Transcriptomic Signatures Associated with Brain-Defined Parkinson's Disease Subtypes

## Overview

This project investigates whether molecular signatures associated with two Parkinson's disease (PD) subtypes identified in post-mortem brain tissue are detectable in peripheral blood and skeletal muscle.

The analysis uses publicly available transcriptomic datasets to explore the representation of brain-derived subtype-associated genes and biological pathways in peripheral tissues, with particular attention to immune-related and mitochondrial processes.

The project is an exploratory cross-tissue transcriptomic analysis. It does not aim to classify patients into PDC1 or PDC2, establish clinical biomarkers, or demonstrate subtype-specific treatment responses.

## Background

The project is based on the study:

> A Hybrid Machine Learning and Network Analysis Approach Reveals Two Parkinson's Disease Subtypes from 115 RNA-Seq Post-Mortem Brain Samples.

The original study identified two molecular subtypes, PDC1 and PDC2, using post-mortem brain RNA-sequencing data.

The study reported:
- 42 differentially expressed genes distinguishing PDC1 from PDC2.
- Subtype-associated gene signatures from comparisons with healthy controls.
- Differences in biological processes associated with the two subtypes, including immune and inflammatory processes in PDC2.

These brain-derived signatures were used as reference gene sets for the peripheral tissue analyses.

**Original study:** https://pmc.ncbi.nlm.nih.gov/articles/PMC8910747/

## Research Question

To what extent are brain-derived molecular signatures associated with Parkinson's disease subtypes detectable in peripheral blood and skeletal muscle, and what tissue-specific transcriptional patterns may influence their correspondence?

## Hypothesis

Brain-derived Parkinson's disease subtype-associated transcriptional signatures are expected to show limited and potentially tissue-dependent representation in peripheral blood and skeletal muscle, reflecting differences in tissue-specific immune and mitochondrial gene-expression patterns.

## Objectives

1. Characterise the brain-derived PDC1 and PDC2 molecular signatures.
2. Evaluate the representation of these signatures in peripheral blood.
3. Evaluate their representation in skeletal muscle.
4. Investigate tissue-dependent immune-related transcriptional patterns.
5. Compare mitochondrial pathway enrichment between blood and muscle.
6. Assess the feasibility and limitations of peripheral transcriptomic signatures as potential molecular surrogates of brain-defined PD subtypes.

## Datasets

### 1. Brain reference signatures

The brain-derived gene sets were obtained from the supplementary material of the original study.

- Supplementary Table S3: PDC1 versus PDC2 differential-expression results, including 42 subtype-associated genes.
- Supplementary Table S4: subtype-associated differential-expression results from comparisons with healthy controls.
- The PDC2-unique gene set contains 482 genes.

These gene sets were used as reference signatures rather than as peripheral subtype labels.

### 2. Peripheral blood

**GEO accession:** GSE165082
**Tissue:** Whole blood
**Technology:** RNA sequencing

The analysed count matrix contained:
- 12 Parkinson's disease samples.
- 14 control samples.
- 26 samples in total.

The count data were filtered, normalised using the trimmed mean of M-values (TMM) method, transformed using voom, and analysed with limma.

### 3. Skeletal muscle

The analysed count matrix contained 36 skeletal muscle samples:

- 12 Parkinson's disease samples.
- 12 older control samples.
- 12 young control samples.

The primary differential-expression comparison included 12 Parkinson's disease samples and 12 older controls. Young controls were excluded to focus on the comparison between PD and age-comparable older adults.

The study describes the older control group as age- and sex-matched to the PD group. However, explicit individual pair identifiers were not available in the GEO sample metadata. The primary analysis therefore used an unpaired group comparison.

**GEO accession:** GSE128177
**Tissue:** Skeletal muscle

## Analytical Workflow

### 1. Data acquisition and quality control

- Obtained publicly available transcriptomic datasets.
- Examined sample metadata, expression matrices and dataset structure.
- Assessed sample numbers, expression distributions and data quality.
- Applied expression-based filtering before differential-expression analysis.

### 2. Differential expression analysis

Differential-expression analyses were conducted independently for blood and skeletal muscle.

**Blood:**
- CPM-based expression filtering.
- TMM normalisation.
- Voom transformation.
- Limma differential-expression modelling.
- Non-TMM analysis used as a sensitivity comparison.

**Muscle:**
- CPM-based expression filtering.
- Log2(CPM + 1) transformation.
- Limma modelling with empirical Bayes moderation.

Adjusted p-values were used to account for multiple testing. An FDR threshold of 0.05 was used to identify statistically significant results.

### 3. Brain signature representation

The brain-derived PDC1, PDC2 and shared gene sets were evaluated against the peripheral transcriptomic datasets.

The analysis included:
- Detection of brain subtype-associated genes after peripheral expression filtering.
- Comparison of differential-expression estimates.
- Gene-set enrichment analysis to investigate whether brain-derived signatures showed significant enrichment in peripheral tissues.

### 4. Biological pathway analysis

Gene-set enrichment analysis was performed using:
- Gene Ontology Biological Process (GO BP).
- Hallmark gene sets.
- Reactome pathways.

The analysis used ranked gene lists based on signed moderated t-statistics and the fgsea framework.

Immune-related pathways were examined separately to investigate the representation of PDC2-associated genes in peripheral immune-related pathway results.

### 5. Cross-tissue mitochondrial comparison

Mitochondrial GO BP pathways showing significant enrichment in both blood and skeletal muscle were examined.

The comparison included:
- Normalised enrichment scores (NES).
- Adjusted p-values.
- Shared leading-edge genes.
- Cross-tissue expression estimates of shared leading-edge genes.

This analysis was used to investigate tissue-dependent differences in enrichment direction.

## Key Findings

### 1. Limited representation of brain-derived subtype signatures

Direct gene-set enrichment analysis did not identify statistically significant enrichment of the brain-derived PDC1, PDC2 or shared signatures in either peripheral tissue after multiple-testing correction.

Of the 42 genes distinguishing PDC1 from PDC2, 11 passed expression filtering in each peripheral tissue. None was individually significant after FDR correction.

The six highlighted PDC2-associated inflammatory genes showed no individually significant differential expression after FDR correction in blood. In skeletal muscle, four genes passed the expression filter, and none was FDR-significant; P2RX7 and C3AR1 did not pass the expression filter and were therefore not included in the muscle differential-expression testing.

### 2. Immune-related pathway findings

Significant immune-related GO BP pathways were identified in both peripheral tissues:

- Blood: 9 significant immune-related pathways.
- Skeletal muscle: 13 significant immune-related pathways.

The leading-edge genes of these pathways were compared with the PDC2-unique brain gene set.

- Blood: 13 overlapping genes.
- Skeletal muscle: 20 overlapping genes.
- Four genes were shared between the two peripheral tissue overlap sets.

None of the overlapping genes was individually significant after FDR correction in either tissue.

These findings represent descriptive gene-set overlap, not statistically established replication of the brain-derived PDC2 immune signature.

Immune-related pathway enrichment also showed heterogeneous directions across the two tissues. The results should not be interpreted as evidence of global immune activation or suppression.

### 3. Contrasting mitochondrial pathway enrichment

Six selected mitochondrial GO BP pathways showed significant enrichment in both blood and skeletal muscle.

All six displayed contrasting enrichment directions:
- Negative NES in blood.
- Positive NES in skeletal muscle.

A total of 42 shared leading-edge genes were identified across the selected pathways. These genes showed negative log2 fold-change estimates in blood and positive estimates in muscle.

The correlation between tissue-specific estimates was weakly negative (Pearson r = -0.240).

None of the 42 shared genes was individually significant after FDR correction in either tissue.

These findings suggest tissue-dependent mitochondrial transcriptomic patterns, but do not establish differences in mitochondrial functional activity.

## Figures

The main figures are available in `results/figures/`.

1. `mitochondrial_pathway_cross_tissue_comparison.png` — Comparison of enrichment scores for six mitochondrial pathways across blood and skeletal muscle.
2. `mitochondrial_shared_gene_logFC_scatter.png` — Comparison of expression estimates for 42 shared mitochondrial leading-edge genes across the two tissues.

## Repository Structure

```text
parkinsons-subtype-cross-tissue-analysis/
│
├── analysis/
│   ├── 01_data_audit.py
│   ├── 02_quality_control.py
│   ├── 03_muscle_differential_expression.py
│   ├── 04_brain_signature.py
│   ├── 05_signature_overlap.py
│   ├── 06_signature_muscle_analysis.py
│   ├── 07_muscle_DE_final.py
│   ├── 08_muscle_DE_limma.R
│   ├── 09_signature_enrichment.py
│   ├── 10_blood_data_audit.py
│   ├── 11_blood_quality_control.py
│   ├── 12_blood_metadata_audit.py
│   ├── 13_blood_DE_limma_voom.R
│   ├── 14_signature_gene_matching_audit.py
│   ├── 15_blood_TMM_sensitivity.R
│   ├── 16_blood_GSEA.R
│   ├── 17_muscle_GSEA.R
│   ├── 18_mitochondrial_leading_edge_integration.R
│   ├── 19_brain_signature_cross_tissue.R
│   └── 20_blood_signature_enrichment.py
│
├── data/
│   └── brain_signature/
│       └── brain_gene_sets.csv
│
├── results/
│   ├── blood_GSEA/
│   ├── muscle_GSEA/
│   ├── figures/
│   ├── GSE128177_muscle_DE_limma.csv
│   └── GSE165082_blood_DE_TMM_limma_voom.csv
│
└── README.md
```

*Note: This is a high-level representation of the repository. Additional data and result files are present in the corresponding directories.*

## Researcher’s Contribution

This project was self-directed as an exploratory bioinformatics investigation. I formulated the research question around the extent to which brain-defined Parkinson’s disease subtype signatures may be represented in peripheral tissues.

I examined the original subtype study and its supplementary gene sets, identified publicly available blood and skeletal muscle transcriptomic datasets, and conducted metadata assessment, expression filtering, differential-expression analysis and gene-set enrichment analysis using Python and R-based workflows.

I then integrated the results across tissues, examining immune-related pathway overlaps and contrasting mitochondrial enrichment patterns. An important part of the work involved evaluating non-significant findings, identifying limitations in dataset metadata and sample comparability, and distinguishing descriptive overlap from statistically supported evidence.

The project documents the analytical decisions, limitations and interpretation of results, including approaches explored during its development.

## Limitations

1. **Absence of peripheral subtype labels:** Neither peripheral dataset provides PDC1/PDC2 classifications. Therefore, direct peripheral subtype classification and validation were not possible.

2. **No matched brain and peripheral samples:** The datasets do not contain matched brain, blood and muscle samples from the same individuals.

3. **Limited sample sizes:** The analysed peripheral cohorts are relatively small, which may limit statistical power.

4. **Blood covariates:** Individual age, sex and RNA integrity number (RIN) information could not be reconstructed for the analysed blood count matrix. Consequently, the reanalysis could not adjust for these covariates, unlike the original study's reported model.

5. **Tissue-specific filtering:** Different expression-filtering and preprocessing approaches were used for blood and muscle. Cross-tissue comparisons therefore require cautious interpretation.

6. **No independent peripheral replication cohort:** The current analysis does not include an independent dataset for replication.

7. **Exploratory pathway overlap:** Overlap between PDC2-unique genes and immune-related pathway leading-edge genes was descriptive. A formal statistical gene-set overlap test was not performed.

8. **Transcriptomic rather than functional evidence:** Gene-expression and pathway-enrichment results do not directly measure mitochondrial function, immune-cell activity or biological mechanisms.

9. **No clinical validation:** The findings do not establish diagnostic biomarkers, subtype-specific treatment responses or clinical utility.

## Conclusion

This exploratory analysis found limited statistically supported representation of brain-derived Parkinson's disease subtype signatures in peripheral blood and skeletal muscle.

Although immune-related pathways showed descriptive overlap with PDC2-unique brain genes, none of the overlapping genes was individually significant after multiple-testing correction.

In contrast, six selected mitochondrial pathways showed significant but opposing enrichment directions across blood and skeletal muscle, accompanied by weakly negative correlation among shared leading-edge gene-expression estimates.

Overall, the findings highlight the importance of tissue context when investigating peripheral transcriptomic patterns in Parkinson's disease. They also demonstrate the limitations of using peripheral datasets without subtype labels or matched brain samples to assess the correspondence of brain-defined molecular signatures.

Further work using larger independent cohorts, matched tissue samples and clinically characterised PD subtypes would be required to establish whether peripheral transcriptomic signatures can serve as reliable molecular surrogates.

## Project Status

**Status: Exploratory analysis completed; documentation and interpretation in progress.**

This repository documents the analytical workflow, methodological decisions, exploratory findings and limitations of an independent bioinformatics research project.

The project is intended to be transparent and reproducible, including documentation of analytical decisions and approaches that were not retained as primary analyses.