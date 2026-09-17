# Parkinson's Subtype Cross-Tissue Analysis

## Overview

This project investigates whether molecular programs that distinguish Parkinson's disease subtypes in the brain are conserved in peripheral tissues, and which features may be tissue-specific.

The project is based on the study:

> A Hybrid Machine Learning and Network Analysis Approach Reveals Two Parkinson’s Disease Subtypes from 115 RNA-Seq Post-Mortem Brain Samples

The original study identified two molecular subtypes of Parkinson's disease, PDC1 and PDC2, using post-mortem brain RNA-seq data.

## Research Question

Which molecular programs that distinguish the two Parkinson's disease subtypes identified in the brain are conserved in peripheral tissues, and which are tissue-specific?

## Starting Hypothesis

We hypothesize that some, but not all, of the molecular programs distinguishing PDC1 and PDC2 in the brain will be detectable in peripheral tissues. The extent of this conservation may vary between tissues, with some subtype-associated programs being shared across tissues and others remaining tissue-specific.

## Data

### Brain

The brain subtype signature was obtained from the supplementary material of the original study.

Supplementary Table S3 contains the differential-expression results comparing PDC1 and PDC2, including the reported 42 subtype-associated genes.

### Peripheral tissue

Publicly available transcriptomic datasets will be evaluated for suitability for cross-tissue subtype analysis.

The first dataset explored was:

- GEO accession: GSE165082
- Tissue: Whole blood
- Technology: RNA-seq

## Analysis Log

### Attempt 1 — Testing the brain signature in whole blood

The 42 genes from the brain PDC1 vs PDC2 comparison were checked against the GSE165082 whole-blood RNA-seq dataset.

This was an exploratory analysis rather than a final test of subtype conservation.

#### What was learned

The brain study compares:

PDC1 vs PDC2

whereas GSE165082 compares:

PD vs healthy controls.

Therefore, these two contrasts do not directly answer the same biological question.

This means that the initial blood comparison cannot be used as the primary evidence for conservation of the PDC1/PDC2 subtype distinction.

#### Decision

This approach was retained as a methodological learning step but will not be treated as the final cross-tissue subtype analysis.

## Current Direction

The next step is to identify peripheral transcriptomic data that allow a more appropriate comparison of molecular variation among Parkinson's disease patients and can therefore be used to investigate the PDC1/PDC2 subtype signature.

The analysis will proceed from:

1. Gene-level comparison
2. Biological pathway/program comparison
3. Network-level comparison, if supported by the available data

## Project Status

**In progress**

The project is being documented as an evolving research process, including unsuccessful approaches, methodological problems, and changes in analytical strategy.
