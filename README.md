# Landscape genomics highlights the adaptive evolution of chickpea across the Silks Roads 

# Landscape Genomics Analysis Pipeline

## Overview

This repository contains the analysis pipeline used to process genomic variant data, characterize population structure, perform landscape genomic analyses, and identify candidate genes and biological functions associated with environmental adaptation.

The workflow is organized into sequential steps, starting from raw VCF filtering and preparation, followed by population structure and phylogenetic analyses, environmental data processing, landscape genomic analyses, and Gene Ontology (GO) enrichment analysis.

## Pipeline structure

```text
01_01_VCF_file_filtering_missing_0.9.txt
        │
        ▼
01_02vcf preparation_landscape.txt
        │
        ├──────────────────────┐
        ▼                      ▼
02_01_plotting_pop_structure.R  02_02plotting_phylogenetic_tree.R
        │                      │
        └──────────┬───────────┘
                   ▼
03_01_readscript_Worldclim.txt
                   │
                   ▼
03_02_landscape_analysis.R
                   │
                   ▼
03_03_GO_analysis.R
```

## 01. VCF preparation

### `01_01_VCF_file_filtering_missing_0.9.txt`

This step performs the initial filtering of the raw VCF dataset.

The main objective is to generate a high-quality SNP dataset by applying missing-data filtering. The script is configured around a **90% genotype call threshold**, retaining variants with sufficient genotype information across individuals.

**Input:**

* Raw VCF file

**Output:**

* Filtered VCF containing SNPs passing the missing-data threshold

---

### `01_02vcf preparation_landscape.txt`

This step prepares the filtered VCF for downstream landscape genomic analyses.

The resulting dataset is used as the genomic input for population structure, phylogenetic analyses, and genotype–environment association analyses.

**Input:**

* Filtered VCF from `01_01`

**Output:**

* Landscape-genomics-ready VCF/genotype dataset

---

## 02. Population structure and phylogeny

### `02_01_plotting_pop_structure.R`

This script analyzes and visualizes the genetic population structure of the sampled individuals.

The results provide information on genetic clustering and population differentiation, which can subsequently be considered when interpreting genotype–environment associations.

**Input:**

* Prepared genomic dataset from Section 01

**Output:**

* Population structure results
* Population structure plots/figures

---

### `02_02plotting_phylogenetic_tree.R`

This script generates a phylogenetic/tree-based representation of the genetic relationships among individuals or populations.

The resulting tree provides a complementary visualization of genetic relationships and population structure.

**Input:**

* Prepared genomic dataset from Section 01

**Output:**

* Phylogenetic tree
* Tree visualization/figure

---

## 03. Environmental and landscape genomic analyses

### `03_01_readscript_Worldclim.txt`

This step prepares the environmental variables used for the landscape genomic analysis.

The script is associated with the extraction/processing of **WorldClim environmental data**, generating environmental predictors corresponding to the sampling locations.

**Input:**

* Sampling-location information
* WorldClim climate data

**Output:**

* Environmental dataset associated with sampling locations

---

### `03_02_landscape_analysis.R`

This is the main landscape genomic analysis step.

The script combines genomic information with environmental variables to identify genetic variants associated with environmental gradients.

The analysis includes the preparation of genomic and environmental matrices and downstream genotype–environment association analyses.

**Input:**

* Prepared VCF/genomic dataset
* Environmental variables generated in `03_01`
* Population structure information, where required

**Output:**

* Candidate SNPs associated with environmental variables
* Environmental association statistics
* Candidate genomic regions/markers for downstream functional analysis

---

## 04. Functional annotation and GO enrichment

### `03_03_GO_analysis.R`

This script performs functional annotation and Gene Ontology enrichment analysis of candidate genes identified from the landscape genomic analysis.

The analysis is used to determine whether candidate genes are significantly enriched for particular biological processes, molecular functions, or cellular components.

The workflow includes:

1. Linking candidate genomic markers to genes.
2. Retrieving or assigning GO annotations.
3. Defining an appropriate genomic/background gene set.
4. Testing GO-term enrichment.
5. Correcting enrichment P-values for multiple testing.
6. Identifying significantly enriched GO categories.
7. Generating summary tables and visualization of enriched biological functions.

**Input:**

* Candidate SNPs/markers from the landscape analysis
* Gene annotation
* GO annotation
* Background gene set

**Output:**

* GO enrichment results
* Significantly enriched GO terms
* Candidate gene functional summaries
* GO enrichment plots

---



## Summary

This pipeline provides a complete workflow from **raw genomic variant filtering to functional interpretation of candidate genes**. It integrates genomic data, population structure, environmental information, genotype–environment association analyses, and GO enrichment to investigate the genetic basis of environmental adaptation.
