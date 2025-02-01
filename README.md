# 1KCP Analysis

This repository contains the analysis scripts used for the **1000 Chinese Pangenome (1KCP) project**, which generated diploid genome assemblies for 1,116 Chinese individuals.

## Script Descriptions

The scripts are organized into the following directories:

### 01_Assembly

- `high_cov_assembly.sh`: Generates genome assemblies using the 1KCP high-coverage dataset.
- `assembly_assessment.sh`: Assesses the quality of the 1KCP assemblies using both short-read and long-read data.
- `variant_assessment.sh`: Evaluates the variant calling and haplotype reconstruction performance of the PIGA assemblies.
- `epigenome_pred_evaluation.sh`: Compares the epigenome annotation performance of the SEI and Enformer models against ATAC-seq data.
- `assembly_annotation.sh`: Annotates genomic elements in the assemblies, including repeats, genes, and epigenome elements.

### 02_Pangenome

- `minigraph.sh`: Constructs the Minigraph pangenome.
- `subgraph_split.sh`: Splits the pangenome, alignments, and sequences into subgraphs.
- `subgraph_construction.sh`: Constructs, normalizes, and clips the base-level pangenome for each subgraph.
- `pangenome_assessment.sh`: Assesses the quality of the pangenome through path concordance and haplotype coverage.
- `pangenome_annotation.sh`: Generates pangenome annotations based on the path-guided annotation method.

### 03_Variant

- `graph_variant_detection.sh`: Generates variant callsets (small variants, structural variants, and nested variants) based on the pangenome decomposition method.
- `tandem_repeat_detection.sh`: Generates the tandem repeat callset using vamos.
- `variant_annotation.sh`: Annotates the function of nested variants and the repetitive content of structural variant sites based on the pangenome.

### 04_Medical Relevance

- `tr_expansion_detection.sh`: Detects tandem repeat expansion events based on the DBSCAN algorithm.
- `pangene.sh`: Constructs the Pangene graph and detects gene cluster variations.
- `hla_allele_detection.sh`: Performs HLA gene haplotype analysis.

### 05_eQTL

- `expression_quantification.sh`: Quantifies gene expression for each individual.
- `eQTL.sh`: Prepares variant and gene expression data and performs eQTL mapping.
- `colocalization.sh`: Conducts eQTL-GWAS colocalization analysis.

### 06_Imputation

- `colocalization.sh`: Performs leave-one-out imputation analysis.

## Cite