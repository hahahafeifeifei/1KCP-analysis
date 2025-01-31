

1KCP analysis pipeline
The document record the pipeline of 1KCP analysis
This repository contains all scripts used for the analysis for single-cell eQTL analyses in the brain. The scripts are organized in the following folders:

Description
1_process_expression contains scripts to pseudobulk the single-cell data and generate the input expression matrices for eqtl discovery.
2_process_genotype contains scripts to process the combined vcf genotype data into a readable format.
3_eQTL_discovery contains scripts to perform the eQTL analysis.
4_interaction_modelling contains scripts to perform re-analysis of eQTL results using LME models (and interaction terms).
5_colocalisation contains scripts used for the colocalisation of eQTL discovery outputs and example GWAS.
6_MendelianRandomisation contains scripts for mendelian randomisation following eQTL / GWAS colocalisation
7_pQTL_analysis contains scripts to perform the post-hoc pQTL analysis with UKB-PPP and Robins et al.,(2021) data.


Cite 1KCP