# 05 eQTL

## Description

This pipeline performs pan-variant eQTL analysis. Specifically, it:

- Quantifies gene expression for each individual
- Extracts the genotype of different variant types
- Performs eQTL mapping, including conditional and fine-mapping analysis

## Requirements

The pipeline is implemented as a Snakemake workflow (`Snakemake >= 7.0`). The required conda environment is described in `envs/env.yml`.

## Inputs

Input files are specified in `config.yaml`.

## Outputs

The pipeline produces the following output files:

- `results/merge.eqtl.txt`: All significant variant-gene pairs
- `results/merge.top.eqtl.txt`: Lead eVariant for significant eGenes
- `results/merge.cond.eqtl.txt`: Conditional eQTL signals for significant eGenes
- `results/merge.susie.eqtl.txt`: Fine-mapped eQTL signals for significant eGenes

## Run

Execute the Snakemake workflow in the current working directory:
```
snakemake -j <cores> --use-conda 
```