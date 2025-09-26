# 06 Imputation

## Description

This pipeline performs genotype imputation based on a reference panel. Specifically, it:

- Performs leave-one-out imputation, where each sample is removed from the reference panel in turn to evaluate imputation accuracy

## Requirements

The pipeline is implemented as a Snakemake workflow (`Snakemake >= 7.0`). The required conda environment is described in `envs/env.yml`.

## Inputs

Input files are specified in `config.yaml`.

## Outputs

The pipeline produces the following output files:

- `results/{sample}/{sample}.{chr}.impute.vcf.gz`: imputed genotypes for the left-out sample (per chromosome {chr})

## Run

Execute the Snakemake workflow in the current working directory:
```
snakemake -j <cores> --use-conda 
```