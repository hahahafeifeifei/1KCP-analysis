# 04 Medical

## Description

This pipeline performs medically relevant genomic analyses. Specifically, it:

- Performs gene cluster analysis using Pangene

## Requirements

The pipeline is implemented as a Snakemake workflow (`Snakemake >= 7.0`). The required conda environment is described in `envs/env.yml`.

## Inputs

Input files are specified in `config.yaml`.

## Outputs

The pipeline produces the following output files:

- `results/merge.pangene.gfa`: Pangene gene-level graph

## Run

Execute the Snakemake workflow in the current working directory:
```
snakemake -j <cores> --use-conda 
```