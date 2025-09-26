# 03 Variant

## Description

This pipeline builds an extensive variant catalog from an assembly panel. Specifically, it:

- Detects graph variants (SNVs, indels, SVs, and nested variants) from the pangenome graph
- Identifies TR alleles using vamos
- Performs HLA typing using HifiHLA
- Annotates SVs and nested variant based on the pangenome annotation

## Requirements

The pipeline is implemented as a Snakemake workflow (`Snakemake >= 7.0`). The required conda environment is described in `envs/env.yml`.

## Inputs

Input files are specified in `config.yaml`.

## Outputs

The pipeline produces the following output files:

- `results/merge.ref.small.allele.reliable.vcf.gz`: small variant callset
- `results/merge.ref.sv.allele.reliable.vcf.gz`: SV callset
- `results/merge.nest.allele.reliable.vcf.gz`: nested variant callset
- `results/merge.tr.site.reliable.vcf`: TR site callset
- `results/merge.tr.length.reliable.vcf`: TR length variant callset
- `results/merge.tr.motif.reliable.vcf`: TR motif variants callset
- `results/merge.hla.reliable.vcf`: HLA allele callset
- `results/merge.ref.sv.repeat.anno`: Repeat composition annotation of SVs
- `results/merge.ref.nest.flanking.anno`: Flanking annotation of nested variants

## Run

Execute the Snakemake workflow in the current working directory:
```
snakemake -j <cores> --use-conda
```