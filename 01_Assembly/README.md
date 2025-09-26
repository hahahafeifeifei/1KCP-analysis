# 01 Assembly

## Description

This pipeline performs diploid genome assembly, annotation, and quality assessment for each sample. Specifically, it:

- Generates the diploid genome assemblies
- Annotates repeat, gene and epigenomic elements in the assemblies
- Assesses the quality of assemblies

## Requirements

The pipeline is implemented as a Snakemake workflow (`Snakemake >= 7.0`). The required conda environment is described in `envs/env.yml`.

## Inputs

Input files are specified in `config.yaml`.

## Outputs

The pipeline produces the following output files:

- `results/{sample}/{sample}.{haplotype}.fasta`: assembly sequence
- `results/{sample}/{sample}.{haplotype}.repeatmasker.out`: RepeatMasker repeat annotation
- `results/{sample}/{sample}.{haplotype}.trf.dat.gz`: TRF repeat annotation
- `results/{sample}/{sample}.{haplotype}.liftoff.gff_polished`: Liftoff gene annotation
- `results/{sample}/{sample}.{haplotype}.stringtie.annotated.gtf`: Stringtie gene annotation
- `results/{sample}/{sample}.{haplotype}.sei.bed`: SEI epigenome annotation
- `results/{sample}/{sample}.{haplotype}.type.final.bed`: integrated annotation from repeat, gene, and epigenomo annotations
- `results/{sample}/{sample}.merqury`: merqury assessment result
- `results/{sample}/{sample}.{haplotype}.inspector`: inspector assessment result

## Run

Execute the Snakemake workflow in the current working directory:
```
snakemake -j <cores> --use-conda 
```