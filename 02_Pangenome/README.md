# 02 Pangenome

## Description

This pipeline constructs, filters, and annotates a pangenome graph. Specifically, it:

- Constructs base-level pangenome subgraph in parallel
- Filters repetitive regions and erroneous nodes/edges in the pangenome graph
- Annotate pangenome nodes guided by the assembly paths

## Requirements

The pipeline is implemented as a Snakemake workflow (`Snakemake >= 7.0`). The required conda environment is described in `envs/env.yml`.

## Inputs

Input files are specified in `config.yaml`.

## Outputs

The pipeline produces the following output files:

- `results/subgraph/subgraph_{id}/minigraph_{id}.seqwish.smoothxg.gfaffix.clip.kmer.chop.ids.gfa`: constructed pangenome subgraphs 
- `results/merge.seqwish.smoothxg.gfaffix.clip.kmer.chop.ids.node.anno`: pangenome node annotations

## Run

Execute the Snakemake workflow in the current working directory:
```
snakemake -j <cores> --use-conda 
```