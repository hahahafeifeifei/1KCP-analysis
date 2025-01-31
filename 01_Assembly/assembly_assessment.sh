#!/bin/bash
###parameter setting
sample=$1
thread=16
hap1_fa=$sample.hap1.fasta
hap2_fa=$sample.hap2.fasta
sr_fq1=${sample}-SR.R1.fastq.gz
sr_fq2=${sample}-SR.R2.fastq.gz
lr_fq=${sample}-LR.fastq.gz

###Meryl (short read) assessment
meryl count k=21 memory=50G threads=${thread} output $sample.R1.meryl $sr_fq1
meryl count k=21 memory=50G threads=${thread} output $sample.R2.meryl $sr_fq2
meryl union-sum output $sample.WGS.meryl $sample.R1.meryl $sample.R2.meryl
merqury.sh $sample.WGS.meryl $hap1_fa $hap2_fa $sample.merqury

###Inspector (long read) assessment
inspector.py -c $hap1_fa -r $lr_fq --skip_base_error_detect --skip_base_error -t $thread -o $sample.hap1.inspector
inspector.py -c $hap2_fa -r $lr_fq --skip_base_error_detect --skip_base_error -t $thread -o $sample.hap2.inspector