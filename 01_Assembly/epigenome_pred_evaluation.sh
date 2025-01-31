#!/bin/bash
###parameter setting
sample=$1
prefix=$sample.hap1
fa=$prefix.fasta
atac_fq1=${sample}-ATAC.R1.fastq.gz
atac_fq1=${sample}-ATAC.R2.fastq.gz

###Truth
bowtie2-build --threads $thread $fa $prefix
bowtie2 --very-sensitive -k 10 -p $thread -x $prefix -1 $atac_fq1 -2 $atac_fq2 | samtools view -Sbh | samtools sort -O bam -@ $thread -o $prefix.atac.bam
gatk MarkDuplicates -I $prefix.atac.bam -O $prefix.mkDup.bam -M $prefix.mkDup.metrics.txt
samtools view -h -q 30 -F 1804 -f 2 -@ $thread $prefix.mkDup.bam | grep -v chrM | samtools sort -@ $thread -O bam -o $prefix.filtered.bam
samtools index -@ $thread $prefix.filtered.bam
macs3 callpeak -f BAMPE -t $prefix.filtered.bam -g hs -n $prefix -B -q 0.01 --keep-dup all

###Prediction
####SEI
python3 generate_sei_config.py
python3 sei_slurm_submit_script.py
####Enformer
python3 enformer_slurm_pipeline.py --fasta_file $fa --sample_save_dir .
####Compare
python3 epi_benchmark.py input_full.hdf5 sei.atac_dnase.list ${prefix}_peaks.narrowPeak ${prefix}_peaks.noref.narrowPeak prefix 128
python3 epi_benchmark.py input_full.hdf5 enformer.atac_dnase.list ${prefix}_peaks.narrowPeak ${prefix}_peaks.noref.narrowPeak prefix 128
