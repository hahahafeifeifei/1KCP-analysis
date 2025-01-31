#!/bin/bash
###parameter setting
sample=$1
thread=16
star_index=STAR_genome_GRCh38_v40
gtf=gencode.v40.GRCh38.annotation.gtf
rna_fq1=${sample}-RNA.R1.fastq.gz
rna_fq2=${sample}-RNA.R2.fastq.gz

###RNA-Seq read alignment
STAR --runMode alignReads \
--runThreadN $thread \
--genomeDir $star_index \
--twopassMode Basic \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.1 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outFilterType BySJout \
--outFilterScoreMinOverLread 0.33 \
--outFilterMatchNminOverLread 0.33 \
--limitSjdbInsertNsj 1200000 \
--readFilesIn ${rna_fq1} ${rna_fq2} \
--readFilesCommand zcat \
--outFileNamePrefix ${sample}_ \
--outSAMstrandField intronMotif \
--outFilterIntronMotifs None \
--alignSoftClipAtReferenceEnds Yes \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMtype BAM Unsorted \
--outSAMunmapped Within \
--genomeLoad NoSharedMemory \
--chimSegmentMin 15 \
--chimJunctionOverhangMin 15 \
--chimOutType WithinBAM SoftClip \
--chimMainSegmentMultNmax 1 \
--outSAMattributes NH HI AS nM NM ch \
--outSAMattrRGline ID:rg1 SM:sm1
samtools sort -@ $thread -o ${sample}_Aligned.sortedByCoord.out.bam ${sample}_Aligned.out.bam
picard MarkDuplicates I=${sample}_Aligned.sortedByCoord.out.bam O=${sample}_Aligned.sortedByCoord.md.out.bam M=${sample}_marked_dup_metrics.txt

###Gene quantification
rnaseqc.v2.4.2.linux $gtf ${sample}_Aligned.sortedByCoord.md.out.bam . -s $sample -vv --stranded=rf
