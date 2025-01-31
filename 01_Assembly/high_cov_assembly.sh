#!/bin/bash
###parameter setting
sample=$1
thread=16
reference=hg38.no_alt.fa
adaptor=pacbio.adapter.fa
refseq_dir=Refseq
nonhuman_accession=non-human.accession.list
hic_fq1=${sample}-HiC.R1.fastq.gz
hic_fq2=${sample}-HiC.R2.fastq.gz
hifi_fq=${sample}-HiFi.fastq.gz

###Hifiasm assembly
hifiasm -o $sample.hifiasm -t $thread --h1 $hic_fq1 --h2 $hic_fq2 $hifi_fq
grep ^S $sample.hifiasm.bp.hap1.p_ctg.gfa | awk '{print ">"$2"\n"$3}' > $sample.hifiasm.hap1.fasta
grep ^S $sample.hifiasm.bp.hap2.p_ctg.gfa | awk '{print ">"$2"\n"$3}' > $sample.hifiasm.hap2.fasta

###Unicycler assembly (chrM)
minimap2 -t $thread -x map-hifi -a -Y -L --eqx --cs $reference $hifi_fq | samtools view -Shb | samtools sort -@ $thread -o $sample.bam -O bam
samtools index -@ $thread $sample.bam
samtools view -F 256 -Shb $sample.bam chrM | samtools fastq > $sample.chrM.fastq
unicycler-runner.py -l $sample.chrM.fastq -o $sample

###Merge chrM assembly
minimap2 -t $thread -x asm5 -a -Y -L --eqx --cs $reference $sample.hifiasm.hap1.fasta | samtools view -S -F 256 | awk '{if($3=="chrM") print $1}' > $sample.hap1.chrM.contig
seqkit grep -v -f $sample.hap1.chrM.contig $sample.hifiasm.hap1.fasta | cat - assembly.fasta > $sample.merge.hap1.fasta
minimap2 -t $thread -x asm5 -a -Y -L --eqx --cs $reference $sample.hifiasm.hap2.fasta | samtools view -S -F 256 | awk '{if($3=="chrM") print $1}' > $sample.hap2.chrM.contig
seqkit grep -v -f $sample.hap2.chrM.contig $sample.hifiasm.hap2.fasta > $sample.merge.hap2.fasta

###Masking adaptor sequence
minimap2 -t $thread -cxsr -f5000 -N2000 -secondary=yes --cs $sample.merge.hap1.fasta $adaptor | awk -v OFS='\t' '{print$6,$8,$9}' | bedtools sort -i - | bedtools merge -i - > $sample.adaptor.hap1.bed
bedtools maskfasta -fi $sample.adaptor.hap1.bed $sample.merge.hap1.fasta > $sample.masked.hap1.fasta
minimap2 -t $thread -cxsr -f5000 -N2000 -secondary=yes --cs $sample.merge.hap2.fasta $adaptor | awk -v OFS='\t' '{print$6,$8,$9}' | bedtools sort -i - | bedtools merge -i - > $sample.adaptor.hap2.bed
bedtools maskfasta -fi $sample.adaptor.hap2.bed $sample.merge.hap2.fasta > $sample.masked.hap2.fasta

###Remove non-human contigs
> $sample.non-human.hap1.contig
> $sample.non-human.hap2.contig
for type in {bacteria,viruses_and_viroids,archaea,fungi}
do
    blastn -query $sample.masked.hap1.fasta -db $refseq_dir/$type/$type -task megablast -word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 0.0001 -min_raw_gapped_score 100 -penalty -5 -perc_identity 98.0 -soft_masking true -out $sample.$type.hap1.blast -num_threads $thread -outfmt 7
    grep -v "^#" $sample.$type.hap1.blast | csvtk -H -t join -f "2;1" - $nonhuman_accession | awk '{print $1}' >> $sample.non-human.hap1.contig
    blastn -query $sample.masked.hap2.fasta -db $refseq_dir/$type/$type -task megablast -word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 0.0001 -min_raw_gapped_score 100 -penalty -5 -perc_identity 98.0 -soft_masking true -out $sample.$type.hap2.blast -num_threads $thread -outfmt 7
    grep -v "^#" $sample.$type.hap2.blast | csvtk -H -t join -f "2;1" - $nonhuman_accession | awk '{print $1}' >> $sample.non-human.hap2.contig
done
seqkit grep -v -f $sample.non-human.hap1.contig $sample.masked.hap1.fasta > $sample.hap1.fasta
seqkit grep -v -f $sample.non-human.hap2.contig $sample.masked.hap2.fasta > $sample.hap2.fasta