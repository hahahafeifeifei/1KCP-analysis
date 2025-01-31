#!/bin/bash
###parameter setting
thread=16
tr=vamos.motif.hg38.v2.1.e0.1.nohomo.bed
reference=GRCh38_no_alt_analysis_set.fasta
assembly_list=assembly.list

###Individual TR detection
> tr.list
cat $assembly_list | while read assembly
do
    prefix=$(ls $assembly | cut -d "." -f 1-2)
    sample=$(ls $assembly | cut -d "." -f 1)
    hap=$(ls $assembly | cut -d "." -f 2)
    minimap2 -t $thread -x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 -a $reference $assembly | samtools view -Shb | samtools sort -@ $thread -o $prefix.bam
    vamos --contig -b $prefix.bam -r $reference -o $prefix.vamos.tr.vcf -L 50000 -s $prefix -t $thread
    ls $prefix.vamos.tr.vcf >> tr.list
done

###Merge TR callset
tryvamos.py combineVCF tr.list 1kcp.vamos.tr.merge.raw.vcf
cat 1kcp.vamos.tr.merge.raw.vcf | python3 graph_vcf_merge_split.py merge > 1kcp.vamos.tr.merge.vcf 
