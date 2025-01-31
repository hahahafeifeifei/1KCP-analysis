#!/bin/bash
###parameter setting
thread=16
reference=GRCh38_no_alt_analysis_set.fasta
assembly_list=assembly.list
samples=1kcp.sample

###Individual HLA calling
cat $assembly_list | while read assembly
do
    prefix=$(ls $assembly | cut -d "." -f 1-2)
    sample=$(ls $assembly | cut -d "." -f 1)
    hap=$(ls $assembly | cut -d "." -f 2)
    minimap2 -t $thread -x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 -a $reference $assembly | samtools view -Shb | samtools sort -@ $thread -o $prefix.bam
    samtools index $prefix.bam 
    hifihla call-contigs --abam $prefix.bam --hap1 $assembly --outdir tmp
    awk 'NR>1{print$0}' tmp/hifihla_summary.tsv > $prefix.hla.typing
done

###Merge HLA callset
cat */*hla.typing | awk '{if($4!="N/A") print $1"\t"$2"\t"$3"\t"$4"\t4"}' | awk '{for(i=1;i<=$5;i++) {split($4,a,"*");split(a[2],b,":"); printf $1"\t"$2"\t"$3"\t"a[1]"*"b[1];max=i;if(i>length(b)) max=length(b); for(j=2;j<=max;j++) printf ":"b[j]; print "\t"i }}' > 1kcp.hifihla.merge.field_decompose
awk '{print $1"."$2"\t"$3"\t"$4"\t"$5}' 1kcp.hifihla.merge.field_decompose | csvtk -H -t fold -f 2,3,4 -v 1 -s ","|  awk '{split($4,a,",");print $1"\t"$2"\t"$3"\t"length(a)"\t"$4}' | csvtk -H -t sort -k 1 -k 3:n -k 4:nr - | csvtk -H -t join -f "1;3" - gene.pos | awk -v OFS='\t' '{print $6,$7,$1,$2,$3,$4,$5}' > 1kcp.hifihla.merge.field_decompose.geno
python3 mhc_vcf.py 1kcp.hifihla.merge.field_decompose.geno $samples | bcftools sort | bcftools plugin fill-tags -t AN,AC,AF > 1kcp.hifihla.vcf

###eLD calculation
bcftools view -i "Field==4 && AF>0.01 && AF<=0.99" 1kcp.hifihla.vcf | bcftools query -f "%Allele[\t%GT]\n" | awk 'BEGIN{for(j=1;j<=5;j++) {printf N;for(i=2;i<=2234;i++) printf "N ";print ""}} {printf "M "$1;for(i=2;i<=NF;i++) {split($i,a,"|"); if(a[1]=="1")printf " P";else printf " A";if(a[2]=="1")printf " P";else printf " A" };print ""}' | grep -v "HLA-Y" > 1kcp.4field.bgl.phased
Rscript eLD.v1.0.R
bcftools view -i "Field==3 && AF>0.01 && AF<=0.99" 1kcp.hifihla.vcf | bcftools query -f "%Allele[\t%GT]\n" | awk 'BEGIN{for(j=1;j<=5;j++) {printf N;for(i=2;i<=2234;i++) printf "N ";print ""}} {printf "M "$1;for(i=2;i<=NF;i++) {split($i,a,"|"); if(a[1]=="1")printf " P";else printf " A";if(a[2]=="1")printf " P";else printf " A" };print ""}' | grep -v "HLA-Y" > 1kcp.3field.bgl.phased
Rscript eLD.v1.0.R