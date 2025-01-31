#!/bin/bash
###parameter setting
thread=16
prefix=1kcp
sample=${prefix}.rna.sample
gencode=gencode.v40.GRCh38.gene.annotation.bed

###Prepare genotype data (Non-TR)
bcftools plugin fill-tags --threads $thread 1kcp.small.vcf | bcftools view -i 'AN>=1786 && HWE>1e-6 && (AF>=0.01 && AF<=0.99)' -o 1kcp.small.QC.vcf.gz
tabix 1kcp.small.QC.vcf.gz
bcftools plugin fill-tags --threads $thread 1kcp.sv.vcf | bcftools view -i 'AN>=1786 && HWE>1e-6 && (AF>=0.01 && AF<=0.99)' -o 1kcp.sv.QC.vcf.gz
tabix 1kcp.sv.QC.vcf.gz
bcftools plugin fill-tags --threads $thread 1kcp.nest.vcf | bcftools view -i 'AN>=1786 && HWE>1e-6 && (AF>=0.01 && AF<=0.99)' -o 1kcp.nest.QC.vcf.gz
tabix 1kcp.nest.vcf.gz

###Prepar genotype data (TR)
awk -v OFS='\t' '{if(substr($1,1,1)=="#"){if(substr($1,1,2)=="#C")print "##INFO=<ID=VID,Number=1,Type=Integer,Description=\"VNTR ID\">"} else {$8=$8";VID="NR};print $0}' 1kcp.vamos.tr.merge.vcf  | sed "s/source=vamos_2.1.4/source=adVNTR/g" > 1kcp.vamos.tr.merge.modified.vcf 
dumpSTR --min-locus-callrate 0.8 --min-locus-hwep 0.000001 --vcf 1kcp.vamos.tr.merge.modified.vcf  --out 1kcp.vamos.tr.merge.trtools
bcftools annotate -x FORMAT/FILTER 1kcp.vamos.tr.merge.trtools.vcf | bcftools view -f PASS -o 1kcp.vamos.tr.merge.trtools.filter.vcf
python3 vamos_tr_dosage.py 1kcp.vamos.tr.merge.trtools.filter.vcf 1kcp.vamos.tr.merge.trtools.filter.length.vcf 1kcp.vamos.tr.merge.trtools.filter.motif.vcf
dumpSTR --min-locus-het 0.02 --vcf 1kcp.vamos.tr.merge.trtools.filter.length.vcf --out 1kcp.vamos.tr.merge.trtools.filter.length.trtools
dumpSTR --min-locus-het 0.02 --vcf 1kcp.vamos.tr.merge.trtools.filter.motif.vcf --out 1kcp.vamos.tr.merge.trtools.filter.motif.trtools
bcftools annotate -x FORMAT/FILTER 1kcp.vamos.tr.merge.trtools.filter.length.trtools.vcf | bcftools view -f PASS | awk -v OFS='\t' '{if(substr($1,1,1)!="#"){split($8,a,";");split(a[length(a)-1],b,"=");$3=$1"-"$2"-TR-Length"};print $0}' | awk -v OFS='\t' '{if(substr($1,1,1)!="#") {for(i=10;i<=NF;i++){split($i,a,":"); $i="1|1:"a[2]":"a[3] }};print $0 }' | bgzip -c > 1kcp.tr_length.QC.vcf.gz 
tabix 1kcp.tr_length.QC.vcf.gz 
bcftools annotate -x FORMAT/FILTER 1kcp.vamos.tr.merge.trtools.filter.motif.trtools.vcf | bcftools view -f PASS | awk -v OFS='\t' '{if(substr($1,1,1)!="#"){split($8,a,";");split(a[length(a)-1],b,"=");$3=$1"-"$2"-TR-Motif-"b[2]};print $0}' | awk -v OFS='\t' '{if(substr($1,1,1)!="#") {for(i=10;i<=NF;i++){split($i,a,":"); $i="1|1:"a[2]":"a[3] }};print $0 }' | bgzip -c > 1kcp.tr_motif.QC.vcf.gz
tabix 1kcp.tr_motif.QC.vcf.gz 

###Merg genotype data
bcftools concat --threads $thread -o 1kcp.QC.vcf.gz -a 1kcp.small.QC.vcf.gz 1kcp.sv.QC.vcf.gz 1kcp.nest.vcf.gz 1kcp.tr_length.QC.vcf.gz 1kcp.tr_motif.QC.vcf.gz 
plink2 --vcf-half-call m --vcf 1kcp.QC.vcf.gz dosage=DS --out 1kcp --chr 1-22 --double-id
plink2 --keep $sample --make-pgen --out 1kcp.eqtl --pfile 1kcp

###Prepare gene expression data
Rscript extract_expression_matrix.R
Rscript gene_PEER.R
csvtk -t transpose 1kcp_Counts_TMM_RINT_AdjCov_60PEER.residuals.rint.txt | awk '{print $1}' | csvtk -H -t join -f "1;4" - $gencode | awk -v OFS='\t' '{print "chr"$2,$3-1,$3,$1}' > tmp1
csvtk -H -t join -f 1 $sample 1kcp_Counts_TMM_RINT_AdjCov_60PEER.residuals.rint.txt | csvtk -t transpose - | awk 'NR>3{print $0}' > tmp2
csvtk -H -t join -f 1 $sample 1kcp_Counts_TMM_RINT_AdjCov_60PEER.residuals.rint.txt | csvtk -t transpose - | awk 'NR==1{print "chr\tstart\tend\tgene_id\t"$0}' > tmp3
paste tmp1 tmp2 | bedtools sort | cat tmp3 - > 1kcp_Counts_TMM_RINT_AdjCov_60PEER.residuals.rint.bed

###eQTL mapping
python3 eQTL_mapping.py 1kcp.eqtl 1kcp.eqtl
