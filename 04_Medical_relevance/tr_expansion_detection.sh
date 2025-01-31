#!/bin/bash
###parameter setting
thread=16

###TR site QC
awk -v OFS='\t' '{if(substr($1,1,1)=="#"){if(substr($1,1,2)=="#C")print "##INFO=<ID=VID,Number=1,Type=Integer,Description=\"VNTR ID\">"} else {$8=$8";VID="NR};print $0}' 1kcp.vamos.tr.merge.vcf  | sed "s/source=vamos_2.1.4/source=adVNTR/g" > 1kcp.vamos.tr.merge.modified.vcf 
dumpSTR --min-locus-callrate 0.8 --min-locus-hwep 0.000001 --vcf 1kcp.vamos.tr.merge.modified.vcf  --out 1kcp.vamos.tr.merge.trtools
bcftools annotate -x FORMAT/FILTER 1kcp.vamos.tr.merge.trtools.vcf | bcftools view -f PASS -o 1kcp.vamos.tr.merge.trtools.filter.vcf

###TR expansion threshold determination
python3 vamos_tr_outlier.py 1kcp.vamos.tr.merge.trtools.filter.vcf 1kcp.vamos.tr.merge.length_range.vcf