#!/bin/bash
###parameter setting
thread=16
subgraph=$1
samples=sample.list

###SV repeat content anntoation
bcftools view -H $prefix.subgraph_${subgraph}.ref.sv.site.vcf.gz | awk '{split($5,a,",");max_len=length($4);max_i=0;for(i=1;i<=length(a);i++) {if(length(a[i])>max_len) {max_len=length(a[i]);max_i=i}};split($8,b,";");split(b[4],c,"=");split(c[2],d,",") ;print $3"\t"d[max_i+1]}' | python3 graph_vcf_repeat_annotate.py - $prefix.subgraph_${subgraph}.min0.8.anno $prefix.subgraph_${subgraph}.hap.count > $prefix.subgraph_${subgraph}.sv.site.anno

###Nested variant annotation
python3 graph_vcf_biallelic_annotate.py $prefix.subgraph_${subgraph}.nest.allele.vcf $prefix.subgraph_${subgraph}.min0.8.anno $prefix.subgraph_${subgraph}.variant.content.anno $prefix.subgraph_${subgraph}.variant.flanking.anno