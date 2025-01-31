#!/bin/bash
###parameter setting
sample=$1
thread=16
reference=hg38.no_alt.fa
hap1_fa=$sample.hap1.fasta
hap2_fa=$sample.hap2.fasta
truth_dir=$sample.truth

###Generating dipcall callset
run-dipcall -t $thread $sample $reference $hap1_fa $hap2_fa > $sample.mak
make -j2 -f $sample.mak
bcftools view -f . $sample.dip.vcf.gz | bcftools norm -f $reference -m -any |  awk '{if(substr($0,1,1)=="#")print$0; else {if(length($5)-length($4)<50 && length($4)-length($5)<50 ) print$0}}' | awk -v OFS='\t' '{if(substr($0,1,1)=="#") print$0;else print $1,$2,$3,toupper($4),toupper($5),$6,$7,$8,$9,$10}' | grep -v "*" > $sample.dipcall.small.vcf
bgzip -f $sample.dipcall.small.vcf
tabix $sample.dipcall.small.vcf.gz
bcftools view -f . $sample.dip.vcf.gz | bcftools norm -f $reference -m -any |  awk '{if(substr($0,1,1)=="#")print$0; else {if(length($5)-length($4)>=50 || length($4)-length($5)>=50 ) print$0}}' | awk -v OFS='\t' '{if(substr($0,1,1)=="#") print$0;else print $1,$2,$3,toupper($4),toupper($5),$6,$7,$8,$9,$10}' | grep -v "*" > $sample.dipcall.sv.vcf
bgzip -f $sample.dipcall.sv.vcf
tabix $sample.dipcall.sv.vcf.gz

###Small variant benchmark
hap.py $truth_dir/$sample.truth.small.vcf.gz $sample.dipcall.small.vcf.gz -f $truth_dir/$sample.confident.bed -r $reference -o $sample.hap --engine=vcfeval --no-roc --threads=${thread}

###SV benchmark
truvari bench --includebed $truth_dir/$sample.confident.bed -f $reference -b $truth_dir/$sample.truth.sv.vcf.gz -c $sample.dipcall.sv.vcf.gz -o $sample.truvari --passonly --pick ac --sizemin 50
cd $sample.truvari
truvari refine -r candidate.refine.bed -f reference ./ -t $thread --recount -u -U

###Haplotype phasing benchmark
whatshap compare --ignore-sample-name --names truth,whatshap --only-snvs --tsv-pairwise $sample.phase.bench.tsv --switch-error-bed $sample.phase.bench.bed $truth_dir/$sample.truth.small.vcf.gz $sample.dipcall.small.vcf.gz