#!/bin/bash
###parameter setting
thread=16
subgraph=$1
reference=GRCh38_no_alt_analysis_set.fasta
samples=sample.list

###Reference variant dectection
vg deconstruct -C -t $thread -P GRCh38 $prefix.subgraph_${subgraph}.1kcp.gfa | python3 graph_vcf_norm.py - $samples | bgzip -c > $prefix.subgraph_${subgraph}.deconstruct.vcf.gz
tabix $prefix.subgraph_${subgraph}.deconstruct.vcf.gz
vcfbub -a 100000 -l 0 --input $prefix.subgraph_${subgraph}.deconstruct.vcf.gz > $prefix.subgraph_${subgraph}.vcfbub.vcf.gz
python graph_vcf_vt_annotate.py $prefix.subgraph_${subgraph}.vcfbub.vcf.gz $prefix.subgraph_${subgraph}.1kcp.gfa | python graph_vcf_id_norm.py - > $prefix.subgraph_${subgraph}.ref.site.vcf

###Reference small variant allele
python3 annotate_vcf_modified.py -vcf $prefix.subgraph_${subgraph}.ref.site.vcf -gfa $prefix.subgraph_${subgraph}.1kcp.gfa -o $prefix.subgraph_${subgraph}.ref.site.merge.decompose.vcf
awk -v FS='\t' '{if(substr($0,1,1)=="#")print $0;else {if($4!="" && ((length($5)-length($4))<50 && (length($4)-length($5))<50)) print $0}}' $prefix.subgraph_${subgraph}.ref.site.merge.decompose.vcf | python3 graph_vcf_id_add_biallelic.py - | awk -v OFS='\t' '{if(substr($0,1,1)!="#") {if(length($5)==length($4)) $3=$3"-SNV";else $3=$3"-INDEL"};print $0}' | bcftools sort -o $prefix.subgraph_${subgraph}.ref.small.allele.vcf 

###Reference SV allele
bcftools view -i "(VT=='Biallelic_SV' || VT=='Multiallelic_SV')" $prefix.subgraph_${subgraph}.ref.site.vcf | bcftools norm -m -any | python3 graph_vcf_merge_split.py split | bgzip -c > $prefix.subgraph_${subgraph}.ref.sv.site.vcf.gz
tabix $prefix.subgraph_${subgraph}.ref.sv.site.vcf.gz
truvari collapse -i $prefix.subgraph_${subgraph}.ref.sv.site.vcf.gz -f $reference -k common --median-info -S 1000000 -s 0 -r 0 -p 0.7 -P 0.7 -o $prefix.subgraph_${subgraph}.ref.sv.site.raw.merge.vcf.gz
bcftools sort $prefix.subgraph_${subgraph}.ref.sv.site.raw.merge.vcf.gz | sed "s/\//\|/g" | bcftools norm -m +any | python3 graph_vcf_merge_split.py merge > $prefix.subgraph_${subgraph}.ref.sv.site.merge.vcf
python3 annotate_vcf_modified.py -vcf $prefix.subgraph_${subgraph}.ref.sv.site.merge.vcf -gfa $prefix.subgraph_${subgraph}.1kcp.gfa -o $prefix.subgraph_${subgraph}.ref.sv.site.merge.decompose.vcf
awk -v FS='\t' '{if(substr($0,1,1)=="#")print $0;else {if($4!="" && ((length($5)-length($4))>=50 || (length($4)-length($5))>=50) ) print $0}}' $prefix.subgraph_${subgraph}.ref.sv.site.merge.decompose.vcf | python3 graph_vcf_id_add_biallelic.py - | awk -v OFS='\t' '{if(substr($0,1,1)!="#")$3=$3"-SV";print $0}' | bcftools sort -o $prefix.subgraph_${subgraph}.ref.sv.allele.vcf 

###Nested variant detection
vg snarls $prefix.subgraph_${subgraph}.1kcp.gfa | vg view -R - > $prefix.subgraph_${subgraph}.1kcp.snarls.txt
python3 snarls_vcf_version5.py $prefix.subgraph_${subgraph}.vcfbub.vcf.gz $prefix.subgraph_${subgraph}.1kcp.gfa $prefix.subgraph_${subgraph}.1kcp.snarls.txt $prefix.subgraph_${subgraph}.nest.raw.vcf
python graph_vcf_vt_annotate.py $prefix.subgraph_${subgraph}.nest.raw.vcf $prefix.subgraph_${subgraph}.1kcp.gfa | python graph_vcf_id_norm.py - > $prefix.subgraph_${subgraph}.nest.site.vcf

###Nested variant allele
bcftools norm -m -any $prefix.subgraph_${subgraph}.nest.site.vcf | python3 graph_vcf_id_add_biallelic.py - | awk -v OFS='\t' '{if(substr($0,1,1)!="#")$3=$3"-NEST";print $0}' | bcftools view -o $prefix.subgraph_${subgraph}.nest.allele.vcf
