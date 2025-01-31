#!/bin/bash
###parameter setting
thread=16
subgraph=$1
prefix=1kcp
assembly_anno=$prefix.assembly.anno.list

cd subgraph/subgraph${subgraph}
mkdir graph_annotate
cd graph_annotate
###Node annotation per assembly
cat $assembly_anno | while read line
do
    sample=$(echo $line | awk '{print $1}')
    id=$(echo $sample | cut -d "." -f 1)
    hap=$(echo $sample | cut -d "." -f 2)
    anno=$(echo $line | awk '{print $2}')
    python3 gfa_node_bed.py ../$prefix.subgraph_${subgraph}.1kcp.gfa $sample.node.bed $id $hap
    bedtools intersect -a $sample.node.bed -b $anno -wao | awk '{if($9>=0.8*($3-$2)) print $4"\t"$8}' | sort | uniq > $sample.node.anno
done

###Merge the node annotation across the population
cat *node.anno | sort --buffer-size=25G | uniq -c > ../$prefix.subgraph_${subgraph}.anno.count
cd ..

###Filter the singleton node annotation and annotation support by <80% of assemblies
cat $prefix.subgraph_${subgraph}.anno.count | sed 's/STR\|VNTR/TR/g' | awk '{gsub(/[0-9]+/, "",$3);print $0}' | sed 's/nc_intron\|pc_intron/intron/g' | grep -v pc_exon | grep -v TN | \
 awk '{print $2"\t"$3"\t"$1}' | csvtk -H -t join -f 1 $prefix.subgraph_${subgraph}.hap.count  - | awk -v OFS='\t' '{if($8/($4+$3)>0.8 && ($4+$3)>1)print $1,$2,$7,$3,$4,$8}' > $prefix.subgraph_${subgraph}.min0.8.anno
