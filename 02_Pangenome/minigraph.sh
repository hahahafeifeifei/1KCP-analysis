#!/bin/bash
###parameter setting
thread=16
reference=CHM13v2.0.fasta
prefix=1kcp
assembly_list=assembly.list

###Mash distance calculation
mash sketch $reference -o CHM13
> mash.dist
cat $assembly_list | while read assembly
do
    sample=$(ls $assembly | cut -d "." -f 1)
    hap=$(ls $assembly | cut -d "." -f 2)
    mash dist $assembly CHM13.msh | awk -v assembly=$assembly -v OFS='\t' '{print assembly,$3}' >> mash.dist
done

###Minigraph pangenome construction
fa=$(cat mash.dist | sort -k 2n | awk '{print $1}')
minigraph -c -x ggs -t $thread $reference $fa > $prefix.rgfa 
sed "s/s//" $prefix.rgfa | vg convert -fWg - > $prefix.gfa