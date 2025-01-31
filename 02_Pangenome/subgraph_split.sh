#!/bin/bash
###parameter setting
thread=16
reference=CHM13v2.0.fasta
prefix=1kcp
assembly_list=assembly.list

###Breaking minigraph pangenome
vg snarls $prefix.gfa | vg view -R - > $prefix.gfa.snarls.txt
odgi build -t $thread -g $prefix.gfa -o $prefix.og
python3 border_node_select.py $prefix.gfa $prefix.gfa.snarls.txt $prefix.og $prefix.node_split.gfa $prefix.node_split.edge_clip.gfa $prefix.subgraph.info
vg convert -g $prefix.node_split.gfa -f -Q chr | awk '{if(substr($0,1,1)=="S" && $6!="SR:i:0")print$0"\tSN:Z:Other\tSO:i:0\tSR:i:1";else print$0}' | awk -v OFS='\t' '{if(substr($0,1,1)=="S")print$1,"s"$2,$3,$4,$5,$6;else {if(substr($0,1,1)=="L")print$1,"s"$2,$3,"s"$4,$5,$6;else print$0} }' > $prefix.node_split.rgfa

###Align contigs to minigraph pangenome with divided nodes
cactus-graphmap ./jobstore $prefix.seqfile $prefix.node_split.rgfa $prefix.paf --outputGAFDir graphmap --outputFasta $prefix.fa.gz --reference CHM13 --mapCores 4 --delFilter 10000000 --defaultPreemptable --maxNodes $thread --logFile graphmap.log
zcat $prefix.gaf.gz  | gaffilter - -r 5.0 -m 0.25 -q 5 -b 250000 -o 0 -i 0.5 > $prefix.filter.gaf

###Generate 20Mb subgraph graph
mkdir subgraph
vg chunk -C -x $prefix.node_split.edge_clip.gfa --prefix subgraph/$prefix.subgraph -O gfa

###Generate graph alignments and sequences for subgraph graph
cd subgraph
seqfile=../../$prefix.seqfile
n=$[$(cat ../$prefix.seqfile | wc -l)]
row=$[$(cat ../$prefix.subgraph.info | wc -l)-1]
for i in `seq 0 $row`
do
        mkdir subgraph${i}
        mv $prefix.subgraph_${i}.* subgraph${i}
        cd subgraph${i}
        mkdir fa
        > $prefix.subgraph_${i}.seqfile
        awk "NR>1 && NR<${n}{print\$0}" $seqfile | while read line
        do
                sample=$(echo $line | awk '{print$1}')
                fa=$(echo $line | awk '{print$2}')
                awk -v sample=${sample} '{if($1==sample)print$2}' $prefix.subgraph_${i}.bed > fa/$sample.subgraph${i}.bed
                samtools faidx -r fa/$sample.subgraph${i}.bed $fa | awk -v sample=${sample} '{if(substr($0,1,1)==">")print ">id="sample"|"substr($0,2,length($0));else print$0}' > fa/$sample.subgraph${i}.fasta
                echo -e ${sample}"\t"fa/$sample.subgraph${i}.fasta >> $prefix.subgraph_${i}.seqfile
        done
        awk '{if($1=="S")print">id=_MINIGRAPH_|s"$2"\n"$3}' $prefix.subgraph_${i}.gfa > fa/_MINIGRAPH_.subgraph${i}.fasta
        echo -e _MINIGRAPH_"\t"fa/_MINIGRAPH_.subgraph${i}.fasta >> $prefix.subgraph_${i}.seqfile
        cd ..
done