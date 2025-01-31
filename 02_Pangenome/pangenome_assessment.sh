#!/bin/bash
###parameter setting
thread=16
subgraph=$1
prefix=1kcp
assembly_anno=$prefix.assembly.anno.list
hprc=hprc.sample
cpc=cpc.sample
1kcp=1kcp.sample
benchmark=benchmark.sample

cd subgraph/subgraph${subgraph}
###Nonreference size and comparison with other assembly panel
grep ^S $prefix.subgraph_${subgraph}.full-sample.gfa | awk '{print $2"\t"length($3)}' > node.len
grep "CHM13\|GRCh38" $prefix.subgraph_${subgraph}.full-sample.gfa | awk '{split($7,a,"[<>]");for(i=2;i<=length(a);i++)print $2"\t"$3"\t"a[i]}' | sort | uniq | awk '{print $3}' | sort | uniq -c | awk '{print $2"\t"$1}' > ref.node
grep -f $1kcp $prefix.subgraph_${subgraph}.full-sample.gfa | awk '{split($7,a,"[<>]");for(i=2;i<=length(a);i++)print $2"\t"$3"\t"a[i]}' | sort | uniq | awk '{print $3}' | sort | uniq -c | awk '{print $2"\t"$1}' > 1kcp.node
grep -f $hprc $prefix.subgraph_${subgraph}.full-sample.gfa | awk '{split($7,a,"[<>]");for(i=2;i<=length(a);i++)print $2"\t"$3"\t"a[i]}' | sort | uniq | awk '{print $3}' | sort | uniq -c | awk '{print $2"\t"$1}' > hprc.node
grep -f $cpc $prefix.subgraph_${subgraph}.full-sample.gfa | awk '{split($7,a,"[<>]");for(i=2;i<=length(a);i++)print $2"\t"$3"\t"a[i]}' | sort | uniq | awk '{print $3}' | sort | uniq -c | awk '{print $2"\t"$1}' > cpc.node
csvtk -H -t join -L --na 0 -f 1 node.len ref.node | csvtk -H -t join -L --na 0 -f 1 - 1kcp.node | csvtk -H -t join -L --na 0 -f 1 - hprc.node | csvtk -H -t join -L --na 0 -f 1 - cpc.node > $prefix.subgraph_${subgraph}.full-sample.hap.count

###Compare the node and edge consistence between PIGA and Hifiasm assembly
> $prefix.subgraph_${subgraph}.graph.compare
cat $benchmark | while read sample
do
        python3 graph_compare.py $prefix.subgraph_${subgraph}.full-sample.gfa ${sample}_low ${sample} ${sample}_low.node ${sample}_low.edge
        n1=$(grep FP ${sample}_low.node | wc -l)
        n2=$(grep TP ${sample}_low.node | wc -l)
        n3=$(grep FP ${sample}_low.edge | wc -l)
        n4=$(grep TP ${sample}_low.edge | wc -l)
        python3 graph_compare.py $prefix.subgraph_${subgraph}.full-sample.gfa ${sample} ${sample}_low ${sample}.node ${sample}.edge
        n5=$(grep FP ${sample}.node | wc -l)
        n6=$(grep TP ${sample}.node | wc -l)
        n7=$(grep FP ${sample}.edge | wc -l)
        n8=$(grep TP ${sample}.edge | wc -l)
        echo -e ${sample}_low"\t"${n1}"\t"${n2}"\t"${n3}"\t"${n4} >> $prefix.subgraph_${subgraph}.graph.compare
        echo -e ${sample}"\t"${n5}"\t"${n6}"\t"${n7}"\t"${n8} >> $prefix.subgraph_${subgraph}.graph.compare
done

###Mapping haplotype coverage to the CHM13
grep ^S $prefix.subgraph_${subgraph}.1kcp.gfa | awk '{print $2"\t"length($3)}' > node.len
grep "CHM13\|GRCh38" $prefix.subgraph_${subgraph}.1kcp.gfa | awk '{split($7,a,"[<>]");for(i=2;i<=length(a);i++)print $2"\t"$3"\t"a[i]}' | sort | uniq | awk '{print $3}' | sort | uniq -c | awk '{print $2"\t"$1}' > ref.node
grep -f $1kcp $prefix.subgraph_${subgraph}.1kcp.gfa | awk '{split($7,a,"[<>]");for(i=2;i<=length(a);i++)print $2"\t"$3"\t"a[i]}' | sort | uniq | awk '{print $3}' | sort | uniq -c | awk '{print $2"\t"$1}' > 1kcp.node
csvtk -H -t join -L --na 0 -f 1 node.len ref.node | csvtk -H -t join -L --na 0 -f 1 - 1kcp.node > $prefix.subgraph_${subgraph}.hap.count
python3 gfa_node_bed.py $prefix.subgraph_${subgraph}.1kcp.gfa CHM13.node.bed CHM13 0
csvtk -H -t join -f "4;1" CHM13.node.bed $prefix.subgraph_${subgraph}.hap.count | awk '{print $2"\t"$3"\t"$4"\t"$6+$7}' > $prefix.subgraph_${subgraph}.1kcp.chm13_cov\

