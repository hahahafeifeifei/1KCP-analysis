#!/bin/bash
###parameter setting
thread=16
subgraph=$1
prefix=1kcp
mask=${prefix}.dna-brnn.bed
node_index=${prefix}.node.index

cd subgraph/subgraph${subgraph}
###Base-level graph initialization
cat fa/*.fasta > $prefix.subgraph_${subgraph}.fasta
seqwish -P -t $thread -s $prefix.subgraph_${subgraph}.fasta -p $prefix.subgraph_${subgraph}.paf -g $prefix.subgraph_${subgraph}.seqwish.origin.gfa
awk -v OFS='\t' '{if(substr($1,1,1)=="P") {split($2,a,"id=");split(a[2],b,"|");if(b[1]=="CHM13" || b[1]=="GRCh38" || b[1]=="_MINIGRAPH_") name=b[1]"#"b[2];else {split(b[1],c,".");name=c[1]"#"c[2]-1"#"b[2]"#0"}  print $1,name,$3,$4 }else  print$0  }' $prefix.subgraph_${subgraph}.seqwish.origin.gfa > $prefix.subgraph_${subgraph}.seqwish.gfa

###Graph refinement
sample_number=$(ls -lh fa/*fasta | awk '{if($5!=0)print $0}' | wc -l)
mkdir smoothxg_tmp
smoothxg -t $thread -g $prefix.subgraph_${subgraph}.seqwish.gfa -r $sample_number --base smoothxg_tmp --chop-to 100 -I 0.98 -R 0 -j 0 -e 0 -l 1400,1800,2200 -q 2800 -w $[1400*${sample_number}] -p 1,19,39,3,81,1 -O 0.001 -Y $[100*${sample_number}] -d 0 -D 0 -V -c 200M -W 1 -o $prefix.subgraph_${subgraph}.seqwish.smoothxg.gfa

###Graph normalization
gfaffix $prefix.subgraph_${subgraph}.seqwish.smoothxg.gfa -o $prefix.subgraph_${subgraph}.seqwish.smoothxg.gfaffix.gfa > $prefix.subgraph_${subgraph}.seqwish.smoothxg.gfaffix.info
vg convert -fg $prefix.subgraph_${subgraph}.seqwish.smoothxg.gfaffix.gfa | awk -v OFS='\t' '{if(substr($0,1,1)=="W" && $2!="_MINIGRAPH_") {split($4,a,":");split(a[2],b,"-");print$1,$2,$3,a[1],$5+b[1]-1,$6+b[1]-1,$7} else print$0 }' > $prefix.subgraph_${subgraph}.seqwish.smoothxg.gfaffix.norm.gfa

###Graph filter
python3 gfa_nominigraph_bed.py $prefix.subgraph_${subgraph}.seqwish.smoothxg.gfaffix.norm.gfa $prefix.subgraph_${subgraph}.nominigraph.bed
cat $mask $prefix.subgraph_${subgraph}.nominigraph.bed | bedtools sort | bedtools merge -d 1000 > $prefix.subgraph_${subgraph}.merge.bed
python3 gfa_bed_filter.py $prefix.subgraph_${subgraph}.seqwish.smoothxg.gfaffix.norm.gfa $prefix.subgraph_${subgraph}.merge.bed CHM13,GRCh38,_MINIGRAPH_ $prefix.subgraph_${subgraph}.seqwish.smoothxg.gfaffix.norm.clip.gfa
grep -v "_MINIGRAPH_" $prefix.subgraph_${subgraph}.seqwish.smoothxg.gfaffix.norm.clip.gfa | vg clip -d 1 -P GRCh38 -P CHM13 - | vg mod -u - | vg view - > $prefix.subgraph_${subgraph}.full-sample.gfa

###Extract 1KCP pangenome
grep "^S\|^L\|GRCh38\|CHM13\|1KCP" $prefix.subgraph_${subgraph}.full-sample.gfa | grep -v "modest_cov" | vg clip -d 1 -P GRCh38 -P CHM13 - | vg mod -u - | vg mod -X 20 | vg view - > $prefix.subgraph_${subgraph}.1kcp.raw.gfa
python3 gfa_ids.py $prefix.subgraph_${subgraph}.1kcp.raw.gfa $prefix.subgraph_${subgraph}.1kcp.gfa $node_index