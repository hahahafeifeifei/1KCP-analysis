#!/bin/bash
###parameter setting
thread=16
gencode=gencode.v40.canonical.annotation.faa
assembly_list=assembly.list

###Individual protein alignment
cat $assembly_list | while read assembly
do
    prefix=$(ls $assembly | cut -d "." -f 1-2)
    sample=$(ls $assembly | cut -d "." -f 1)
    hap=$(ls $assembly | cut -d "." -f 2)
    miniprot --outs=0.97 --no-cs -Iut3 $assembly $gencode > $prefix.pangene.paf
done

###Pangene graph construction
pangene *pangene.paf -p 0 > 1kcp.pangene.gfa

###Gene cluster variation detection
k8 pangene.js call 1kcp.pangene.gfa > 1kcp.pangene.bubble