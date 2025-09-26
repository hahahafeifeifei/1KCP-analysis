#Rule: Annotate nodes for each paths
rule annotate_sample:
    input:
        ids_gfa="results/subgraph/subgraph_{id}/minigraph_{id}.seqwish.smoothxg.gfaffix.clip.kmer.chop.ids.gfa",
        anno_bed=config["anno"]
    output:
        node_bed=temp("results/subgraph/subgraph_{id}/anno/{sample}.{haplotype}.node.bed"),
        node_anno=temp("results/subgraph/subgraph_{id}/anno/{sample}.{haplotype}.node.anno")
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=60,
        runtime_hrs=10
    shell:
        """
            hap=$(echo {wildcards.haplotype} | awk '{{split($1,a,"hap");print a[2]-1}}')
            python3 scripts/graph_node_bed.py {input.ids_gfa} {output.node_bed} {wildcards.sample} $hap
            bedtools intersect -a {output.node_bed} -b {input.anno_bed} -wao | awk -v OFS='\\t' '{{if($9>=0.8*($3-$2)) print$4,$3-$2,$8}}' | sort | uniq > {output.node_anno}
        """

#Rule: Merging node annotations
rule merge_subgraph_annotate:
    input:
        node_annos=expand("results/subgraph/subgraph_{{id}}/anno/{sample}.{haplotype}.node.anno", sample=samples_list, haplotype=["hap1","hap2"]),
        ids_gfa="results/subgraph/subgraph_{id}/minigraph_{id}.seqwish.smoothxg.gfaffix.clip.kmer.chop.ids.gfa"
    output:
        node_count="results/subgraph/subgraph_{id}/minigraph_{id}.seqwish.smoothxg.gfaffix.clip.kmer.chop.ids.node.count",
        node_anno="results/subgraph/subgraph_{id}/minigraph_{id}.seqwish.smoothxg.gfaffix.clip.kmer.chop.ids.node.anno"
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=60,
        runtime_hrs=10
    shell:
        """
            grep -v "CHM13\|GRCh38\|MINIGRAPH" {input.ids_gfa} | grep ^W | \
            awk -v OFS='\\t' '{{split($7,a,"[><]");for(i=1;i<=length(a);i++) print $2,$3,a[i]}}' | \
            sort --buffer-size={resources.mem_gb}G | uniq | awk '{{print $3}}' | \
            sort --buffer-size={resources.mem_gb}G | uniq -c | awk -v OFS='\\t' '{{print $2,$1}}' > {output.node_count}

            if [ $(cat {input.node_annos} | wc -l ) -eq 0 ]
            then
                > {output.node_anno}
            else
                cat {input.node_annos} | sort --buffer-size={resources.mem_gb}G | uniq -c | awk -v OFS='\\t' '{{print $2,$3,$4,$1}}' |
                csvtk -H -t join -f 1 {output.node_count} - | awk -v OFS='\\t' '{{if($5/$2>0.8 && $5>1) print $1,$3,$4}}' > {output.node_anno}
            fi
        """

#Rule: Merge subgraph annotations
rule merge_annotate:
    input:
        annos=lambda wildcards: expand(
            "results/subgraph/subgraph_{id}/minigraph_{id}.seqwish.smoothxg.gfaffix.clip.kmer.chop.ids.node.anno",
            id=glob_wildcards(f"{checkpoints.minigraph_aln_partition.get(**wildcards).output[0]}/minigraph_{{id}}.gfa").id
        )
    output:
        anno="results/merge.seqwish.smoothxg.gfaffix.clip.kmer.chop.ids.node.anno"
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=60,
        runtime_hrs=10
    shell:
        """
            cat {input.annos} > {output.anno}
        """
