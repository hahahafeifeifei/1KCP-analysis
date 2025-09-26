# Rule: BRNN
rule brnn_reference:
    input:
        chm13_fa='results/fasta/CHM13.fasta',
        grch38_fa='results/fasta/GRCh38.fasta',
        brnn_knm=config['brnn_knm']
    output:
        chm13_brnn='results/brnn/CHM13.brnn.bed',
        grch38_brnn='results/brnn/GRCh38.brnn.bed'
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=60,
        runtime_hrs=10
    threads: 8
    shell:
        """
            scripts/dna-brnn {input.chm13_fa} -A -i {input.brnn_knm} -t {threads} | awk '{{if($3-$2 > 100000) print $0}}' | \
            bedtools sort -i - | bedtools merge -i - -d 100000 > {output.chm13_brnn}
            scripts/dna-brnn {input.grch38_fa} -A -i {input.brnn_knm} -t {threads} | awk '{{if($3-$2 > 100000) print $0}}' | \
            bedtools sort -i - | bedtools merge -i - -d 100000 > {output.grch38_brnn}
        """

rule brnn_assembly:
    input:
        fa='results/fasta/{sample}.{haplotype}.fasta',
        brnn_knm=config['brnn_knm']
    output:
        brnn='results/brnn/{sample}.{haplotype}.brnn.bed'
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=60,
        runtime_hrs=10
    threads: 8
    shell:
        """
        scripts/dna-brnn {input.fa} -A -i {input.brnn_knm} -t {threads} | awk '{{if($3-$2 > 100000) print $0}}' | \
        bedtools sort -i - | bedtools merge -i - -d 100000 > {output.brnn}
        """

#Rule: brnn merge
rule brnn_merge:
    input:
        brnns=expand("results/brnn/{sample}.{haplotype}.brnn.bed", sample=samples_list, haplotype=["hap1","hap2"]),
        chm13_brnn='results/brnn/CHM13.brnn.bed',
        grch38_brnn='results/brnn/GRCh38.brnn.bed'
    output:
        final_brnn="results/brnn.bed"
    params:
        brnn_dir=directory("results/brnn")
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=60,
        runtime_hrs=10
    shell:
        """
        cat {input.brnns} {input.chm13_brnn} {input.grch38_brnn} > {output.final_brnn}
        rm -r {params.brnn_dir}
        """

#Rule: region filtering
rule region_filter:
    input:
        final_brnn="results/brnn.bed",
        gfaffix_gfa="results/subgraph/subgraph_{id}/minigraph_{id}.seqwish.smoothxg.gfaffix.gfa"
    output:
        minigraph_filter_bed="results/subgraph/subgraph_{id}/minigraph_{id}.minigraph_filter.bed",
        merge_filter_bed="results/subgraph/subgraph_{id}/minigraph_{id}.merge_filter.bed",
        clip_gfa="results/subgraph/subgraph_{id}/minigraph_{id}.seqwish.smoothxg.gfaffix.clip.gfa"
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=100,
        runtime_hrs=30
    threads: 8
    shell:
        """
            python3 scripts/gfa_nominigraph_bed.py {input.gfaffix_gfa} {output.minigraph_filter_bed}
            cat {input.final_brnn} {output.minigraph_filter_bed} | bedtools sort | bedtools merge -d 1000 > {output.merge_filter_bed}
            python3 scripts/gfa_bed_filter.py {input.gfaffix_gfa} {output.merge_filter_bed} CHM13,GRCh38,_MINIGRAPH_ {output.clip_gfa}
        """

#Rule: kmc generating sample kmer
rule sample_kmer:
    input:
        fq1=config['sr_fq1'],
        fq2=config['sr_fq2']
    output:
        file='results/kmc/{sample}.read',
        histo=temp('results/kmc/{sample}.kmer.histo'),
        scope=temp(directory('results/kmc/{sample}.genomescope')),
        conf='results/kmc/{sample}.kmer.conf',
        kmc_list=temp('results/kmc/{sample}.kmer.list')
    params:
        kmc_tmp='results/kmc/{sample}.tmp',
        prefix='results/kmc/{sample}.kmer'
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=60,
        runtime_hrs=10
    threads: 4
    shell:
        """
            echo {input.fq1} > {output.file}
            echo {input.fq2} >> {output.file}
            mkdir -p {params.kmc_tmp}
            kmc -k29 -m20 -t{threads} -hp -ci1 -cs10000 @{output.file} {params.prefix} {params.kmc_tmp}
            kmc_tools transform {params.prefix} histogram {output.histo} -cx10000
            genomescope2 -i {output.histo} -o {output.scope} -k 29
            python3 scripts/confident_range.py {output.scope}/model.txt > {output.conf}
            echo -e {wildcards.sample}"\\t"{params.prefix} > {output.kmc_list}
        """

#Rule: merge kmer file
rule merge_kmer:
    input:
        kmers=expand("results/kmc/{sample}.kmer.list", sample=samples_list)
    output:
        kmc_list='results/kmer.list'
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=60,
        runtime_hrs=10
    shell:
        """
            cat {input.kmers} > {output.kmc_list}
        """

#Rule: kmer filtering
rule kmer_filter:
    input:
        clip_gfa="results/subgraph/subgraph_{id}/minigraph_{id}.seqwish.smoothxg.gfaffix.clip.gfa",
        kmc_list='results/kmer.list'
    output:
        kmer_supp_node=temp("results/subgraph/subgraph_{id}/kmer_supp.node"),
        kmer_supp_edge=temp("results/subgraph/subgraph_{id}/kmer_supp.edge"),
        kmer_filter_node="results/subgraph/subgraph_{id}/kmer_filter.node",
        kmer_filter_edge="results/subgraph/subgraph_{id}/kmer_filter.edge",
        kmer_filter_gfa="results/subgraph/subgraph_{id}/minigraph_{id}.seqwish.smoothxg.gfaffix.clip.kmer.gfa",
        chop_gfa="results/subgraph/subgraph_{id}/minigraph_{id}.seqwish.smoothxg.gfaffix.clip.kmer.chop.gfa",
        chop_node_count="results/subgraph/subgraph_{id}/chop.node.count"
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=100,
        runtime_hrs=10
    threads: 8
    shell:
        """
            python3 scripts/graph_kmer_mismatch_filter.py {input.clip_gfa} {input.kmc_list} CHM13,GRCh38 {output.kmer_supp_node} {output.kmer_supp_edge} 0.2
            awk '{{if($2<0.5) print $1}}' {output.kmer_supp_node} > {output.kmer_filter_node}
            awk '{{if($2<0.5) print $1}}' {output.kmer_supp_edge} > {output.kmer_filter_edge}
            python3 scripts/gfa_node_edge_filtering.py {input.clip_gfa} {output.kmer_filter_node} {output.kmer_filter_edge} GRCh38,CHM13 {output.kmer_filter_gfa}
            grep -v "_MINIGRAPH_" {output.kmer_filter_gfa} | vg clip -d 1 -P GRCh38 -P CHM13 - | vg mod -u - | vg mod -X 20 - | vg view - > {output.chop_gfa}
            grep ^S {output.chop_gfa} | wc -l | awk -v id={wildcards.id} '{{print id"\\t"$1}}' > {output.chop_node_count}
        """

#Rule: subgraph node count
rule node_count_merge:
    input:
        counts=lambda wildcards: expand(
            "results/subgraph/subgraph_{id}/chop.node.count",
            id=glob_wildcards(f"{checkpoints.minigraph_aln_partition.get(**wildcards).output[0]}/minigraph_{{id}}.gfa").id
        )
    output:
        node_counts=temp("results/chop.node.count"),
        node_sum_counts="results/chop.node.sum.count"
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=60,
        runtime_hrs=10
    shell:
        """
            cat {input.counts} > {output.node_counts}
            awk 'BEGIN{{sum=0}}{{print $1"\\t"sum;sum+=$2}}' {output.node_counts} > {output.node_sum_counts}
        """

#Rule: subgraph ids
rule node_ids:
    input:
        chop_gfa="results/subgraph/subgraph_{id}/minigraph_{id}.seqwish.smoothxg.gfaffix.clip.kmer.chop.gfa",
        node_sum_counts="results/chop.node.sum.count"
    output:
        ids_gfa="results/subgraph/subgraph_{id}/minigraph_{id}.seqwish.smoothxg.gfaffix.clip.kmer.chop.ids.gfa",
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=60,
        runtime_hrs=10
    shell:
        """
            index=$(awk -v id={wildcards.id} '{{if($1==id) print $2}}' {input.node_sum_counts})
            python3 scripts/gfa_ids.py {input.chop_gfa} {output.ids_gfa} $index
        """

#Rule: extract subgraph with reference path
rule subgraph_ref:
    input:
        ids_gfa="results/subgraph/subgraph_{id}/minigraph_{id}.seqwish.smoothxg.gfaffix.clip.kmer.chop.ids.gfa",
    output:
        ref=temp("results/subgraph/subgraph_{id}/gfa_grch38.list")
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=60,
        runtime_hrs=10
    shell:
        """
            awk '{{if($2=="GRCh38") print $0}}' {input.ids_gfa} | wc -l | awk -v id={wildcards.id} '{{if($1>0) print id}}' > {output.ref}
        """

rule merge_ref:
    input:
        refs=lambda wildcards: expand(
            "results/subgraph/subgraph_{id}/gfa_grch38.list",
            id=glob_wildcards(f"{checkpoints.minigraph_aln_partition.get(**wildcards).output[0]}/minigraph_{{id}}.gfa").id
        )
    output:
        ref_list="results/gfa_grch38.list"
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=60,
        runtime_hrs=10
    shell:
        """
            cat {input.refs} > {output.ref_list}
        """