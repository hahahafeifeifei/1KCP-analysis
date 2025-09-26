#Rule: Fasta rename
rule rename_reference:
    input:
        chm13_fa=config["chm13"],
        grch38_fa=config["grch38"]
    output:
        chm13_rename_fa='results/fasta/CHM13.fasta',
        chm13_rename_fai='results/fasta/CHM13.fasta.fai',
        grch38_rename_fa='results/fasta/GRCh38.fasta',
        grch38_rename_fai='results/fasta/GRCh38.fasta.fai'
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=20,
        runtime_hrs=10
    shell:
        """
        awk '{{if(substr($0,1,1)==">") print ">CHM13."substr($1,2,length($1));else print $0}}' {input.chm13_fa} > {output.chm13_rename_fa}
        awk '{{if(substr($0,1,1)==">") print ">GRCh38."substr($1,2,length($1));else print $0}}' {input.grch38_fa} > {output.grch38_rename_fa} 
        samtools faidx {output.chm13_rename_fa}
        samtools faidx {output.grch38_rename_fa}
        """

rule rename_assembly:
    input:
        fa='../01_Assembly/results/{sample}/{sample}.{haplotype}.fasta',
    output:
        rename_fa='results/fasta/{sample}.{haplotype}.fasta',
        rename_fai='results/fasta/{sample}.{haplotype}.fasta.fai'
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=20,
        runtime_hrs=10
    shell:
        """
        awk -v sample={wildcards.sample} -v hap={wildcards.haplotype} '{{if(substr($0,1,1)==">") print ">"sample"."hap"."substr($1,2,length($1));else print $0}}' {input.fa} > {output.rename_fa}
        samtools faidx {output.rename_fa}
        """

# Rule: mash for chm13 reference
rule mash_reference:
    input:
        chm13_fa='results/fasta/CHM13.fasta',
    output:
        chm13_msh='results/mash/chm13.msh'
    params:
        prefix="results/mash/chm13"
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=60,
        runtime_hrs=10
    shell:
        """
        mash sketch {input.chm13_fa} -o {params.prefix}
        """

# Rule: mash for assemblies
rule mash_assembly:
    input:
        fa='results/fasta/{sample}.{haplotype}.fasta',
        chm13_msh='results/mash/chm13.msh'
    output:
        dist=temp('results/mash/{sample}.{haplotype}.dist')
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=60,
        runtime_hrs=10
    shell:
        """
        mash dist {input.chm13_msh} {input.fa} | awk -v fa={input.fa} -v OFS='\\t' '{{print fa,$3}}' > {output.dist}
        """

#mash merge
rule mash_merge:
    input:
        dists=expand("results/mash/{sample}.{haplotype}.dist", sample=samples_list, haplotype=["hap1","hap2"]),
    output:
        final_dist="results/mash.dist"
    params:
        mash_dir=directory("results/mash")
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=60,
        runtime_hrs=10
    shell:
        """
        cat {input.dists} > {output.final_dist}
        rm -r {params.mash_dir}
        """

#minigraph construction
rule minigraph_construct:
    input:
        chm13_fa='results/fasta/CHM13.fasta',
        grch38_fa='results/fasta/GRCh38.fasta',
        dists="results/mash.dist"
    output:
        gfa="results/minigraph.gfa",
        snarls=temp("results/minigraph.snarls"),
        og=temp("results/minigraph.og"),
        node_clip_gfa="results/minigraph.node_clip.gfa",
        node_clip_rgfa="results/minigraph.node_clip.rgfa",
        node_edge_clip_gfa="results/minigraph.node_edge_clip.gfa",
        node_len=temp("results/minigraph.node_len.tsv"),
        subgraph_info="results/subgraph.info"
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=lambda wildcards, attempt: 200 * attempt,
        runtime_hrs=lambda wildcards, attempt: 24 * attempt
    threads: 16
    shell:
        """
        minigraph -c -x ggs -t {threads} {input.chm13_fa} {input.grch38_fa} $(cat {input.dists} | sort -k 2n | awk '{{print $1}}') | \
        sed "s/s//g" | vg convert -fWg - > {output.gfa}
        vg snarls {output.gfa} | vg view -R - > {output.snarls}
        odgi build -t {threads} -g {output.gfa} -o {output.og}
        export LD_PRELOAD=$CONDA_PREFIX/lib/libjemalloc.so.2
        python3 scripts/border_node_select.py {output.gfa} {output.snarls} {output.og} {output.node_clip_gfa} {output.node_edge_clip_gfa} {output.subgraph_info}
        vg convert -g {output.node_clip_gfa} -f -Q chr | awk '{{if(substr($0,1,1)=="S" && $6!="SR:i:0")print$0"\\tSN:Z:Other\\tSO:i:0\\tSR:i:1";else print$0}}' | \
        awk -v OFS='\\t' '{{if(substr($0,1,1)=="S")print$1,"s"$2,$3,$4,$5,$6;else {{if(substr($0,1,1)=="L")print$1,"s"$2,$3,"s"$4,$5,$6;else print$0}} }}' > {output.node_clip_rgfa}
        grep ^S {output.node_clip_rgfa} | awk '{{print $2"\\t"length($3)}}' > {output.node_len}
        """

#Rule: minigraph alignment
rule minigraph_alignment_reference:
    input:
        node_clip_rgfa="results/minigraph.node_clip.rgfa",
        fa='results/fasta/{sample}.{haplotype}.fasta'
    output:
        gaf="results/gaf/{sample}.{haplotype}.gaf"
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=100,
        runtime_hrs=30
    threads: 8
    shell:
        """
        minigraph -x asm -t {threads} --vc -c {input.node_clip_rgfa} {input.fa} | scripts/gaffilter - -r 5.0 -m 0.25 -q 5 -b 250000 -o 0 -i 0.5 > {output.gaf}
        """
        
rule minigraph_alignment_assembly:
    input:
        node_clip_rgfa="results/minigraph.node_clip.rgfa",
        chm13_fa='results/fasta/CHM13.fasta',
        grch38_fa='results/fasta/GRCh38.fasta'
    output:
        chm13_gaf="results/gaf/CHM13.gaf",
        grch38_gaf="results/gaf/GRCh38.gaf"
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=100,
        runtime_hrs=30
    threads: 8
    shell:
        """
        minigraph -x asm -t {threads} --vc -c {input.node_clip_rgfa} {input.chm13_fa} | scripts/gaffilter - -r 5.0 -m 0.25 -q 5 -b 250000 -o 0 -i 0.5 > {output.chm13_gaf}
        minigraph -x asm -t {threads} --vc -c {input.node_clip_rgfa} {input.grch38_fa} | scripts/gaffilter - -r 5.0 -m 0.25 -q 5 -b 250000 -o 0 -i 0.5 > {output.grch38_gaf}
        """

#Rule: minigraph alignment partition
checkpoint minigraph_aln_partition:
    input:
        gafs=expand("results/gaf/{sample}.{haplotype}.gaf", sample=samples_list, haplotype=["hap1","hap2"]),
        chm13_gaf="results/gaf/CHM13.gaf",
        grch38_gaf="results/gaf/GRCh38.gaf",
        node_len="results/minigraph.node_len.tsv",
        node_edge_clip_gfa="results/minigraph.node_edge_clip.gfa"
    output:
        subgraph_dir=directory("results/subgraph"),
        gaf="results/minigraph.gaf",
        paf="results/minigraph.paf"
    params:
        prefix="results/subgraph/minigraph"
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=100,
        runtime_hrs=30
    threads: 8
    shell:
        """
            mkdir {output.subgraph_dir}
            cat {input.gafs} {input.chm13_gaf} {input.grch38_gaf} > {output.gaf}
            scripts/gaf2paf -l {input.node_len} {output.gaf} > {output.paf}
            vg chunk -C -x {input.node_edge_clip_gfa} --prefix {params.prefix} -O gfa
            python3 scripts/subgraph_paf_fa.py {params.prefix} {output.paf} {output.gaf} {output.subgraph_dir}
        """

rule check_subgraph:
    input:
        lambda wildcards: expand(
            "results/subgraph/minigraph_{id}.gfa",
            id=glob_wildcards(f"{checkpoints.minigraph_aln_partition.get(**wildcards).output.subgraph_dir}/minigraph_{{id,[0-9]+}}.gfa").id
        )

#Rule: minigraph sequence partition
rule minigraph_seq_partition:
    input:
        bed="results/subgraph/minigraph_{id}.bed",
        gfa="results/subgraph/minigraph_{id}.gfa",
        paf="results/subgraph/minigraph_{id}.paf",
        fa_dir=directory("results/fasta")
    output:
        subgraph_dir=directory("results/subgraph/subgraph_{id}"),
        subgraph_fa_dir=directory("results/subgraph/subgraph_{id}/fasta"),
        fa="results/subgraph/subgraph_{id}/minigraph_{id}.fasta",
        bed="results/subgraph/subgraph_{id}/minigraph_{id}.bed",
        paf="results/subgraph/subgraph_{id}/minigraph_{id}.paf"
    params:
        prefix="result/subgraph/minigraph"
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=100,
        runtime_hrs=30
    threads: 1
    shell:
        """
            mkdir -p {output.subgraph_dir}
            mkdir -p {output.subgraph_fa_dir}
            awk '{{print $1}}' {input.bed} | while read sample
            do
                awk -v sample=$sample '{{if($1==sample)print$2}}' {input.bed} > {output.subgraph_fa_dir}/$sample.subgraph_{wildcards.id}.bed
                samtools faidx -r {output.subgraph_fa_dir}/$sample.subgraph_{wildcards.id}.bed {input.fa_dir}/$sample.fasta > {output.subgraph_fa_dir}/$sample.subgraph_{wildcards.id}.fasta
            done
            awk '{{if($1=="S")print">_MINIGRAPH_.s"$2"\\n"$3}}' {input.gfa} > {output.subgraph_fa_dir}/_MINIGRAPH_.subgraph_{wildcards.id}.fasta
            cat {output.subgraph_fa_dir}/*fasta > {output.fa}
            mv {input.bed} {output.bed}
            awk -v OFS='\\t' '{{$6="_MINIGRAPH_."$6; print $0}}' {input.paf} > {output.paf}
            rm {input.paf}
        """

#Rule: Seqwish graph inducing
rule seqwish:
    input:
        fa="results/subgraph/subgraph_{id}/minigraph_{id}.fasta",
        paf="results/subgraph/subgraph_{id}/minigraph_{id}.paf"
    output:
        raw_gfa="results/subgraph/subgraph_{id}/minigraph_{id}.seqwish.origin.gfa",
        seqwish_gfa="results/subgraph/subgraph_{id}/minigraph_{id}.seqwish.gfa"
    params:
        subgraph_dir="results/subgraph/subgraph_{id}"
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=lambda wildcards, attempt: 100 * attempt,
        runtime_hrs=lambda wildcards, attempt: 20 * attempt
    threads: 8
    shell:
        """
            seqwish -P -t {threads} -s {input.fa} -p {input.paf} -g {output.raw_gfa} -b {params.subgraph_dir}
            awk -v OFS='\\t' '{{if(substr($1,1,1)=="P") {{split($2,a,".");if(a[1]=="CHM13" || a[1]=="GRCh38" || a[1]=="_MINIGRAPH_") name=a[1]"#"a[1]"."a[2];else {{split(a[2],b,"hap");name=a[1]"#"b[2]-1"#"a[1]".hap"b[2]"."a[3]"#0"}}  print $1,name,$3,$4 }}else  print$0 }}' {output.raw_gfa} > {output.seqwish_gfa}
        """

#Rule: smoothxg graph smoothing
rule smoothxg:
    input:
        fa="results/subgraph/subgraph_{id}/minigraph_{id}.fasta",
        seqwish_gfa="results/subgraph/subgraph_{id}/minigraph_{id}.seqwish.gfa"
    output:
        smoothxg_gfa="results/subgraph/subgraph_{id}/minigraph_{id}.seqwish.smoothxg.gfa"
    params:
        smoothxg_dir="results/subgraph/subgraph_{id}/smoothxg_tmp"
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=lambda wildcards, attempt: 100 * attempt,
        runtime_hrs=lambda wildcards, attempt: 20 * attempt
    threads: 8
    shell:
        """
            mkdir -p {params.smoothxg_dir}
            sample_number=$(grep ">" {input.fa} | awk '{{split($1,a,".");if(length(a)==2) print a[1];else print a[1]"."a[2] }}' | sort -u | wc -l)
            smoothxg -t {threads} -g {input.seqwish_gfa} -r $sample_number --base {params.smoothxg_dir} --chop-to 100 -I 0.98 -R 0 -j 0 -e 0 -l 1400,1800,2200 -p 1,19,39,3,81,1 -O 0.001 -Y $[100*$sample_number] -d 0 -D 0 -V -c 200M -W 1 -o {output.smoothxg_gfa}
        """

#Rule: gfaffix norm
rule gfaffix:
    input:
        smoothxg_gfa="results/subgraph/subgraph_{id}/minigraph_{id}.seqwish.smoothxg.gfa"
    output:
        gfaffix_gfa="results/subgraph/subgraph_{id}/minigraph_{id}.seqwish.smoothxg.gfaffix.gfa"
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=lambda wildcards, attempt: 100 * attempt,
        runtime_hrs=lambda wildcards, attempt: 20 * attempt
    threads: 8
    shell:
        """
            gfaffix {input.smoothxg_gfa} | vg convert -fg - | \
            awk -v OFS='\\t' '{{if(substr($0,1,1)=="W" && $2!="_MINIGRAPH_") {{split($4,a,":");split(a[2],b,"-");print$1,$2,$3,a[1],$5+b[1]-1,$6+b[1]-1,$7}} else print$0 }}' > {output.gfaffix_gfa}
        """

