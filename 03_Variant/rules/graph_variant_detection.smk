#Rule: graph deconstruct
rule deconstruct:
    input:
        gfa=config['gfa'],
        samples=config['sample_list_file']
    output:
        deconstruct_vcf="results/subgraph/subgraph_{id}/subgraph_{id}.vcf.gz",
        deconstruct_vcf_tbi="results/subgraph/subgraph_{id}/subgraph_{id}.vcf.gz.tbi",
        vcfbub_vcf=temp("results/subgraph/subgraph_{id}/subgraph_{id}.vcfbub.vcf"),
        ref_site_vcf="results/subgraph/subgraph_{id}/subgraph_{id}.ref.site.vcf"
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=lambda wildcards, attempt: 100 * attempt,
        runtime_hrs=lambda wildcards, attempt: 20 * attempt
    threads: 8
    shell:
        """
            vg deconstruct -C -a -t {threads} -P GRCh38 {input.gfa} | bcftools view -s ^CHM13 | \
            bcftools norm -m -any | bcftools view -a | bcftools view -i "ALT!='.'" | bcftools norm -m +any | \
            python3 scripts/graph_vcf_norm.py - {input.samples} | sed "s/GRCh38.chr/chr/g" | bcftools view -o {output.deconstruct_vcf}
            tabix {output.deconstruct_vcf}
            vcfbub -a 100000 -l 0 --input {output.deconstruct_vcf} > {output.vcfbub_vcf}
            python3 scripts/graph_vcf_vt_annotate.py {output.vcfbub_vcf} {input.gfa} | \
            python3 scripts/graph_vcf_id_norm.py - > {output.ref_site_vcf}
        """

#Rule: small variant dection
rule small_varaint_call:
    input:
        gfa=config['gfa'],
        ref_site_vcf="results/subgraph/subgraph_{id}/subgraph_{id}.ref.site.vcf"
    output:
        decompose_bi_vcf=temp("results/subgraph/subgraph_{id}/subgraph_{id}.ref.decompose.biallelic.vcf"),
        decompose_multi_vcf=temp("results/subgraph/subgraph_{id}/subgraph_{id}.ref.decompose.vcf"),
        ref_small_vcf="results/subgraph/subgraph_{id}/subgraph_{id}.ref.small.allele.vcf.gz",
        ref_small_vcf_tbi="results/subgraph/subgraph_{id}/subgraph_{id}.ref.small.allele.vcf.gz.tbi"
    params:
        decompose_prefix="results/subgraph/subgraph_{id}/subgraph_{id}.ref.decompose"
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=60,
        runtime_hrs=10
    threads: 4
    shell:
        """
            python3 scripts/annotate_vcf_modified.py -vcf {input.ref_site_vcf} -gfa {input.gfa} -o {params.decompose_prefix}
            awk -v FS='\\t' '{{if(substr($0,1,1)=="#")print $0; else {{if($4!="" && ((length($5)-length($4))<50 && (length($4)-length($5))<50)) print $0}}}}' {output.decompose_bi_vcf} | \
            python3 scripts/graph_vcf_id_add_biallelic.py - | \
            awk -v OFS='\\t' '{{if(substr($0,1,1)!="#") {{if(length($5)==length($4)) $3=$3"-SNV"; else $3=$3"-INDEL"}}; print $0}}' | \
            bcftools sort -o {output.ref_small_vcf}
            tabix {output.ref_small_vcf}
        """

#Rule: SV dection
rule sv_call:
    input:
        reference=config['grch38'],
        gfa=config['gfa'],
        ref_site_vcf="results/subgraph/subgraph_{id}/subgraph_{id}.ref.site.vcf"
    output:
        ref_sv_raw_vcf="results/subgraph/subgraph_{id}/subgraph_{id}.ref.sv.raw.vcf.gz",
        ref_sv_raw_vcf_tbi="results/subgraph/subgraph_{id}/subgraph_{id}.ref.sv.raw.vcf.gz.tbi",
        ref_sv_raw_merge_vcf=temp("results/subgraph/subgraph_{id}/subgraph_{id}.ref.sv.merge.raw.vcf.gz"),
        ref_sv_site_merge_vcf=temp("results/subgraph/subgraph_{id}/subgraph_{id}.ref.sv.site.merge.vcf"),
        ref_sv_site_decompose_bi_vcf=temp("results/subgraph/subgraph_{id}/subgraph_{id}.ref.sv.site.decompose.biallelic.vcf"),
        ref_sv_site_decompose_multi_vcf=temp("results/subgraph/subgraph_{id}/subgraph_{id}.ref.sv.site.decompose.vcf"),
        ref_sv_vcf="results/subgraph/subgraph_{id}/subgraph_{id}.ref.sv.allele.vcf.gz",
        ref_sv_vcf_tbi="results/subgraph/subgraph_{id}/subgraph_{id}.ref.sv.allele.vcf.gz.tbi"
    params:
        decompose_prefix="results/subgraph/subgraph_{id}/subgraph_{id}.ref.sv.site.decompose"
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=60,
        runtime_hrs=10
    threads: 4
    shell:
        """
            bcftools view -i "(VT=='Biallelic_SV' || VT=='Multiallelic_SV')" {input.ref_site_vcf} | \
            bcftools norm -m -any | python3 scripts/graph_vcf_merge_split.py split | bcftools view -o {output.ref_sv_raw_vcf}
            tabix {output.ref_sv_raw_vcf}
            truvari collapse -i {output.ref_sv_raw_vcf} -f {input.reference} -k common --median-info -S 1000000 -s 0 -r 0 -p 0.7 -P 0.7 -o {output.ref_sv_raw_merge_vcf}
            bcftools sort {output.ref_sv_raw_merge_vcf} | sed "s/\\//\\|/g" | bcftools norm -m +any | python3 scripts/graph_vcf_merge_split.py merge > {output.ref_sv_site_merge_vcf}
            python3 scripts/annotate_vcf_modified.py -vcf {output.ref_sv_site_merge_vcf} -gfa {input.gfa} -o {params.decompose_prefix}
            awk -v FS='\\t' '{{if(substr($0,1,1)=="#")print $0; else {{if($4!="" && ((length($5)-length($4))>=50 || (length($4)-length($5))>=50)) print $0}}}}' {output.ref_sv_site_decompose_bi_vcf} | \
            python3 scripts/graph_vcf_id_add_biallelic.py - | awk -v OFS='\\t' '{{if(substr($0,1,1)!="#") $3=$3"-SV"; print $0}}' | \
            bcftools sort -o {output.ref_sv_vcf}
            tabix {output.ref_sv_vcf}
        """

#Rule: nested variant dection
rule nested_call:
    input:
        reference=config['grch38'],
        gfa=config['gfa'],
        ref_site_vcf="results/subgraph/subgraph_{id}/subgraph_{id}.ref.site.vcf",
        vcfbub_vcf="results/subgraph/subgraph_{id}/subgraph_{id}.vcfbub.vcf"
    output:
        snarls="results/subgraph/subgraph_{id}/subgraph_{id}.seqwish.smoothxg.gfaffix.clip.kmer.chop.ids.snarls",
        nest_raw_vcf="results/subgraph/subgraph_{id}/subgraph_{id}.nest.raw.vcf",
        nest_site_vcf="results/subgraph/subgraph_{id}/subgraph_{id}.nest.site.vcf.gz",
        nest_site_vcf_tbi="results/subgraph/subgraph_{id}/subgraph_{id}.nest.site.vcf.gz.tbi",
        nest_allele_vcf="results/subgraph/subgraph_{id}/subgraph_{id}.nest.allele.vcf.gz",
        nest_allele_vcf_tbi="results/subgraph/subgraph_{id}/subgraph_{id}.nest.allele.vcf.gz.tbi"
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=60,
        runtime_hrs=10
    threads: 4
    shell:
        """
            vg snarls {input.gfa} | vg view -R - > {output.snarls}
            python3 scripts/snarls_vcf_version5.py {input.vcfbub_vcf} {input.gfa} {output.snarls} {output.nest_raw_vcf}
            python3 scripts/graph_vcf_vt_annotate.py {output.nest_raw_vcf} {input.gfa} | \
            python3 scripts/graph_vcf_id_norm.py - | bcftools view -o {output.nest_site_vcf}
            tabix {output.nest_site_vcf}
            bcftools norm -m -any {output.nest_site_vcf} | python3 scripts/graph_vcf_id_add_biallelic.py - | \
            awk -v OFS='\\t' '{{if(substr($0,1,1)!="#")$3=$3"-NEST"; print $0}}' | bcftools view -o {output.nest_allele_vcf}
            tabix {output.nest_allele_vcf}
        """

#Rule: merge variant
rule merge_variant:
    input:
        small_vcfs=expand("results/subgraph/subgraph_{id}/subgraph_{id}.ref.small.allele.vcf.gz", id=id_list),
        sv_vcfs=expand("results/subgraph/subgraph_{id}/subgraph_{id}.ref.sv.allele.vcf.gz", id=id_list),
        nest_vcfs=expand("results/subgraph/subgraph_{id}/subgraph_{id}.nest.allele.vcf.gz", id=id_list)
    output:
        small_vcf="results/merge.ref.small.allele.vcf.gz",
        small_vcf_tbi="results/merge.ref.small.allele.vcf.gz.tbi",
        sv_vcf="results/merge.ref.sv.allele.vcf.gz",
        sv_vcf_tbi="results/merge.ref.sv.allele.vcf.gz.tbi",
        nest_vcf="results/merge.nest.allele.vcf.gz",
        nest_vcf_tbi="results/merge.nest.allele.vcf.gz.tbi"
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=100,
        runtime_hrs=10
    threads: 16
    shell:
        """
            bcftools concat --threads {threads} {input.small_vcfs} -a -o {output.small_vcf}
            tabix {output.small_vcf}
            bcftools concat --threads {threads} {input.sv_vcfs} -a -o {output.sv_vcf}
            tabix {output.sv_vcf}
            bcftools concat --threads {threads} {input.nest_vcfs} -a -o {output.nest_vcf}
            tabix {output.nest_vcf}
        """

#Filter variant
rule filter_variant:
    input:
        small_vcf="results/merge.ref.small.allele.vcf.gz",
        sv_vcf="results/merge.ref.sv.allele.vcf.gz",
        nest_vcf="results/merge.nest.allele.vcf.gz",
        unreliable=config['unreliable_bed']
    output:
        small_filter_vcf="results/merge.ref.small.allele.reliable.vcf.gz",
        small_filter_vcf_tbi="results/merge.ref.small.allele.reliable.vcf.gz.tbi",
        small_filter_id="results/merge.ref.small.allele.filter.id",
        sv_filter_vcf="results/merge.ref.sv.allele.reliable.vcf.gz",
        sv_filter_vcf_tbi="results/merge.ref.sv.allele.reliable.vcf.gz.tbi",
        sv_filter_id="results/merge.ref.sv.allele.filter.id",
        nest_filter_vcf="results/merge.nest.allele.reliable.vcf.gz",
        nest_filter_vcf_tbi="results/merge.nest.allele.reliable.vcf.gz.tbi",
        nest_filter_id="results/merge.nest.allele.filter.id"
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=100,
        runtime_hrs=10
    threads: 16
    shell:
        """
            bcftools view {input.small_vcf} -H | awk -v OFS='\\t' '{{print $1,$2-1,$2-1+length($4),$3}}' | \
            bedtools intersect -a - -b {input.unreliable} -wao | csvtk -H -t summary -g 1,2,3,4 -f 8:sum | \
            awk '{{if($5/($3-$2)>0.2) print $4}}' > {output.small_filter_id}
            bcftools view -e "ID=@{output.small_filter_id}" {input.small_vcf} -o {output.small_filter_vcf}
            tabix {output.small_filter_vcf}

            bcftools view {input.sv_vcf} -H | awk -v OFS='\\t' '{{print $1,$2-1,$2-1+length($4),$3}}' | \
            bedtools intersect -a - -b {input.unreliable} -wao | csvtk -H -t summary -g 1,2,3,4 -f 8:sum | \
            awk '{{if($5/($3-$2)>0.2) print $4}}' > {output.sv_filter_id}
            bcftools view -e "ID=@{output.sv_filter_id}" {input.sv_vcf} -o {output.sv_filter_vcf}
            tabix {output.sv_filter_vcf}

            bcftools view {input.nest_vcf} -H | awk -v OFS='\\t' '{{print $1,$2-1,$2,$3}}' | \
            bedtools intersect -a - -b {input.unreliable} -wao | csvtk -H -t summary -g 1,2,3,4 -f 8:sum | \
            awk '{{if($5/($3-$2)>0.2) print $4}}' > {output.nest_filter_id}
            bcftools view -e "ID=@{output.nest_filter_id}" {input.nest_vcf} -o {output.nest_filter_vcf}
            tabix {output.nest_filter_vcf}
        """

#Annotate variant
rule annotate_variant:
    input:
        gfas=expand(config['gfa'], id=id_list),
        sv_filter_vcf="results/merge.ref.sv.allele.reliable.vcf.gz",
        nest_filter_vcf="results/merge.nest.allele.reliable.vcf.gz",
        node_anno=config['node_anno']
    output:
        node_len=temp("results/merge.node.len"),
        sv_repeat_anno="results/merge.ref.sv.repeat.anno",
        nest_content_anno=temp("results/merge.ref.nest.content.anno"),
        nest_flanking_anno="results/merge.ref.nest.flanking.anno"
    conda:
        "../envs/env.yml"
    resources:
        mem_gb=100,
        runtime_hrs=10
    threads: 16
    shell:
        """
            cat {input.gfas} | grep ^S | awk '{{print $2"\\t"length($3)}}' > {output.node_len}
            bcftools view -H {input.sv_filter_vcf} | \
            awk '{{split($5,a,",");max_len=length($4);max_i=0;for(i=1;i<=length(a);i++) {{if(length(a[i])>max_len) {{max_len=length(a[i]);max_i=i}}}};split($8,b,";");split(b[1],c,"=");split(c[2],d,",") ;print $3"\\t"d[max_i+1]}}'  | \
            python3 scripts/variant_longest_annotate.py - {input.node_anno} {output.node_len} > {output.sv_repeat_anno}
            python3 scripts/variant_biallelic_annotate.py {input.nest_filter_vcf} {input.node_anno} {output.nest_content_anno} {output.nest_flanking_anno}
        """