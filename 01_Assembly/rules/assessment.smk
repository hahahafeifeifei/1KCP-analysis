# Rule: Merqury assessment
rule merqury:
    input:
        sr_fq1 = config["sr_fq1"],
        sr_fq2 = config["sr_fq2"],
        fa_hap1='results/{sample}/{sample}.hap1.fasta',
        fa_hap2='results/{sample}/{sample}.hap2.fasta'
    output:
        merqury_dir = directory("results/{sample}/{sample}.merqury"),
        merqury_qv = "results/{sample}/{sample}.merqury/{sample}.merqury.qv",
        r1_meryl = temp(directory("results/{sample}/{sample}.R1.meryl")),
        r2_meryl = temp(directory("results/{sample}/{sample}.R2.meryl")),
        wgs_meryl = temp(directory("results/{sample}/{sample}.WGS.meryl"))
    params:
        merqury_prefix = "{sample}.merqury"
    conda:
        "../envs/env.yml"
    threads: 8
    resources:
        mem_gb=60,
        runtime_hrs=10
    shell:
        """
        meryl count k=21 memory={resources.mem_gb}G threads={threads} output {output.r1_meryl} {input.sr_fq1}
        meryl count k=21 memory={resources.mem_gb}G threads={threads} output {output.r2_meryl} {input.sr_fq2}
        meryl union-sum output {output.wgs_meryl} {output.r1_meryl} {output.r2_meryl}
        mkdir -p {output.merqury_dir}
        cd {output.merqury_dir}
        merqury.sh {output.wgs_meryl} {input.fa_hap1} {input.fa_hap2} {params.merqury_prefix}
        """

# Rule: Inspector assessment
rule inspector:
    input:
        hifi_fq = config["hifi_fq"],
        fa_hap1='results/{sample}/{sample}.hap1.fasta',
        fa_hap2='results/{sample}/{sample}.hap2.fasta'
    output:
        inspector_hap1_stat = "results/{sample}/{sample}.hap1.inspector/summary_statistics",
        inspector_hap2_stat = "results/{sample}/{sample}.hap2.inspector/summary_statistics"
    params:
        inspector_hap1_dir = "results/{sample}/{sample}.hap1.inspector",
        inspector_hap2_dir = "results/{sample}/{sample}.hap2.inspector"
    conda:
        "../envs/env.yml"
    threads: 8
    resources:
        mem_gb=60,
        runtime_hrs=10
    shell:
        """
        inspector.py -c {input.fa_hap1} -r {input.hifi_fq} --skip_base_error_detect --skip_base_error -t {threads} -o {params.inspector_hap1_dir}
        inspector.py -c {input.fa_hap2} -r {input.hifi_fq} --skip_base_error_detect --skip_base_error -t {threads} -o {params.inspector_hap2_dir}
        """