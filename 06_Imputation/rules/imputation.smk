#Rule: build imputation panel
rule panel_build:
    input:
        small=config['small'],
        sv=config['sv'],
        nest=config['nest'],
        tr=config['tr'],
        hla=config['hla'],
        reference_tr=config['reference_tr']
    output:
        small_panel=temp("results/merge.small.panel.vcf.gz"),
        sv_panel=temp("results/merge.sv.panel.vcf.gz"),
        nest_panel=temp("results/merge.nest.panel.vcf.gz"),
        hla_panel=temp("results/merge.hla.panel.vcf.gz"),
        tr_norm=temp("results/merge.tr.norm.vcf"),
        tr_trtools=temp("results/merge.tr.trtools.vcf"),
        tr_filter=temp("results/merge.tr.trtools.filter.vcf"),
        tr_addref=temp("results/merge.tr.addref.vcf"),
        tr_addref_trtools=temp("results/merge.tr.addref.trtools.vcf"),
        tr_panel=temp("results/merge.tr.panel.vcf.gz"),
        panel="results/merge.panel.vcf.gz",
        panel_tbi="results/merge.panel.vcf.gz.tbi"
    params:
        tr_trtools_prefix="results/merge.tr.trtools",
        tr_addref_trtools_prefix="results/merge.tr.addref.trtools"
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=30
    threads: 8
    shell:
        """
            bcftools plugin fill-tags {input.small} --threads {threads} | bcftools view -i "AC>1 && AN-AC>1 && F_MISSING<0.1 && HWE>1e-6" -o {output.small_panel}
            tabix {output.small_panel} 
            bcftools plugin fill-tags {input.sv} --threads {threads} | bcftools view -i "AC>1 && AN-AC>1 && F_MISSING<0.1 && HWE>1e-6" -o {output.sv_panel}
            tabix {output.sv_panel}
            bcftools plugin fill-tags {input.nest} --threads {threads} | bcftools view -i "AC>1 && AN-AC>1 && F_MISSING<0.1 && HWE>1e-6" -o {output.nest_panel}
            tabix {output.nest_panel}
            bcftools plugin fill-tags {input.hla} --threads {threads} | bcftools view -i "AC>1 && AN-AC>1 && F_MISSING<0.1 && HWE>1e-6" -o {output.hla_panel}
            tabix {output.hla_panel}

            awk -v OFS='\\t' '{{if(substr($1,1,1)=="#"){{if(substr($1,1,2)=="#C")print "##INFO=<ID=VID,Number=1,Type=Integer,Description=\\"VNTR ID\\">"}} else {{$3=$1"-"$2"-TR";$8=$8";VID="NR}};print $0}}' {input.tr} | \
            sed "s/source=vamos_2.1.7/source=adVNTR/g" > {output.tr_norm}
            dumpSTR --min-locus-callrate 0.9 --min-locus-hwep 0.000001 --vcf {output.tr_norm} --out {params.tr_trtools_prefix}
            bcftools annotate -x FORMAT/FILTER {output.tr_trtools} | bcftools view -f PASS -o {output.tr_filter}
            python3 scripts/tr_ref.py {output.tr_filter} {input.reference_tr} {output.tr_addref}
            dumpSTR --vcf {output.tr_addref} --out {params.tr_addref_trtools_prefix}
            bcftools view -i "AC>1 && REFAC>1" {output.tr_addref_trtools} | bcftools annotate -x FORMAT/FILTER | sed "s/ALTANNO,Number=A/ALTANNO,Number=R/g" | \
            bcftools norm -m -any | python3 scripts/graph_vcf_id_add_biallelic.py - | bcftools view -o {output.tr_panel}
            tabix {output.tr_panel}

            bcftools concat -a --threads {threads} {output.small_panel} {output.sv_panel} {output.nest_panel} {output.hla_panel} {output.tr_panel} -o {output.panel}
            tabix {output.panel}
        """

#Rule: perform leave-one-out imputation
rule leave_one_imputation:
    input:
        panel="results/merge.panel.vcf.gz",
        omni=config['omni']
    output:
        panel_leave_one=temp("results/{sample}/{sample}_out.{chr}.vcf.gz"),
        msav_leave_one=temp("results/{sample}/{sample}_out.{chr}.msav"),
        sample_omni=temp("results/{sample}/{sample}.{chr}.array.vcf.gz"),
        sample_omni_tbi=temp("results/{sample}/{sample}.{chr}.array.vcf.gz.tbi"),
        sampel_impute="results/{sample}/{sample}.{chr}.impute.vcf.gz"
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=30
    threads: 2
    shell:
        """
            bcftools view --threads {threads} -s ^{wildcards.sample} {input.panel} {wildcards.chr} -o {output.panel_leave_one}
            minimac4 -t {threads} --compress-reference {output.panel_leave_one} -o {output.msav_leave_one}
            bcftools view -v snps --threads {threads} -s {wildcards.sample} -R {input.omni} {input.panel} {wildcards.chr} -o {output.sample_omni}
            tabix {output.sample_omni}
            minimac4 --format DS,GT -t {threads} {output.msav_leave_one} {output.sample_omni} -o {output.sampel_impute}
        """
