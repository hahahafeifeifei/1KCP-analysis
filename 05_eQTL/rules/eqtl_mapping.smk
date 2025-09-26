#Rule: tensorqtl genotype
rule tensorqtl_genotype:
    input:
        small=config['small'],
        sv=config['sv'],
        nest=config['nest'],
        tr=config['tr'],
        reference_tr=config['reference_tr']
    output:
        small_genotype=temp("results/merge.small.genotype.vcf.gz"),
        sv_genotype=temp("results/merge.sv.genotype.vcf.gz"),
        nest_genotype=temp("results/merge.nest.genotype.vcf.gz"),
        tr_norm=temp("results/merge.tr.norm.vcf"),
        tr_trtools=temp("results/merge.tr.trtools.vcf"),
        tr_filter=temp("results/merge.tr.trtools.filter.vcf"),
        tr_length=temp("results/merge.tr.length.vcf"),
        tr_motif=temp("results/merge.tr.motif.vcf"),
        tr_length_trtools=temp("results/merge.tr.length.trtools.vcf"),
        tr_motif_trtools=temp("results/merge.tr.motif.trtools.vcf"),
        tr_length_genotype=temp("results/merge.tr.length.genotype.vcf.gz"),
        tr_motif_genotype=temp("results/merge.tr.motif.genotype.vcf.gz"),
        genotype="results/merge.genotype.vcf.gz",
        genotype_tbi="results/merge.genotype.vcf.gz.tbi"
    params:
        tr_length_trtools_prefix="results/merge.tr.length.trtools",
        tr_motif_trtools_prefix="results/merge.tr.motif.trtools"
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=30
    threads: 8
    shell:
        """
            bcftools plugin fill-tags {input.small} --threads {threads} | bcftools view -i "AF>0.01 && AF<0.99 && F_MISSING<0.1 && HWE>1e-6" -o {output.small_genotype}
            tabix {output.small_genotype} 
            bcftools plugin fill-tags {input.sv} --threads {threads} | bcftools view -i "AF>0.01 && AF<0.99 && F_MISSING<0.1 && HWE>1e-6" -o {output.sv_genotype}
            tabix {output.sv_genotype}
            bcftools plugin fill-tags {input.nest} --threads {threads} | bcftools view -i "AF>0.01 && AF<0.99 && F_MISSING<0.1 && HWE>1e-6" -o {output.nest_genotype}
            tabix {output.nest_genotype}

            awk -v OFS='\\t' '{{if(substr($1,1,1)=="#"){{if(substr($1,1,2)=="#C")print "##INFO=<ID=VID,Number=1,Type=Integer,Description=\\"VNTR ID\\">"}} else {{$3=$1"-"$2"-TR";$8=$8";VID="NR}};print $0}}' {input.tr} | \
            sed "s/source=vamos_2.1.7/source=adVNTR/g" > {output.tr_norm}
            dumpSTR --min-locus-callrate 0.9 --min-locus-hwep 0.000001 --vcf {output.tr_norm} --out {params.tr_trtools_prefix}
            bcftools annotate -x FORMAT/FILTER {output.tr_trtools} | bcftools view -f PASS -o {output.tr_filter}
            python3 scripts/vamos_tr_dosage.py {output.tr_filter} {output.tr_length} {output.tr_motif}
            dumpSTR --min-locus-het 0.02 --vcf {output.tr_length} --out {tr_length_trtools_prefix}
            dumpSTR --min-locus-het 0.02 --vcf {output.tr_motif} --out {tr_motif_trtools_prefix}

            bcftools annotate -x FORMAT/FILTER {output.tr_length_trtools} | sed "s/ALTANNO,Number=A/ALTANNO,Number=R/g" | \
            bcftools view -f PASS | awk -v OFS='\\t' '{{if(substr($1,1,1)!="#"){split($8,a,";");split(a[length(a)-1],b,"=");$3=$1"-"$2"-TR-Length"}};print $0}}' | \
            awk -v OFS='\\t' '{{if(substr($1,1,1)!="#") {{for(i=10;i<=NF;i++){{split($i,a,":"); $i="1|1:"a[2]":"a[3] }}}};print $0 }}' | bgzip -c > {output.tr_length_genotype} 
            tabix {output.tr_length_genotype} 

            bcftools annotate -x FORMAT/FILTER {output.tr_motif_trtools} | sed "s/ALTANNO,Number=A/ALTANNO,Number=R/g" | \
            bcftools view -f PASS | awk -v OFS='\\t' '{{if(substr($1,1,1)!="#"){{split($8,a,";");split(a[length(a)-1],b,"=");$3=$1"-"$2"-TR-Motif-"b[2]}};print $0}}' | \
            awk -v OFS='\\t' '{{if(substr($1,1,1)!="#") {{for(i=10;i<=NF;i++){{split($i,a,":"); $i="1|1:"a[2]":"a[3] }}}};print $0 }}' | bgzip -c > {output.tr_motif_genotype} 
            tabix {output.tr_motif_genotype} 

            bcftools concat -a --threads {threads} {output.small_genotype} {output.sv_genotype} {output.nest_genotype} {output.hla_genotype} {output.tr_length_genotype} {output.tr_motif_genotype} -o {output.genotype}
            tabix {output.genotype}
        """

#Rule: tensorqtl plink2
rule plink2:
    input:
        genotype="results/merge.genotype.vcf.gz",
    output:
        pvar="results/merge.genotype.pvar",
    params:
        plink2_prefix="results/merge.genotype"
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=100,
        runtime_hrs=30
    threads: 16
    shell:
        """
            plink2 --vcf-half-call m --vcf {input.genotype} dosage=DS --out {params.plink2_prefix} --chr 1-22 --double-id
        """

#Rule: tensorqtl
rule tensorqtl:
    input:
        pvar="results/merge.genotype.pvar",
        merge_list="results/merge.gene.list",
        cov=contig['cov']
    output:
        sig="results/merge.eqtl.txt",
        cond="results/merge.cond.eqtl.txt",
        susie="results/merge.susie.eqtl.txt",
        top="results/merge.top.eqtl.txt"
    params:
        plink2_prefix="results/merge.genotype",
        prefix="results/merge"
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=100,
        runtime_hrs=30
    threads: 16
    shell:
        """
            python3 scripts/eQTL_mapping.py {params.plink2_prefix} {input.merge_list} {input.cov} {params.prefix}
        """