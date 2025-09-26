# Rule: hla typing
rule hla_typing:
    input:
        fa=config['fa'],
        bam="results/bam/{sample}.{haplotype}.bam",
        hla_len=config['hla_len']
    output:
        hla_dir=directory("results/hla/{sample}.{haplotype}"),
        hla="results/hla/{sample}.{haplotype}.hla"
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=30
    threads: 8
    shell:
        """
            mkdir -p {output.hla_dir}
            hifihla call-contigs --abam {input.bam} --hap1 {input.fa} --out_prefix {output.hla_dir}
            csvtk -H -t join -f "4;1" {output.hla_dir}/hifihla_summary.tsv {input.hla_len} | \
            awk -v OFS='\\t' -v hap={wildcards.haplotype} -v sample={wildcards.sample} '{{if($4!="N/A") {{split($4,a,"*");div=100-$7/$NF*100; print sample,hap,a[1],$4,$6,div}}}}' | \
            csvtk -H -t sort -k 5:nr -k 6:nr | csvtk -H -t uniq -f 3  > {output.hla}
        """

#Rule: merge hla typing results
rule hla_merge:
    input:
        hlas=expand("results/hla/{sample}.{haplotype}.hla", sample=samples_list, haplotype=["hap1","hap2"]),
        hla_pos=config['hla_pos'],
        samples=config['sample_list_file']
    output:
        raw_hla_list=temp("results/merge.hla.list"),
        hla_list="results/merge.hla.reliable.list",
        hla_decompose_list=temp("results/merge.hla.reliable.decompose.list"),
        hla_decompose_geno="results/merge.hla.reliable.decompose.geno",
        hla_vcf="results/merge.hla.reliable.vcf"
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=200,
        runtime_hrs=30
    threads: 8
    shell:
        """
            cat {input.hlas} > {output.raw_hla_list}
            awk '{{if($5>=95 && $6>=99 && $6!="N/A") print $0}}' {output.raw_hla_list} | grep -v "HLA-U\|HLA-K\|HLA-W" > {output.hla_list}
            awk -v OFS='\\t' '{{print $1,$2,$3,$4,"4"}}' {output.hla_list} | \
            awk -v OFS='\\t' '{{for(i=1;i<=$5;i++) {{split($4,a,"*");split(a[2],b,":"); printf $1"\\t"$2"\\t"$3"\\t"a[1]"*"b[1];max=i;if(i>length(b)) max=length(b); for(j=2;j<=max;j++) printf ":"b[j]; print "\\t"i }}}}' > {output.hla_decompose_list}
            awk -v OFS='\\t' '{{print $1"."$2,$3,$4,$5}}' {output.hla_decompose_list} | csvtk -H -t fold -f 2,3,4 -v 1 -s ","| \
            awk -v OFS='\\t' '{{split($4,a,",");print $1,$2,$3,length(a),$4}}' | csvtk -H -t sort -k 1 -k 3:n -k 4:nr - | \
            csvtk -H -t join -f "1;3" - {input.hla_pos} | awk -v OFS='\\t' '{{print $6,$7,$1,$2,$3,$4,$5}}' > {output.hla_decompose_geno}
            python3 scripts/mhc_vcf.py {output.hla_decompose_geno} {input.samples} | bcftools sort -o {output.hla_vcf}
        """