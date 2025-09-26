# Rule: extract aa sequence
rule gtf_aa:
    input: 
        gencode=config['gencode'],
        reference=config['grch38'],
        unreliable=config['unreliable_bed']
    output:
        gencode_coding_canonical=temp("results/gencode.canonical.coding.gtf"),
        gencode_coding_canonical_id=temp("results/gencode.canonical.coding.id"),
        gencode_raw_faa=temp("results/gencode.raw.faa"),
        gencode_faa="results/gencode.faa",
        gencode_filter_gene="results/gencode.filter.gene"
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=30
    threads: 8
    shell:
        """
            grep "Ensembl_canonical" {input.gencode} | grep "protein_coding" > {output.gencode_coding_canonical}
            awk '{{if($3=="transcript") print $10"\\t"$12"\\t"$16}}' {output.gencode_coding_canonical} | sed "s/;\|\\"//g" | awk '{{print $2"\\t"$3":"$1}}' > {output.gencode_coding_canonical_id} 
            gffread -g {input.reference} -y {output.gencode_raw_faa} {output.gencode_coding_canonical}
            seqkit fx2tab {output.gencode_raw_faa} | awk '{{print $1"\\t"$NF}}' | csvtk -H -t join -f 1 {output.gencode_coding_canonical_id} - | awk '{{print ">"$2"\\n"$3}}' > {output.gencode_faa}
            grep protein_coding {input.gencode} | awk -v OFS='\\t' '{{if($3=="gene") print $1,$4-1,$5,$14}}' | \
            sed "s/;\|\\"//g" | bedtools intersect -a - -b {input.unreliable} -wao | csvtk -H -t summary -g 1,2,3,4 -f 8:sum | \
            awk '{{if($5/($3-$2)>0.05)print $4}}' | sort -u > {output.gencode_filter_gene}
        """


# Rule: miniprot alignment
rule miniprot:
    input:
        fa=config['fa'],
        gencode_faa="results/gencode.faa"
    output:
        miniprot_paf=temp("results/{sample}.{haplotype}.paf")
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=30
    threads: 8
    shell:
        """
            miniprot --outs=0.97 --no-cs -Iu -t {threads} {input.fa} {input.gencode_faa} > {output.miniprot_paf}
        """

#Rule: pangene analysis
rule pangene:
    input:
        miniprot_pafs=expand("results/{sample}.{haplotype}.paf", sample=samples_list, haplotype=["hap1","hap2"]),
        gencode_filter_gene="results/gencode.filter.gene"
    output:
        pangene_gfa="results/merge.pangene.gfa"
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=30
    threads: 8
    shell:
        """
            pangene {input.miniprot_pafs} -p 0 -X @{input.gencode_filter_gene} > {output.pangene_gfa}
        """
