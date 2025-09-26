# Rule: assembly alignment
rule assembly_alignment:
    input:
        fa=config['fa'],
        reference=config['grch38']
    output:
        bam="results/bam/{sample}.{haplotype}.bam",
        bai="results/bam/{sample}.{haplotype}.bam.bai"
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=80,
        runtime_hrs=30
    threads: 8
    shell:
        """
            minimap2 -t {threads} -x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 -a {input.reference} {input.fa} | \
            samtools view -Shb | samtools sort -@ {threads} -o {output.bam}
            samtools index -@ {threads} {output.bam}
        """

# Rule: tr detection
rule tr_detection:
    input:
        bam="results/bam/{sample}.{haplotype}.bam",
        tr_reference=config['tr_reference']
    output:
        tr_raw_vcf=temp("results/tr/{sample}.{haplotype}.tr.raw.vcf"),
        tr_edit_vcf=temp("results/tr/{sample}.{haplotype}.tr.edit.vcf"),
        tr_vcf=temp("results/tr/{sample}.{haplotype}.tr.vcf"),
        tr_edit=temp("results/tr/{sample}.{haplotype}.tr.editdist")
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=30
    threads: 8
    shell:
        """
            vamos --contig -S -b {input.bam} -r {input.tr_reference} -o {output.tr_raw_vcf} -L 50000 -s {wildcards.sample}.{wildcards.haplotype} -t {threads}
            python3 scripts/vamos_mod_editdist.py {output.tr_raw_vcf} {output.tr_edit_vcf}
            bcftools annotate -x FORMAT/SS {output.tr_edit_vcf} | bcftools sort -o {output.tr_vcf}
            bcftools query -f "%CHROM\\t%POS\\t%END\\t[%SS]\\n" {output.tr_edit_vcf} | \
            awk -v name={wildcards.sample}.{wildcards.haplotype} -v OFS='\\t' '{{print $1,$2,$3,name,$4}}' > {output.tr_edit}
        """

#Rule: merge tr detection results
rule tr_merge:
    input:
        tr_vcfs=expand("results/tr/{sample}.{haplotype}.tr.vcf", sample=samples_list, haplotype=["hap1","hap2"]),
        tr_edits=expand("results/tr/{sample}.{haplotype}.tr.editdist", sample=samples_list, haplotype=["hap1","hap2"]),
        unreliable=config['unreliable_bed']
    output:
        tr_vcf_list=temp("results/tr.vcf.list"),
        tr_raw_vcf=temp("results/merge.tr.site.raw.vcf"),
        tr_vcf="results/merge.tr.site.vcf",
        tr_filter_vcf="results/merge.tr.site.reliable.vcf"
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=200,
        runtime_hrs=30
    threads: 8
    shell:
        """
            ls {input.tr_vcfs} > {output.tr_vcf_list}
            scripts/tryvamos.py combineVCF {output.tr_vcf_list} {output.tr_raw_vcf}
            cat {output.tr_raw_vcf} | python3 ~/software/script/graph-annotation/graph_vcf_merge_split.py merge | bcftools sort -o {output.tr_vcf}
            grep "^#" {output.tr_vcf} > {output.tr_filter_vcf}
            cat {input.tr_edits} | csvtk -H -t summary -g 1,2,3 -f 5:mean -w 4 | \
            awk -v OFS='\\t' '{{if($4>0.8) print $1,$2,$3+1}}' | bedtools intersect -a - -b {input.unreliable} -wao | \
            csvtk -H -t summary -g 1,2,3 -f 7:sum | awk -v OFS='\\t' '{{if($4/($3-$2)<0.2) print $1,$2}}' | \
            csvtk -H -t join -f 1,2 {output.tr_vcf} - >> {output.tr_filter_vcf}
        """

#Rule: identify tr outliers
rule tr_outlier:
    input:
        tr_filter_vcf="results/merge.tr.site.reliable.vcf"
    output:
        tr_outlier_range_vcf="results/merge.tr.site.reliable.outlier_range.vcf"
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=200,
        runtime_hrs=30
    threads: 8
    shell:
        """
            python3 scripts/vamos_tr_outlier.py {input.tr_filter_vcf} {output.tr_outlier_range_vcf}
        """

#Rule: decompose tr site into length and motifs
rule tr_decompose:
    input:
        tr_filter_vcf="results/merge.tr.site.reliable.vcf"
    output:
        tr_length_filter_vcf="results/merge.tr.length.reliable.vcf",
        tr_motif_filter_vcf="results/merge.tr.motif.reliable.vcf"
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=30
    threads: 8
    shell:
        """
           python3 scripts/vamos_tr_dosage.py {input.tr_filter_vcf} {output.tr_length_filter_vcf} {output.tr_motif_filter_vcf}
        """