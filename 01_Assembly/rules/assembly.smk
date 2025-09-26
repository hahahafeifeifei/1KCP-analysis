# Define a list of microbial types for contig removal
nonhuman_list = ['bacteria', 'viruses_and_viroids', 'archaea', 'fungi']

# Rule: Hifiasm Assembly
rule hifiasm_assembly:
    input:
        hic_fq1=config['hic_fq1'],
        hic_fq2=config['hic_fq2'],
        hifi_fq=config['hifi_fq']
    output:
        fa_hap1=temp('results/{sample}/{sample}.hifiasm.hap1.fasta'),
        fa_hap2=temp('results/{sample}/{sample}.hifiasm.hap2.fasta')
    params:
        prefix='results/{sample}/{sample}.hifiasm'
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=lambda wildcards, attempt: 100 * attempt,
        runtime_hrs=lambda wildcards, attempt: 48 * attempt
    threads: 16
    shell:
        """
        hifiasm -o {params.prefix} -t {threads} --h1 {input.hic_fq1} --h2 {input.hic_fq2} {input.hifi_fq}
        grep ^S {output.gfa_hap1} | awk '{{print \">\"$2\"\\n\"$3}}' > {output.fa_hap1}
        grep ^S {output.gfa_hap2} | awk '{{print \">\"$2\"\\n\"$3}}' > {output.fa_hap2}
        rm {params.prefix}.*gfa {params.prefix}.*bin
        """

# Rule: Align HiFi reads to reference
rule filter_chrm:
    input:
        hifi_fq=config['hifi_fq'],
        reference=config['reference']
    output:
        bam=temp('results/{sample}/{sample}.bam'),
        bai=temp('results/{sample}/{sample}.bam.bai')
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=10
    threads: 8
    shell:
        """
        minimap2 -t {threads} -x map-hifi -a -Y -L --eqx --cs {input.reference} {input.hifi_fq} | \
        samtools view -Shb | samtools sort -@ {threads} -o {output.bam} -
        samtools index -@ {threads} {output.bam}
        """

# Rule: Unicycler assembly for chrM
rule unicycler_chrm:
    input:
        bam='results/{sample}/{sample}.bam',
    output:
        chrm_fq=temp('results/{sample}/{sample}.chrM.fastq'),
        chrm_fasta=temp('results/{sample}/{sample}.unicycler/assembly.fasta')
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=30,
        runtime_hrs=4
    shell:
        """
        samtools view -F 256 -Shb {input.bam} chrM | samtools fastq - > {output.chrm_fq}
        unicycler -l {output.chrm_fq} -o results/{wildcards.sample}/{wildcards.sample}.unicycler
        > {output.chrm_fasta}
        """

# Rule: Merge chrM assembly into Hifiasm outputs
rule merge_chrm:
    input:
        fa_hap1='results/{sample}/{sample}.hifiasm.hap1.fasta',
        fa_hap2='results/{sample}/{sample}.hifiasm.hap2.fasta',
        chrm_fasta='results/{sample}/{sample}.unicycler/assembly.fasta',
        reference=config['reference']
    output:
        chrm_hap1_origin=temp('results/{sample}/{sample}.hap1.chrM.contig'),
        chrm_hap2_origin=temp('results/{sample}/{sample}.hap2.chrM.contig'),
        fa_hap1_merged=temp('results/{sample}/{sample}.merge.hap1.fasta'),
        fa_hap2_merged=temp('results/{sample}/{sample}.merge.hap2.fasta')
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=10
    threads: 8
    shell:
        """
        minimap2 -t {threads} -x asm5 -a -Y -L --eqx --cs {input.reference} {input.fa_hap1} | \
        samtools view -S -F 256 | awk '{{if($3=="chrM") print $1}}' > {output.chrm_hap1_origin} 
        seqkit grep -v -f {output.chrm_hap1_origin} {input.fa_hap1} | cat - {input.chrm_fasta} > {output.fa_hap1_merged}

        minimap2 -t {threads} -x asm5 -a -Y -L --eqx --cs {input.reference} {input.fa_hap2} | \
        samtools view -S -F 256 | awk '{{if($3=="chrM") print $1}}' > {output.chrm_hap2_origin} 
        seqkit grep -v -f {output.chrm_hap2_origin}  {input.fa_hap2} > {output.fa_hap2_merged}
        """

# Rule: Mask adapter sequences
rule mask_adaptor:
    input:
        fa_hap1_merged='results/{sample}/{sample}.merge.hap1.fasta',
        fa_hap2_merged='results/{sample}/{sample}.merge.hap2.fasta',
        adaptor=config['adaptor']
    output:
        fa_hap1_adaptor='results/{sample}/{sample}.adaptor.hap1.bed',
        fa_hap2_adaptor='results/{sample}/{sample}.adaptor.hap2.bed',
        fa_hap1_masked=temp('results/{sample}/{sample}.masked.hap1.fasta'),
        fa_hap2_masked=temp('results/{sample}/{sample}.masked.hap2.fasta')
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=10
    threads: 8
    shell:
        """
        minimap2 -t {threads} -cxsr -f5000 -N2000 -secondary=yes --cs {input.fa_hap1_merged} {input.adaptor} | \
        awk -v OFS='\t' '{{print$6,$8,$9}}' | bedtools sort -i - | bedtools merge -i - >  {output.fa_hap1_adaptor}
        bedtools maskfasta -fi {input.fa_hap1_merged} -bed {output.fa_hap1_adaptor} -fo {output.fa_hap1_masked}

        minimap2 -t {threads} -cxsr -f5000 -N2000 -secondary=yes --cs {input.fa_hap2_merged} {input.adaptor} | \
        awk -v OFS='\t' '{{print$6,$8,$9}}' | bedtools sort -i - | bedtools merge -i - >  {output.fa_hap2_adaptor}
        bedtools maskfasta -fi {input.fa_hap2_merged} -bed {output.fa_hap2_adaptor} -fo {output.fa_hap2_masked}
        """

# Rule: Detect non-human contigs
rule blast_contigs:
    input:
        fa_masked='results/{sample}/{sample}.masked.{haplotype}.fasta',
        refseq=config['refseq_dir']+'/{type}/{type}.ndb'
    output:
        blast=temp('results/{sample}/{sample}.{type}.{haplotype}.blast')
    params:
        refseq_prefix=config['refseq_dir']+'/{type}/{type}'
    wildcard_constraints:
        type = '|'.join(nonhuman_list)
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=10
    threads: 8
    shell:
        """
        blastn -query {input.fa_masked} -db {params.refseq_prefix} \
               -task megablast -word_size 28 -best_hit_overhang 0.1 \
               -best_hit_score_edge 0.1 -dust yes -evalue 0.0001 \
               -min_raw_gapped_score 100 -penalty -5 -perc_identity 98.0 \
               -soft_masking true -out {output.blast} -num_threads {threads} -outfmt 7
        """

# Rule: Remove non-human contigs
rule remove_contigs:
    input:
        blast=expand('results/{{sample}}/{{sample}}.{type}.{{haplotype}}.blast', type=nonhuman_list),
        fa_masked='results/{sample}/{sample}.masked.{haplotype}.fasta',
        accession=config['nonhuman_accession']
    output:
        contig_list='results/{sample}/{sample}.{haplotype}.non-human.contig',
        fa_final='results/{sample}/{sample}.{haplotype}.fasta'
    wildcard_constraints:
        sample = '|'.join(samples_list)
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=10
    shell:
        """
        cat {input.blast} | grep -v "^#" | (grep -f {input.accession} - || true) | awk '{{print $1}}' | sort -u > {output.contig_list}
        seqkit grep -v -f {output.contig_list} {input.fa_masked} > {output.fa_final} --quiet
        """

