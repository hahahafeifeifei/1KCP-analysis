# Rule: build STAR index
rule index:
    input:
        reference=config['grch38'],
        gencode=config['gencode']
    output:
        index_dir=directory("results/star_index")
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=30
    threads: 16
    shell:
        """
            STAR --runMode genomeGenerate --genomeDir {output.index_dir} --runThreadN {threads} \
            --genomeFastaFiles {input.reference} --sjdbGTFfile {input.gencode} --sjdbOverhang 149
        """

# Rule: collapse transcripts annotation into gene annotation
rule collapse_anno:
    input:
        gencode=config['gencode']
    output:
        collapse_gencode="results/gencode.collapse.gtf",
        gencode_bed="results/gencode.gene.bed"
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=30
    threads: 8
    shell:
        """
            python3 scripts/collapse_annotation.py {input.gencode} {output.collapse_gencode}
            awk -v OFS='\\t' '{{if($3=="gene") print $1,$4-1,$5,$10}}' {input.gencode} | sed "s/\"\|;//g" > {output.gencode_bed}
        """

# Rule: gene expression quantification
rule quantification:
    input:
        index_dir="results/star_index",
        gencode="results/gencode.collapse.gtf",
        rna_fq1 = config['rna_fq1'],
        rna_fq2 = config['rna_fq2']
    output:
        star_bam=temp("results/quantification/{sample}_Aligned.out.bam"),
        star_sort_bam=temp("results/quantification/{sample}_Aligned.sort.out.bam"),
        star_dedup_bam=temp("results/quantification/{sample}_Aligned.dedup.out.bam"),
        star_dedup_txt=temp("results/quantification/{sample}_Aligned.dedup.txt"),
        count="results/quantification/{sample}.gene_reads.gct",
        tpm="results/quantification/{sample}.gene_tpm.gct"
    params:
        star_dir="results/quantification"
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=30
    threads: 8
    shell:
        """
            STAR --runMode alignReads \
            --runThreadN {threads} \
            --genomeDir {input.index_dir}\
            --twopassMode Basic \
            --outFilterMultimapNmax 20 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverLmax 0.1 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --outFilterType BySJout \
            --outFilterScoreMinOverLread 0.33 \
            --outFilterMatchNminOverLread 0.33 \
            --limitSjdbInsertNsj 1200000 \
            --readFilesIn {input.rna_fq1} {input.rna_fq2} \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.star_dir}/{wildcards.sample}_ \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs None \
            --alignSoftClipAtReferenceEnds Yes \
            --quantMode TranscriptomeSAM GeneCounts \
            --outSAMtype BAM Unsorted \
            --outSAMunmapped Within \
            --genomeLoad NoSharedMemory \
            --chimSegmentMin 15 \
            --chimJunctionOverhangMin 15 \
            --chimOutType WithinBAM SoftClip \
            --chimMainSegmentMultNmax 1 \
            --outSAMattributes NH HI AS nM NM ch \
            --outSAMattrRGline ID:rg1 SM:sm1
            samtools sort -@ {threads} -o {output.star_sort_bam} {output.star_bam}
            picard MarkDuplicates I={output.star_sort_bam} O={output.star_dedup_bam} M={output.star_dedup_txt}
            rnaseqc {input.gencode} {output.star_dedup_bam} {params.star_dir} -s {wildcards.sample} -vv --stranded=rf
        """

# Rule: merge expression
rule summary:
    input:
        counts=expand("results/quantification/{sample}.gene_reads.gct", sample=samples_list),
        tpms=expand("results/quantification/{sample}.gene_tpm.gct", sample=samples_list)
    output:
        merge_count="results/merge.gene_tpm.matrix",
        merge_tpm="results/merge.gene_count.matrix"
    params:
        star_dir="results/quantification"
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=30
    threads: 8
    shell:
        """
            paste {input.counts} | awk 'NR>2 {{printf $1"\\t"$2;for(i=1;i<=NF/3;i++) {{j=3*i;printf "\\t"$j}};print ""}}' > {output.merge_count}
            paste {input.tpms} | awk 'NR>2 {{printf $1"\\t"$2;for(i=1;i<=NF/3;i++) {{j=3*i;printf "\\t"$j}};print ""}}' > {output.merge_tpm}
        """

# Rule: TMM normalization
rule tmm_normalization:
    input:
        merge_count="results/merge.gene_tpm.matrix",
        merge_tpm="results/merge.gene_count.matrix",
        cov=contig['cov']
    output:
        merge_count_tmm="results/merge.gene_count_tmm.matrix"
        merge_count_tmm_peer="results/merge.gene_count_peer.matrix"
    conda:
        '../envs/peer.yml'
    resources:
        mem_gb=60,
        runtime_hrs=30
    threads: 8
    shell:
        """
            Rscript gene_tmm.R {input.merge_count} {input.merge_tpm} {output.merge_count_tmm}
            Rscript gene_peer.R {input.cov} {output.merge_count_tmm} {output.merge_count_tmm_peer}
        """

# Rule: tensorqtl format
rule tensorqtl_format:
    input:
        merge_count_tmm_peer="results/merge.gene_count_peer.matrix",
        sample_list=config['sample_list_file'],
        gencode_bed="results/gencode.gene.bed"
    output:
        bed_tmp=temp("results/bed.tmp"),
        sample_tmp=temp("results/sample.tmp"),
        title_tmp=temp("results/title.tmp"),
        merge_list="results/merge.gene.list"
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=30
    threads: 8
    shell:
        """
            csvtk -t transpose {input.merge_count_tmm_peer} | awk '{print $1}' | csvtk -H -t join -f "1;4" - {input.gencode_bed} | awk -v OFS='\\t' '{{print $2,$3-1,$3,$1}}' > {output.bed_tmp}
            csvtk -H -t join -f 1 {input.sample_list} {input.merge_count_tmm_peer} | csvtk -t transpose - | awk 'NR>3{{print $0}}' > {output.sample_tmp}
            csvtk -H -t join -f 1 {input.sample_list} {input.merge_count_tmm_peer} | csvtk -t transpose - | awk 'NR==1{{print "chr\\tstart\\tend\\tgene_id\\t"$0}}' > {output.title_tmp}
            paste {output.bed_tmp} {output.sample_tmp} | bedtools sort | cat {output.title_tmp} - > {output.merge_list}
        """
