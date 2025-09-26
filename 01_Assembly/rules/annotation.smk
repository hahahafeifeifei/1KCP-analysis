# Rule: RepeatMasker annotation
# RepeatMasker database should contain fam39_full.7.h5.gz download from https://www.dfam.org/releases/Dfam_3.9/families/FamDB/dfam39_full.7.h5.gz
rule repeatmasker:
    input:
        fa='results/{sample}/{sample}.{haplotype}.fasta'
    output:
        rmsk_raw_out=temp('results/{sample}/{sample}.{haplotype}.fasta.out'),
        rmsk_out='results/{sample}/{sample}.{haplotype}.repeatmasker.out',
        rmsk_bed='results/{sample}/{sample}.{haplotype}.repeatmasker.bed'
    params:
        rmsk_dir='results/{sample}'
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=30
    threads: 8
    shell:
        """
        RepeatMasker -s -xsmall -e ncbi -pa 2 -species human -dir {params.rmsk_dir} {input.fa}
        cp {output.rmsk_raw_out} {output.rmsk_out} 
        cat {output.rmsk_out} | awk -v OFS='\\t' 'NR>3{{print $5,$6-1,$7,$11}}' > {output.rmsk_bed}
        """

# Rule: TRF annotation
rule trf:
    input:
        fa='results/{sample}/{sample}.{haplotype}.fasta'
    output:
        trf_dat='results/{sample}/{sample}.{haplotype}.trf.dat.gz',
        trf_bed='results/{sample}/{sample}.{haplotype}.trf.bed'
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=30
    shell:
        """
        trf {input.fa} 2 7 7 80 10 50 15 -l 25 -h -ngs | bgzip > {output.trf_dat}
        zcat {output.trf_dat} | awk -v OFS='\\t' '{{if(substr($1,1,1)==\"@\") contig=substr($1,2,length($1)); else print contig,$2-1,$3,\"TR\"}}' > {output.trf_bed}
        """

# Rule: Liftoff annotation
rule liftoff:
    input:
        fa='results/{sample}/{sample}.{haplotype}.fasta',
        gencode=config['gencode'],
        reference=config['reference']
    output:
        liftoff_gff='results/{sample}/{sample}.{haplotype}.liftoff.gff_polished'
    params:
        liftoff_dir='results/{sample}/{sample}.{haplotype}.liftoff',
        prefix='results/{sample}/{sample}.{haplotype}.liftoff.gff'
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=10
    threads: 8
    shell:
        """
        liftoff -g {input.gencode} -sc 0.95 -copies -polish -o {params.prefix} -dir {params.liftoff_dir} -p {threads} {input.fa} {input.reference}
        """

# Rule: Hisat2 alignment
rule hisat2:
    input:
        fa='results/{sample}/{sample}.{haplotype}.fasta',
        rna_fq1 = config['rna_fq1'],
        rna_fq2 = config['rna_fq2']
    output:
        bam=temp('results/{sample}/{sample}.{haplotype}.hisat2.bam'),
        bai=temp('results/{sample}/{sample}.{haplotype}.hisat2.bam.bai')
    params:
        prefix='results/{sample}/{sample}.{haplotype}'
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=10
    threads: 8
    shell:
        """
        hisat2-build -p {threads} {input.fa} {params.prefix}
        hisat2 -x {params.prefix} -p {threads} -1 {input.rna_fq1} -2 {input.rna_fq2} | \
        samtools view -Shb | samtools sort -@ {threads} -o {output.bam}
        samtools index -@ {threads} {output.bam}
        rm *ht2
        """

# Rule: Stringtie transcripts assembly and filtering
rule stringtie:
    input:
        bam='results/{sample}/{sample}.{haplotype}.hisat2.bam',
        bai='results/{sample}/{sample}.{haplotype}.hisat2.bam.bai',
        liftoff_gff='results/{sample}/{sample}.{haplotype}.liftoff.gff_polished'
    output:
        stringtie_gtf=temp('results/{sample}/{sample}.{haplotype}.stringtie.gtf'),
        stringtie_filter_gtf=temp('results/{sample}/{sample}.{haplotype}.stringtie.filter.gtf')
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=10
    threads: 8
    shell:
        """
        stringtie -p {threads} -o {output.stringtie_gtf} {input.bam} -G {input.liftoff_gff}
        python3 scripts/gtf_filter.py {output.stringtie_gtf} {output.stringtie_filter_gtf}
        """

# Rule: Novel transcripts extraction
rule novel_transcripts_extraction:
    input:
        fa='results/{sample}/{sample}.{haplotype}.fasta',
        liftoff_gff='results/{sample}/{sample}.{haplotype}.liftoff.gff_polished',
        stringtie_gtf='results/{sample}/{sample}.{haplotype}.stringtie.filter.gtf'
    output:
        tmap=temp('results/{sample}/{sample}.{haplotype}.gffcompare.{sample}.{haplotype}.stringtie.filter.gtf.tmap'),
        novel_transcripts_list=temp('results/{sample}/{sample}.{haplotype}.novel.transcripts'),
        transcripts_fa=temp('results/{sample}/{sample}.{haplotype}.transcripts.fa'),
        novel_transcripts_fa=temp('results/{sample}/{sample}.{haplotype}.novel.transcripts.fa')
    params:
        prefix='results/{sample}/{sample}.{haplotype}.gffcompare'
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=10
    shell:
        """
        gffcompare -r {input.liftoff_gff} {input.stringtie_gtf} -G -o {params.prefix}
        awk '{{if($3==\"i\" || $3==\"u\" || $3==\"x\")print $5}}' {output.tmap} > {output.novel_transcripts_list}
        gffread -g {input.fa} -w {output.transcripts_fa} {input.stringtie_gtf}
        seqkit grep -f {output.novel_transcripts_list} {output.transcripts_fa} > {output.novel_transcripts_fa} 
        """

# Rule: CPC2 classification
rule cpc2:
    input:
        novel_transcripts_fa='results/{sample}/{sample}.{haplotype}.novel.transcripts.fa'
    output:
        cpc_txt=temp('results/{sample}/{sample}.{haplotype}.cpc.txt'),
        cpc_known_txt=temp('results/{sample}/{sample}.{haplotype}.cpc.known.txt'),
        novel_nc_transcripts='results/{sample}/{sample}.{haplotype}.novel.nc.transcripts',
        novel_pc_raw_transcripts=temp('results/{sample}/{sample}.{haplotype}.novel.pc_raw.transcripts')
    params:
        prefix='results/{sample}/{sample}.{haplotype}.cpc'
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=10
    shell:
        """
        scripts/cpc2/bin/CPC2.py -i {input.novel_transcripts_fa} -o {params.prefix}
        awk 'NR>1{{split($1,a,".");print a[1]"."a[2]"\\t"$8}}' {output.cpc_txt} | sort | uniq | awk '{{print $1}}' | sort | uniq -c | awk '{{if($1==1)print $2}}' > {output.cpc_known_txt}
        awk '{{if($2>200 && $3*3<100 && $8=="noncoding"){{split($1,a,".");print a[1]"."a[2]"\\t"$1}}}}' {output.cpc_txt} | csvtk -H -t join -f 1 - {output.cpc_known_txt} | awk '{{print $2}}' > {output.novel_nc_transcripts}
        awk '{{if($8=="coding")print $1}}' {output.cpc_txt} > {output.novel_pc_raw_transcripts}
        """

# Rule: GeneMarker-ST for filtering novel protein-coding transcripts
rule gmst:
    input:
        transcripts_fa='results/{sample}/{sample}.{haplotype}.transcripts.fa',
        novel_pc_raw_transcripts='results/{sample}/{sample}.{haplotype}.novel.pc_raw.transcripts'
    output:
        gmst_gff=temp('results/{sample}/{sample}.{haplotype}.gmst.gff'),
        gmst_gff_faa=temp('results/{sample}/{sample}.{haplotype}.gmst.gff.faa'),
        novel_pc_filter1_transcripts=temp('results/{sample}/{sample}.{haplotype}.novel.pc_filter1.transcripts'),
        gmst_filter1_gff=temp('results/{sample}/{sample}.{haplotype}.gmst.filter1.gff')
    params:
        prefix="results/{sample}.{haplotype}/gmst/gmst"
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=10
    shell:
        """
        scripts/gmst/gmst.pl --format GFF --output {output.gmst_gff} --strand direct {input.transcripts_fa} --faa
        grep gene_id {output.gmst_gff} | awk '{{print $1}}' | uniq -c | awk '{{if($1==1)print $2}}' | csvtk -H -t join -f 1 - {input.novel_pc_raw_transcripts} > {output.novel_pc_filter1_transcripts}
        grep gene_id {output.gmst_gff} | csvtk -H -t join -f 1 - {output.novel_pc_filter1_transcripts} > {output.gmst_filter1_gff}
        """

# Rule: Repeatmasker filtering novel protein-coding transcripts
rule rmsk_filter:
    input:
        stringtie_filter_gtf='results/{sample}/{sample}.{haplotype}.stringtie.filter.gtf',
        gmst_filter1_gff='results/{sample}/{sample}.{haplotype}.gmst.filter1.gff',
        rmsk_bed = 'results/{sample}/{sample}.{haplotype}.repeatmasker.bed',
        gmst_gff_faa='results/{sample}/{sample}.{haplotype}.gmst.gff.faa'
    output:
        stringtie_filter1_cdsplus_gtf=temp('results/{sample}/{sample}.{haplotype}.stringtie.filter1.cds_plus.gtf'),
        novel_pc_filter2_transcripts=temp('results/{sample}/{sample}.{haplotype}.novel.pc_filter2.transcripts'),
        novel_pc_filter2_transcripts_faa=temp('results/{sample}/{sample}.{haplotype}.novel.pc_filter2.transcripts.faa')
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=10
    shell:
        """
        python3 scripts/gtf_cds_annotate.py {input.stringtie_filter_gtf} {input.gmst_filter1_gff} {output.stringtie_filter1_cdsplus_gtf}
        cat {output.stringtie_filter1_cdsplus_gtf} | awk -v FS='\\t' -v OFS='\\t' '{{if($3==\"CDS\"){{split($9,a,\";\");print $1,$4-1,$5,a[1],a[2]}}}}' | \
        bedtools intersect -a - -b {input.rmsk_bed} -wao | awk -v OFS='\\t' '{{print $7,$12}}' | csvtk -H -t summary -g 1 -f 2:sum | \
        awk '{{if($2==0)print $1}}' > {output.novel_pc_filter2_transcripts}
        seqkit grep -f {output.novel_pc_filter2_transcripts} {input.gmst_gff_faa} > {output.novel_pc_filter2_transcripts_faa}
        """

# Rule: Interproscan alignment
rule interproscan:
    input:
        stringtie_filter_gtf='results/{sample}/{sample}.{haplotype}.stringtie.filter.gtf',
        novel_pc_filter2_transcripts_faa='results/{sample}/{sample}.{haplotype}.novel.pc_filter2.transcripts.faa',
        gmst_gff='results/{sample}/{sample}.{haplotype}.gmst.gff'
    output:
        interproscan_out=temp('results/{sample}/{sample}.{haplotype}.interproscan.out'),
        novel_pc_transcripts='results/{sample}/{sample}.{haplotype}.novel.pc.transcripts',
        gmst_filter_gff=temp('results/{sample}/{sample}.{haplotype}.gmst.filter.gff'),
        stringtie_filter_cdsplus_gtf=temp('results/{sample}/{sample}.{haplotype}.stringtie.filter.cds_plus.gtf')
    params:
        temp_dir="results/{sample}/temp"
    threads: 8
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=30
    shell:
        """
        interproscan.sh -cpu {threads} -i {input.novel_pc_filter2_transcripts_faa} -o {output.interproscan_out} -appl CDD,Gene3D,HAMAP,PANTHER,Pfam,PIRSF,PIRSR,PRINTS,PROSITEPATTERNS,PROSITEPROFILES,SFLD,SMART,SUPERFAMILY,TIGRFAM -f TSV -T {params.temp_dir}
        cat {output.interproscan_out} | awk '{{print $1}}' | sort | uniq > {output.novel_pc_transcripts}
        grep gene_id {input.gmst_gff} | csvtk -H -t join -f 1 - {output.novel_pc_transcripts} > {output.gmst_filter_gff}
        python3 scripts/gtf_cds_annotate.py {input.stringtie_filter_gtf} {output.gmst_filter_gff} {output.stringtie_filter_cdsplus_gtf}
        """

# Rule: Merge transcripts annotation
rule create_transcript_types:
    input:
        novel_pc_transcripts='results/{sample}/{sample}.{haplotype}.novel.pc.transcripts',
        novel_nc_transcripts='results/{sample}/{sample}.{haplotype}.novel.nc.transcripts',
        novel_transcripts='results/{sample}/{sample}.{haplotype}.novel.transcripts',
        stringtie_filter_cdsplus_gtf='results/{sample}/{sample}.{haplotype}.stringtie.filter.cds_plus.gtf'
    output:
        raw_types=temp('results/{sample}/{sample}.{haplotype}.novel.transcripts.raw.type'),
        types='results/{sample}/{sample}.{haplotype}.novel.transcripts.type',
        stringtie_annotated_gtf='results/{sample}/{sample}.{haplotype}.stringtie.annotated.gtf'
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=10
    shell:
        """
        awk '{{print $1"\\tprotein_coding"}}' {input.novel_pc_transcripts} > {output.raw_types}
        awk '{{print $1"\\tlncRNA"}}' {input.novel_nc_transcripts} >> {output.raw_types}
        csvtk -H -t join -L -f 1 {input.novel_transcripts} {output.raw_types} --na unknown > {output.types}
        python3 scripts/gtf_annotate.py {input.stringtie_filter_cdsplus_gtf} {output.types} {output.stringtie_annotated_gtf}
        """

# Rule: Generate transcripts annotation bed
rule generate_gene_bed:
    input:
        liftoff_gff='results/{sample}/{sample}.{haplotype}.liftoff.gff_polished',
        stringtie_gtf='results/{sample}/{sample}.{haplotype}.stringtie.annotated.gtf',
    output:
        gene_bed = 'results/{sample}/{sample}.{haplotype}.gene.bed',
        pc_exon_liftoff_gff = temp('results/{sample}/{sample}.{haplotype}.pc_exon.liftoff.gff'),
        pc_intron_liftoff_bed =  temp('results/{sample}/{sample}.{haplotype}.pc_intron.liftoff.bed'),
        nc_exon_liftoff_gff =  temp('results/{sample}/{sample}.{haplotype}.nc_exon.liftoff.gff'),
        nc_intron_liftoff_bed =  temp('results/{sample}/{sample}.{haplotype}.nc_intron.liftoff.bed'),
        pc_exon_stringtie_gtf =  temp('results/{sample}/{sample}.{haplotype}.pc_exon.stringtie.gtf'),
        pc_intron_stringtie_bed =  temp('results/{sample}/{sample}.{haplotype}.pc_intron.stringtie.bed'),
        nc_exon_stringtie_gtf =  temp('results/{sample}/{sample}.{haplotype}.nc_exon.stringtie.gtf'),
        nc_intron_stringtie_bed =  temp('results/{sample}/{sample}.{haplotype}.nc_intron.stringtie.bed')
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=10
    shell:
        """
        cat {input.liftoff_gff} | grep "gene_type=protein_coding" | grep -v "pseudogene" | \
            awk '{{if($3!="gene" && $3!="start_codon" && $3!="stop_codon" && $3!="transcript" && $3!="exon")print $1"\\t"$4-1"\\t"$5"\\t"$3}}' > {output.gene_bed}

        cat {input.liftoff_gff} | grep "gene_type=protein_coding" | grep -v "pseudogene" | \
            awk '{{if($3=="exon")print $0}}' > {output.pc_exon_liftoff_gff}

        python3 scripts/gff_intron.py {output.pc_exon_liftoff_gff} {output.pc_intron_liftoff_bed}
        cat {output.pc_intron_liftoff_bed} | awk '{{print $1"\\t"$2"\\t"$3"\\tintron"}}' >> {output.gene_bed}

        cat {input.liftoff_gff} | grep "gene_type=miRNA\\|gene_type=misc_RNA\\|gene_type=rRNA\\|gene_type=sRNA\\|gene_type=snRNA\\|gene_type=snoRNA\\|gene_type=vault_RNA\\|gene_type=lncRNA" | grep -v "pseudogene" | \
            awk '{{if($3=="exon")print $1"\\t"$4-1"\\t"$5"\\tnc_exon"}}' >> {output.gene_bed}

        cat {input.liftoff_gff} | grep "gene_type=miRNA\\|gene_type=misc_RNA\\|gene_type=rRNA\\|gene_type=sRNA\\|gene_type=snRNA\\|gene_type=snoRNA\\|gene_type=vault_RNA\\|gene_type=lncRNA" | grep -v "pseudogene" | \
            awk '{{if($3=="exon")print $0}}' > {output.nc_exon_liftoff_gff}
        
        python3 scripts/gff_intron.py {output.nc_exon_liftoff_gff} {output.nc_intron_liftoff_bed}
        cat {output.nc_intron_liftoff_bed} | awk '{{print $1"\\t"$2"\\t"$3"\\tintron"}}' >> {output.gene_bed}

        grep "protein_coding" {input.stringtie_gtf} | grep -v "pseudogene" | \
            awk '{{if($3!="gene" && $3!="start_codon" && $3!="stop_codon" && $3!="transcript")print $1"\\t"$4-1"\\t"$5"\\t"$3"\\tStringtie"}}' >> {output.gene_bed}

        grep "protein_coding" {input.stringtie_gtf} | awk '{{if($3=="exon")print $0}}' > {output.pc_exon_stringtie_gtf}
        python3 scripts/gtf_intron.py {output.pc_exon_stringtie_gtf} {output.pc_intron_stringtie_bed}
        cat {output.pc_intron_stringtie_bed} | awk '{{print $1"\\t"$2"\\t"$3"\\tintron"}}' >> {output.gene_bed}

        grep "lncRNA" {input.stringtie_gtf} | awk '{{if($3=="exon")print $1"\\t"$4-1"\\t"$5"\\tnc_exon"}}' >> {output.gene_bed}

        grep "lncRNA" {input.stringtie_gtf} | awk '{{if($3=="exon")print $0}}' > {output.nc_exon_stringtie_gtf}
        python3 scripts/gtf_intron.py  {output.nc_exon_stringtie_gtf} {output.nc_intron_stringtie_bed}
        cat {output.nc_intron_stringtie_bed} | awk '{{print $1"\\t"$2"\\t"$3"\\tintron"}}' >> {output.gene_bed}
        """
    
# Rule: SEI annotation
rule sei:
    input:
        fa='results/{sample}/{sample}.{haplotype}.fasta'
    output:
        sei_h5='results/{sample}/{sample}.{haplotype}_300_embeddings.h5',
        sei_bed='results/{sample}/{sample}.{haplotype}.sei.bed'
    params:
        sei_model_path = config['sei_model_path'],
        clustervfeat_path = config['clustervfeat_path'],
        target_length = 300,
        save_dir = 'results/{sample}'
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=30
    shell:
        """
        python scripts/sei_inference.py \
            --sei_model_path {params.sei_model_path} \
            --clustervfeat_path {params.clustervfeat_path} \
            --target_length {params.target_length} \
            --save_dir {params.save_dir} \
            --fasta_path {input.fa}
        python3 scripts/readh5.py {output.sei_h5} | awk '{{gsub(/[0-9]+/, "",$3);print $0}}' > {output.sei_bed}
        """

# Rule: Annotation merge
rule merge_annotation:
    input:
        rmsk_bed='results/{sample}/{sample}.{haplotype}.repeatmasker.bed',
        trf_bed='results/{sample}/{sample}.{haplotype}.trf.bed',
        gene_bed = 'results/{sample}/{sample}.{haplotype}.gene.bed',
        sei_bed="results/{sample}/{sample}.{haplotype}.sei.bed"
    output:
        merged_raw_bed = temp('results/{sample}/{sample}.{haplotype}.type.raw.bed'),
        merged_merge_bed = temp('results/{sample}/{sample}.{haplotype}.type.merge.bed'),
        merged_final_bed = 'results/{sample}/{sample}.{haplotype}.type.final.bed'
    conda:
        '../envs/env.yml'
    resources:
        mem_gb=60,
        runtime_hrs=10
    shell:
        """
        > {output.merged_raw_bed}
        cat {input.gene_bed} | awk -v OFS='\\t' '{{print $1,$2,$3,$4}}' >> {output.merged_raw_bed}
        cat {input.sei_bed} | awk -v OFS='\\t' '{{if($5!="L" && $5!="O" && $5!="TN")print $1,$2,$3,$4}}' >> {output.merged_raw_bed}
        cat {input.rmsk_bed} | awk -v OFS='\\t' '{{split($4,a,"/");if(a[1]=="DNA" || a[1]=="LINE" || a[1]=="SINE" || a[1]=="LTR" || a[1]=="Low_complexity" || a[1]=="Retroposon")print $1,$2,$3,a[1]}}' >> {output.merged_raw_bed}
        cat {input.trf_bed} >> {output.merged_raw_bed}

        > {output.merged_merge_bed}
        cat {output.merged_raw_bed} | awk '{{print $4}}' | sort | uniq | while read line; do
            awk -v type=$line '{{if($4==type)print $0}}' {output.merged_raw_bed} | bedtools sort -i - | bedtools merge -i - | \
            awk -v type=$line '{{print $0"\\t"type}}' >> {output.merged_merge_bed}; 
        done

        bedtools sort -i {output.merged_merge_bed} > {output.merged_final_bed}
        """
