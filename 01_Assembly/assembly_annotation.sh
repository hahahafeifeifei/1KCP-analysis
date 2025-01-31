#!/bin/bash
###parameter setting
sample=$1
thread=16
gtf=gencode.v40.GRCh38.annotation.gff3
prefix=$sample.hap1
fa=$prefix.fasta
rna_fq1=${sample}-RNA.R1.fastq.gz
rna_fq2=${sample}-RNA.R2.fastq.gz

###Repeat annotation
RepeatMasker -s -xsmall -e ncbi -pa 2 -species human -dir . $fa
zcat *.repeatmasker.out.gz | awk -v OFS='\t' '{print $5,$6-1,$7,$11}' > $prefix.repeatmasker.bed
trf $fa 2 7 7 80 10 50 15 -l 25 -h -ngs | bgzip > $prefix.trf.dat.gz
zcat $prefix.trf.dat.gz | awk -v OFS='\t' '{if(substr($1,1,1)=="@") contig=substr($1,2,length($1));print contig,$1,$2,"TR"}}' > $prefix.trf.bed
biser --gc-heap 10G --keep-contigs -o $prefix.biser.bed -t 8 $fa

###Gene annotation
####Liftoff
liftoff -g $gtf -sc 0.95 -copies -polish -o $prefix.liftoff.gff -dir $prefix -p $thread $fa $ref -u $prefix.unmapped.txt
####Stringtie
hisat2-build -p $thread $fa $prefix
hisat2 -x $prefix -p $thread -1 $rna_fq1 -2 $rna_fq2 2> $prefix.mapping_stat.txt | samtools sort /dev/stdin -@ $thread -o $prefix.hisat2.bam
samtools index -@ $thread $prefix.hisat2.bam
stringtie -p $thread -o $prefix.stringtie.gtf $prefix.hisat2.bam -G $prefix.liftoff.gff_polished
python3 gtf_filter.py $prefix.stringtie.gtf $prefix.stringtie.filter.gtf
####Stringtie novel genes
gffcompare -r $prefix.liftoff.gff_polished $prefix.stringtie.filter.gtf -G -o $prefix.gffcompare
awk '{if($3=="i" || $3=="u" || $3=="x")print $5}' $prefix.gffcompare.$prefix.stringtie.filter.gtf.tmap > $prefix.novel.transcripts
cat $prefix.stringtie.filter.gtf | gffread -g $fa -w $prefix.transcripts.fa
seqkit grep -f $prefix.novel.transcripts $prefix.transcripts.fa > $prefix.novel.transcripts.fa
####Classifying novel genes into non-coding and protein-coding genes
CPC2.py -i $prefix.novel.transcripts.fa -o $prefix.cpc
awk 'NR>1{split($1,a,".");print a[1]"."a[2]"\t"$8}' $prefix.cpc.txt | sort | uniq | awk '{print $1}' | sort | uniq -c | awk '{if($1==1)print $2}' > $prefix.cpc.known
####Non-coding genes
awk '{if($2>200 && $3*3<100 && $8=="noncoding"){split($1,a,".");print a[1]"."a[2]"\t"$1}}' $prefix.cpc.txt | csvtk -H -t join -f 1 - $prefix.cpc.known | awk '{print $2}' > $prefix.novel.nc.transcripts
####Protein-coding genes 
awk '{if($8=="coding")print $1}'  $prefix.cpc.txt > $prefix.novel.pc.raw.transcripts
#####Filtering step1: GeneMarkS-T validation
gmst.pl --format GFF --output $prefix.gmst.gff --strand direct $prefix.transcripts.fa --faa
grep gene_id $prefix.gmst.gff | awk '{print $1}' | uniq -c | awk '{if($1==1)print $2}' | csvtk -H -t join -f 1 - $prefix.novel.pc.raw.transcripts > $prefix.novel.pc.filter1.transcripts
grep gene_id $prefix.gmst.gff | csvtk -H -t join -f 1 -  $prefix.novel.pc.filter1.transcripts >  $prefix.gmst.filter1.gff
#####Filtering step2: Not overlap with Repeats
python3 gtf_cds_annotate.py $prefix.stringtie.filter.gtf.gz $prefix.gmst.filter1.gff $prefix.stringtie.filter.filter1.cds_plus.gtf
cat $prefix.stringtie.filter.filter1.cds_plus.gtf | awk -v FS='\t' -v OFS='\t' '{if($3=="CDS"){split($9,a,";");print $1,$4-1,$5,a[1],a[2]}}' | bedtools intersect -a - -b $prefix.repeatmasker.bed -wao | awk -v OFS='\t' '{print $7,$12}' | csvtk -H -t summary -g 1 -f 2:sum | awk '{if($2==0)print $1}' > $prefix.novel.pc.filter2.transcripts
seqkit grep -f $prefix.novel.pc.filter2.transcripts $prefix.gmst.gff.faa > $prefix.novel.pc.filter2.transcripts.faa
#####Filtering step3: Interproscan alignment
interproscan.sh -cpu $thread -i $prefix.novel.pc.filter2.transcripts.faa -o $prefix.interproscan.out -appl CDD,Gene3D,HAMAP,PANTHER,Pfam,PIRSF,PIRSR,PRINTS,PROSITEPATTERNS,PROSITEPROFILES,SFLD,SMART,SUPERFAMILY,TIGRFAM -f TSV
cat $prefix.interproscan.out | awk '{print $1}' | sort | uniq > $prefix.novel.pc.transcripts
grep gene_id $prefix.gmst.gff | csvtk -H -t join -f 1 -  $prefix.novel.pc.transcripts > $prefix.gmst.filter.gff
python3 gtf_cds_annotate.py $prefix.stringtie.filter.gtf.gz $prefix.gmst.filter.gff $prefix.stringtie.filter.cds_plus.gtf
####Merging novel transcripts
awk '{print $1"\tprotein_coding"}' $prefix.novel.pc.transcripts > $prefix.novel.transcripts.raw.type
awk '{print $1"\tlncRNA"}' $prefix.novel.nc.transcripts >> $prefix.novel.transcripts.raw.type
csvtk -H -t join -L -f 1 $prefix.novel.transcripts $prefix.novel.transcripts.raw.type --na unknown > $prefix.novel.transcripts.type
python3 gtf_annotate.py $prefix.stringtie.filter.cds_plus.gtf $prefix.novel.transcripts.type $prefix.stringtie.annotated.gtf

###Epigenome annotation
####SEI

####Enformer
