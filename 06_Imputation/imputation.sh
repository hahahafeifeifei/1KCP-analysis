#!/bin/bash
###parameter setting
thread=16
vcf=1kcp.imputation.high.vcf.gz
omni=1kcp.OMNI.vcf.gz
asa=1kcp.ASA.vcf.gz
samples=high.sample.list 
bcftools view -v snps --threads 16 -R ASA.site $vcf -o $asa
bcftools view -v snps --threads 16 -R OMNI.site $vcf -o $omni

###Leave-one-out imputation
cat $samples | while read sample
do
    for chr in chr{1..22}
    do
        bcftools view --threads $thread -s ^${sample} $vcf $chr -o 1kcp.imputation.${sample}_out.$chr.vcf.gz
        bcftools view --threads $thread -s $sample $omni $chr -o $sample.OMNI.$chr.vcf.gz
        tabix $sample.OMNI.$chr.vcf.gz
        bcftools view --threads $thread -s $sample $ASA $chr -o $sample.ASA.$chr.vcf.gz
        tabix $sample.ASA.$chr.vcf.gz
        minimac4 -t $thread --compress-reference 1kcp.imputation.${sample}_out.$chr.vcf.gz -o 1kcp.imputation.${sample}_out.$chr.msav
        minimac4 --format DS,GT -t $thread 1kcp.imputation.${sample}_out.$chr.msav $sample.OMNI.$chr.vcf.gz -o $sample.OMNI.impute.$chr.vcf.gz
        minimac4 --format DS,GT -t $thread 1kcp.imputation.${sample}_out.$chr.msav $sample.ASA.$chr.vcf.gz -o $sample.ASA.impute.$chr.vcf.gz
    done

    > $sample.summary.info
    for chr in chr{1..22}
    do
        bcftools view --threads $thread -s $sample $vcf $chr -o $sample.truth.$chr.vcf.gz
        bcftools query -f "%ID\t[%GT]\n" $sample.truth.$chr.vcf.gz > $sample.truth.$chr.info
        bcftools query -f "%ID\t[%GT\t%DS]\n" $sample.OMNI.impute.$chr.vcf.gz > $sample.OMNI.$chr.info
        bcftools query -f "%ID\t[%GT\t%DS]\n" $sample.ASA.impute.$chr.vcf.gz > $sample.ASA.$chr.info
        csvtk -H -t join -f 1 $sample.truth.$chr.info $sample.ASA.$chr.info $sample.OMNI.$chr.info | awk -v sample=$sample -v OFS='\t' '{print sample,$0}' >> $sample.summary.info
    done
done

###Calculate R square
for array in {OMNI,ASA}
do
    cat *summary.info | grep -v "\.|" | grep -v "|\." | grep -v TR | awk -v OFS='\t' '{split($3,a,"|");split($4,b,"|");print $1,$2,a[1]+a[2],b[1]+b[2],$5}' > 1kcp.$array.nonTR.summary
    python3 nonTR_stat.py 1kcp.$array.nonTR.summary > 1kcp.$array.nonTR.r2
    cat *summary.info | grep -v "\.|" | grep -v "|\." | grep TR | awk -v OFS='\t' '{split($3,a,"|");split($4,b,"|");print $1,$2,a[1]+a[2],b[1]+b[2],$5}' > 1kcp.$array.TR.summary
    awk -v OFS='\t' '{split($2,a,"-");id=a[1];for(i=2;i<=length(a)-1;i++) id=id"-"a[i]; print $1,id,$2,$3,$4,$5}' 1kcp.$array.TR.summary | csvtk -H -t join -f "3;1" - ../1kcp.tr.len | awk -v OFS='\t' '{print $1,$2,$4*$7,$5*$7,$6*$7}' | csvtk -H -t summary -g 1,2 -f 3:sum,4:sum,5:sum -n 4 > 1kcp.$array.TR.loci.summary
    python3 TR_stat.py 1kcp.$array.TR.loci.summary > 1kcp.$array.TR.r2
done