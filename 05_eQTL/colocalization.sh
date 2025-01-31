#!/bin/bash
###parameter setting
thread=16
trait_index=$1
BBJ_summary=BBJ_phenotype.csv
gwas_name=$(sed -n "${trait_index}p" ${BBJ_summary} | awk -F "," '{print $2}' | tr -d '\r\n')
gwas_file="GWASsummary_${gwas_name}_Japanese_SakaueKanai2020.COJO.ma"

###Colocalization
prefix=1kcp.eqtl
for chr in $(seq 1 22); do
    Rscript COLOC.R --gwas ${gwas_file} --reference_freq ${prefix}.afreq --qtl_query ${prefix}.chr${chr}.query --qtl_number 1000 --out ${gwas_name}_chr${chr}.coloc
    smr --bfile ${prefix} --gwas-summary ${gwas_file} --beqtl-summary ${prefix}.chr${chr}.besd --probe-chr ${chr} --maf 0.01 --peqtl-smr 0.003 --smr-multi --out ${gwas_name}_chr${chr}
done
