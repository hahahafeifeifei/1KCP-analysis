library(data.table)
library(stringr)
BASE_DIRT="./"
QC_matrix=read.csv("1kcp_sample_rnaseq_metrics.tsv",head=T,stringsAsFactors=F,sep="\t")
index=which(!QC_matrix$is_rm_low_quality_samples & QC_matrix$D_statistics>=0.85 & QC_matrix$Mapped.Unique.Reads > 10000000)
RNA_ID=QC_matrix$sample_id[index]
rna_id=RNA_ID[1]
gene_count_infle=fread(paste0(BASE_DIRT,rna_id,"-RNA.gene_reads.gct"),head=T,stringsAsFactors=F,data.table=F);
gene_tpm_infle=fread(paste0(BASE_DIRT,rna_id,"-RNA.gene_tpm.gct"),head=T,stringsAsFactors=F,data.table=F);
Counts_gene = gene_count_infle[,1:2]
TPM_gene = gene_tpm_infle[,1:2]
for(ID in RNA_ID){
        print(ID)
        rna_id=ID
        gene_count_infle=read.table(paste0(BASE_DIRT,rna_id,"-RNA.gene_reads.gct"),head=T,skip=2);
        gene_tpm_infle=fread(paste0(BASE_DIRT,rna_id,"-RNA.gene_tpm.gct"),head=T,stringsAsFactors=F,data.table=F);
        Counts_gene=cbind(Counts_gene,gene_count_infle[,3])
        TPM_gene=cbind(TPM_gene,gene_tpm_infle[,3])
}
colnames(Counts_gene)[-c(1:2)]=RNA_ID
colnames(TPM_gene)[-c(1:2)]=RNA_ID
fwrite(Counts_gene,"1kcp_gene_count.matrix",row=F,col=T,quo=F,sep="\t")
fwrite(TPM_gene,"1kcp_gene_TPM.matrix",row=F,col=T,quo=F,sep="\t")