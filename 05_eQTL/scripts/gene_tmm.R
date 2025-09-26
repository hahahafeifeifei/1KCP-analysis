args <- commandArgs(trailingOnly = TRUE)
gene_count_file <- args[1]
gene_TPM_file <- args[2]
gene_count_out_file <- args[3]
library(edgeR);library(data.table);library(stringr)
count_data=fread(gene_count_file,head=T,stringsAsFactors=F,data.table=F)
TPM_data=fread(gene_TPM_file,head=T,stringsAsFactors=F,data.table=F)
# genes were selected based on expression thresholds of ≥0.1 TPM in ≥20% of samples and ≥6 reads (unnormalized) in ≥20% of samples
ratio_TPM_calculator <- function(expr) {
  expr <- as.numeric(expr)
  length(which(expr >= 0.1)) / length(expr)
}
ratio_count_calculator <- function(expr) {
  expr <- as.numeric(expr)
  length(which(expr >= 6)) / length(expr)
}
ratio_TPM <- apply(TPM_data[,-c(1,2)], 1, ratio_TPM_calculator)
ratio_count <- apply(count_data[,-c(1,2)], 1, ratio_count_calculator)
keep_index <- which(ratio_TPM >= 0.2 & ratio_count >= 0.2)
gene_counts_QC = count_data[keep_index,];
y<-DGEList(counts=as.matrix(gene_counts_QC[,-c(1,2)]),genes=gene_counts_QC$Name)
dge=calcNormFactors(y) 
logcpm <- cpm(dge, log=TRUE)
gene_id=gene_counts_QC$Name
TMM_counts=cbind(gene_id,logcpm)
write.table(TMM_counts,gene_count_out_file,row=F,col=T,quo=F,sep="\t")

