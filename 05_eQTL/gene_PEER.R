num=60
library(peer);library(data.table);library(impute);library('fastDummies');library(stringr)
source("CommonFunc.r")
base_dirt="./"
cov_file="1kcp.sex_age.txt"
expr=fread(paste0(base_dirt,"1kcp_gene_counts_TMM_QC.matrix"),head=T,stringsAsFactors=F,data.table=F)
expr_data=as.data.frame(t(expr[,-1]))
expr_data3=expr_data
cov_data = read.table(cov_file,head=T,stringsAsFactors=F,sep="\t");
index=match(cov_data$Participant_id,str_split_fixed(rownames(expr_data3),"-RNA",Inf)[,1],nomatch=0)
expr_data4 = expr_data3[index,]
covs=cov_data[which(index!=0),]
Gender <- model.matrix(~as.factor(covs$Gender))[,-1]
#Batch <- model.matrix(~as.factor(covs$Batch))[,-1]
#cov_data2=data.frame(covs[,c("Age")],Gender,Batch)
cov_data2=data.frame(covs[,c("Age")],Gender)
model = PEER()
PEER_setCovariates(model, as.matrix(cov_data2))
PEER_setPhenoMean(model,as.matrix(expr_data4))
PEER_setAdd_mean(model, TRUE)
PEER_setNk(model,num)
PEER_getNk(model)
PEER_update(model)
PEER_getResidualVars(model)
ID=covs$Participant_id
pheno=PEER_getResiduals(model)
expr_data5=apply(pheno,2,RINT)
out2=cbind(ID,ID,expr_data5)
colnames(out2)=c("FID","IID",expr$gene_id)
write.table(out2,paste0(base_dirt,"1kcp_Counts_TMM_QC_cov_",num,"PEER.residuals.RINT.txt"),row=F,col=T,quo=F,sep="\t")