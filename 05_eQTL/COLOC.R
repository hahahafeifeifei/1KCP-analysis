suppressMessages(library("optparse"))
option_list = list(
        make_option("--gwas", action="store", default=NA, type='character',
              help="Path to GWAS summary statistics [required]"),
        make_option("--reference_freq", action="store", default=NA, type='character',
              help="Path to Reference freq file [required]"),
  make_option("--qtl_query", action="store", default=NA, type='character',
              help="Path to QTL query data [required]"),
        make_option("--qtl_number", action="store", default=NA, type='integer',
              help="Sample size of QTL data [required]"),
        make_option("--out", action="store", default=NA, type='character',
              help="Path to result file")
)


opt = parse_args(OptionParser(option_list=option_list))

suppressMessages({
library("coloc");
library("dplyr");
library("data.table")
})

# ------------------ reference frequency ------------
reference_freq=fread(opt$reference_freq,head=T,stringsAsFactors=F,data.table=F)

# ------------------ gwas input ---------------------
gwas=fread(opt$gwas,head=T,stringsAsFactors=F,data.table=F)
colnames(gwas)=c("SNP","A1","A2","Freq","b","se","P","N")
rm_index=which(gwas$se==0)
if(length(rm_index)>0){
  gwas=gwas[-rm_index,]
}

# ------------------ qtl summary (cis) input ---------------------
qtl_number=as.numeric(opt$qtl_number)
qtl_query=fread(opt$qtl_query,  head=T,stringsAsFactors=F,data.table=F)
colnames(qtl_query)=c("SNP","Chr","BP","A1","A2","Freq","Probe","Probe_Chr","Probe_bp","Gene","Orientation","b","SE","p")
rm_index=which(qtl_query$SE == 0)
if(length(rm_index)>0){
  qtl_query=qtl_query[-rm_index,]
}

if(any(is.na(qtl_query$Freq))){
  print("All/Some freqs are NA in xQTL. We will use the frequency in the reference.")
  freq_na_index=which(is.na(qtl_query$Freq))

  reference_index1=match(paste0(qtl_query$SNP[freq_na_index], "_",qtl_query$A1[freq_na_index]), paste0(reference_freq$ID, "_", reference_freq$ALT), nomatch=0)
  qtl_query$Freq[freq_na_index][which(reference_index1!=0)]=reference_freq$ALT_FREQS[reference_index1]
  reference_index2=match(paste0(qtl_query$SNP[freq_na_index], "_",qtl_query$A1[freq_na_index]), paste0(reference_freq$ID, "_", reference_freq$REF), nomatch=0)
  qtl_query$Freq[freq_na_index][which(reference_index2!=0)]=1-reference_freq$ALT_FREQS[reference_index2]
}


pbs=unique(qtl_query$Probe)

# ------------------ coloc analysis -------------------------------
# GWAS trait and QTL coloc analysis based on QTL probe ------------

summary=data.frame()

for(probe in pbs){

  print(probe)
  # query snps based on probe ----
        index=which(qtl_query$Probe == probe)
        query=qtl_query[index,]

  gene_name=query$Gene[1]

  index=match(query$SNP,gwas$SNP,nomatch=0)
  query_com=query[which(index!=0),]
  gwas_com=gwas[index,]

  if(nrow(gwas_com)<10){next}

  gwas_com[,2]=toupper(gwas_com[,2]);gwas_com[,3]=toupper(gwas_com[,3])
  index1=which(query_com$A1==gwas_com[,3] & query_com$A2==gwas_com[,2])

  if(length(index1)>0){
    gwas_com[index1,2]=query_com$A1[index1];
    gwas_com[index1,3]=query_com$A2[index1];
    gwas_com[index1,5]=-gwas_com[index1,5];
    gwas_com[index1,4]=1-gwas_com[index1,4]
  }

  if(all(is.na(gwas_com[,4]))){
    print("All freqs are NA in GWAS.")
  }

  index2=which(query_com$A1==gwas_com[,2] & query_com$A2==gwas_com[,3])
  if(length(index2)>=0){
    query_com=query_com[index2,]
    gwas_com=gwas_com[index2,]
  }

  snpbuf=query_com$SNP

  if(length(snpbuf)<10){
    next
  }

  freq1=as.numeric(gwas_com[,4]);
  pval1=as.numeric(gwas_com[,7]);
  n1=gwas_com[,8];
  bzy_hat=gwas_com[,5];
  bzy_se=gwas_com[,6]

  if(any(is.na(freq1))){
        index=which(is.na(freq1) | freq1<=0 | freq1>=1)
        if(length(index)>0){
                freq1[index]=query_com$Freq[index]
        }
    std_eff1=calcu_std_b_se(bzy_hat/bzy_se,freq1,n1)
    bzy_hat=std_eff1[,1];
    bzy_se=std_eff1[,2];
  }

  data1 = data.frame(pval1, n1, freq1, bzy_hat, var_bzy=bzy_se^2,
                    type="quant", snpbuf)

  pval2=query_com$p;
        n2=qtl_number;
  freq2=query_com$Freq;
        bzx_hat=query_com$b;
  bzx_se=query_com$SE;
  z2=bzx_hat/bzx_se;
  if(length(which(is.na(freq2)))>0){
    index=which(is.na(freq2))
    freq2[index]=freq1[index]
    std_eff2=calcu_std_b_se(z2,freq2,n2)
    bzx_hat=std_eff2[,1];
    bzx_se=std_eff2[,2];
  }


  data2 = data.frame(pval2, n2, freq2, bzx_hat, var_bzx=bzx_se^2,
                   type="quant", snpbuf)
  colnames(data1) = colnames(data2) = c("pvalues", "N", "MAF",
                                      "beta", "varbeta", "type",
                                      "snp")

  index=which(!is.na(data1$MAF) & !is.na(data2$MAF))
  data1=data1[index,];data2=data2[index,]

  index=which(data1$MAF > 0 & data2$MAF >0)
  data1=data1[index,];data2=data2[index,]

  index=which(data1$MAF < 1 & data2$MAF <1)
  data1=data1[index,];data2=data2[index,]

  if(nrow(data1) < 10){
    next
  }
  # # data1 and data2 p,n,snp inforation
  # test1=data1[,c("pvalues", "MAF","N","type","snp")]
  # test2=data2[,c("pvalues", "MAF","N","type","snp")]
  # test1=as.list(test1);test2=as.list(test2)
  # test1$type="quant";test2$type="quant"
  # coloc_res=coloc.abf(test1, test2)

  # data1 and data2 all information
  data1=as.list(data1);data2=as.list(data2)
  data1$type="quant";data2$type="quant"
  resbuf=coloc.abf(data1, data2)

  pp0 = as.numeric(resbuf$summary[2])
  pp1 = as.numeric(resbuf$summary[3])
  pp2 = as.numeric(resbuf$summary[4])
  pp3 = as.numeric(resbuf$summary[5])
  pp4 = as.numeric(resbuf$summary[6])
  PP3_PP4=pp3+pp4
  res = data.frame(probe,gene_name, pp0, pp1, pp2, pp3, pp4, PP3_PP4)
  summary=rbind(summary,res)

}
write.table(summary, opt$out, row=F,col=T,quo=F,sep="\t")