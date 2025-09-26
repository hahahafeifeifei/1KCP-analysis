import pandas as pd
import torch
import tensorqtl
from tensorqtl import pgen, cis, trans, post, genotypeio, susie
import sys

plink_prefix_path = sys.argv[1]
expression_bed = sys.argv[2]
covariates_file = sys.argv[3]
prefix = sys.argv[4]

phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0)
pgr = pgen.PgenReader(plink_prefix_path,  select_samples=phenotype_df.columns)
genotype_df = pgr.load_dosages()
variant_df = pgr.variant_df
variant_df["chrom"] = "chr"+variant_df["chrom"]

cis.map_nominal(genotype_df, pgr.variant_df, phenotype_df, phenotype_pos_df, prefix ,covariates_df=covariates_df,write_top=True, write_stats=True)
pairs_df = pd.DataFrame()
for i in range(1,23):
    chr = "chr" + str(i)
    pairs_df2 = pd.read_parquet(f'{prefix}.cis_qtl_pairs.{chr}.parquet')
    pairs_df = pd.concat([pairs_df,pairs_df2])

cis_df = cis.map_cis(genotype_df, pgr.variant_df, phenotype_df, phenotype_pos_df, covariates_df=covariates_df, seed=42, warn_monomorphic=False)
post.calculate_qvalues(cis_df, fdr=0.05, qvalue_lambda=0.85)
pval_dict = cis_df['pval_nominal_threshold'].to_dict()
qval_dict = cis_df['qval'].to_dict()
susie_df = susie.map(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df)
susie_df.to_csv(prefix + ".susie.eqtl.txt",sep="\t")
indep_df = cis.map_independent(genotype_df, pgr.variant_df, cis_df, phenotype_df, phenotype_pos_df, covariates_df)
indep_df.to_csv(prefix + ".cond.eqtl.txt",index=False,sep="\t")

pairs_df['pval_nominal_threshold'] = pairs_df["phenotype_id"].map(pval_dict)
pairs_df['qval'] = pairs_df["phenotype_id"].map(qval_dict)
sig_pairs_df = pairs_df[(pairs_df['pval_nominal'] <= pairs_df['pval_nominal_threshold']) & (pairs_df['qval']<=0.05)]
sig_pairs_df.to_csv(prefix + ".eqtl.txt",index=False,sep="\t")
sig_top_df = sig_pairs_df.sort_values(by=['pval_nominal'], ascending=[True]).groupby('phenotype_id').first()
sig_top_df.to_csv(prefix + ".top.eqtl.txt",sep="\t")