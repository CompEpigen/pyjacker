import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
import os
import sys
import pickle
from scipy.stats import fisher_exact

import config
config.load_config(sys.argv[1])
from breakpoints import find_distance_breakpoint, find_fusion, find_enhancers, get_CN_gene
from outlier_expression import nbinom_pvals_gene, nbinom_pvals_gene_sample, OHE_pvals_gene
from monoallelic_expression import find_SNPs_DNA_RNA_sample_gene,compute_monoallelic_llrs_sample_gene
from html_reports import generate_main_report, generate_individual_reports




def find_enhancer_hijacking():
    d = {"gene_id":[],"gene_name":[],"score":[],"chr":[],"start":[],"end":[],"sample":[],"pval":[],"n_SNPs":[],"pval_monoallelic":[],"llr_monoallelic":[],"distance_to_breakpoint":[],"fusion":[],"enhancers":[]}
    for gene_id in config.df_TPM.index:
        gene = config.data.gene_by_id(gene_id)
        gene_name = gene.gene_name

         # Exclude genes with low expression for all samples, non-protein coding genes, genes on Y chromosome and in the mitochondrial genome
        if np.sum(config.df_TPM.loc[gene_id,:])<5: continue
        if gene.biotype in ["snRNA","pseudogene","snoRNA","antisense","lincRNA"]: continue
        if gene.contig in ["Y","MT"]: continue
        
        pvals = OHE_pvals_gene(gene_id)
        if gene_name=="EBLN1": print(pvals)
        for sample in pvals:
            print((sample,gene_name))
            d["gene_id"].append(gene_id)
            d["gene_name"].append(gene_name)
            d["chr"].append(gene.contig)
            d["start"].append(gene.start)
            d["end"].append(gene.end)
            d["sample"].append(sample)
            d["pval"].append(pvals[sample])
            #monoallelic,n_SNPs,pval_monoallelic=monoallelic_expression(sample,gene)
            #d["monoallelic"].append(monoallelic)
            if config.SNV_RNA_dir is not None:
                d["n_SNPs"].append(len(config.SNPs[sample][gene_id]))
                d["pval_monoallelic"].append(config.pvalues_monoallelic[sample][gene_id])
                d["llr_monoallelic"].append(config.llrs_monoallelic[sample][gene_id])
            else:
                if pvals[sample]<0.05:
                    SNPs = find_SNPs_DNA_RNA_sample_gene(sample,gene)
                    d["n_SNPs"].append(len(SNPs))
                    llr = compute_monoallelic_llrs_sample_gene(SNPs)
                    d["llr_monoallelic"].append(llr)
                    config.llrs_monoallelic[sample][gene_id]=llr
                else:
                    d["n_SNPs"].append(0)
                    d["llr_monoallelic"].append(0)
                if gene_name=="MNX1":print(SNPs)
                d["pval_monoallelic"].append(1)
            d["distance_to_breakpoint"].append(find_distance_breakpoint(sample,gene))
            d["fusion"].append(find_fusion(sample,gene))
            d["enhancers"].append(str(find_enhancers(sample,gene)))
            score = -5*np.log(pvals[sample])
            llr = config.llrs_monoallelic[sample][gene_id]
            if llr<=0:
                score+= llr - 10
            elif llr>=10:
                score+=10+(llr-10)*0.1
            else:
                score+=llr
            d["score"].append(score)

    df_result = pd.DataFrame(d)

    df_result = df_result.sort_values("score",ascending=False)
    print(df_result)

    df_result["padj"] = multipletests(df_result["pval"], alpha=0.05, method='fdr_bh',is_sorted=False)[1]
    df_result = df_result.loc[df_result["pval"]<0.1,:]
    df_result = df_result.reset_index(drop=True)
    
    df_result.to_csv(os.path.join(config.output_dir,"enhancer_hijacking.tsv"),index=False,sep="\t")
    return df_result



df_result = find_enhancer_hijacking()


generate_main_report(df_result,config.output_dir,100,filter_monoallelic=False)
generate_individual_reports(df_result,config.output_dir,n_events=200)

