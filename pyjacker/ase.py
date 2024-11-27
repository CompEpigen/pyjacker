import os
import numpy as np
import pandas as pd
from scipy.stats import betabinom
import vcfpy
from tqdm import tqdm
import sys

#from collections import namedtuple
#SNP = namedtuple('SNP', 'id chr pos ref alt RO AD')



from pyjacker.genes import index_genes_by_pos,genes_at_locus


def compute_ase_matrix(samples,ase_dir,genes,genes_index=None,prior_coef=2.0,imprinted_genes_file=None,CNAs=None):
    """
    Inputs:
        samples: list of sample names
        ase_dir: directory containing the ASEReadCounter output
        genes: directory of gene IDs to genes (gene_id,gene_name,chr,start,end,strand)
        imprinted_genes_file: text file where each line is the name of an imprinted gene. Imprinted genes will be ignored.
    Output:
        dataframe of gene_ids by samples giving the ASE score of each pair (gene_id,sample).
    """
    if genes_index is None: genes_index = index_genes_by_pos(genes)
    imprinted_genes = []
    if imprinted_genes_file is not None: 
        with open(imprinted_genes_file,"r") as infile:
            for line in infile:
                imprinted_genes.append(line.rstrip("\n"))
    
    df = pd.DataFrame(np.zeros((len(genes),len(samples))))
    df.index = [genes[g].gene_id for g in genes]
    df.columns = samples
    if ase_dir is None: return df
    for sample in tqdm(samples,file=sys.stdout):
        if CNAs is not None and sample in CNAs: CNAs_sample = CNAs[sample]
        else: CNAs_sample = {}
        ase_file=os.path.join(ase_dir,sample+".tsv")
        geneIDs2llrs={}
        df_sample = pd.read_csv(ase_file,sep="\t",dtype={"contig":str})
        df_sample["contig"] = [x.lstrip("chr") for x in df_sample["contig"]]
        for x in df_sample.index:

            # Exclude positions which have copy number <2 or >=5, because they are expected to deviate from biallelic expression.
            chr_snp = df_sample.loc[x,"contig"]
            pos_snp = df_sample.loc[x,"position"]
            snp_near_diploid=True
            if chr_snp in CNAs_sample:
                for (start,end,cn) in CNAs_sample[chr_snp]:
                    if start<=pos_snp and pos_snp<=end and (cn<=2 or cn>=5): snp_near_diploid = False #if cn==2, must be because ploidy was higher.
            if not snp_near_diploid: continue

            llr = llr_betabinom(df_sample.loc[x,"altCount"], df_sample.loc[x,"altCount"]+df_sample.loc[x,"refCount"])
            for g in genes_at_locus(genes_index,df_sample.loc[x,"contig"].lstrip("chr"),int(df_sample.loc[x,"position"])):
                if not g.gene_id in genes: continue
                if not g.gene_id in geneIDs2llrs: geneIDs2llrs[g.gene_id]=[]
                geneIDs2llrs[g.gene_id].append(llr)
        for gene_id in geneIDs2llrs:
            if (not genes[gene_id].gene_name in imprinted_genes) and (not genes[gene_id].gene_name in ["KCNJ12"]): 
                if genes[gene_id].chr!="Y" and (genes[gene_id].chr!="X" or genes[gene_id].start<2700000): # Exclude chrY and chrX, except the PAR.
                    df.loc[gene_id,sample] = compute_ase_score_from_llrs(geneIDs2llrs[gene_id])
    return df


def llr_betabinom(k,n,alpha=10,beta=10):
    """Compute log likelihood ratio between monoallelic and biallelic and expression."""
    loglik_biallelic = np.log(betabinom.pmf(k,n,alpha,beta)*0.9+0.1/n)

    lik_monoallelic1 = betabinom.pmf(k,n,1,49)
    lik_monoallelic2 = betabinom.pmf(k,n,49,1)
    loglik_monoallelic = np.log(0.499 * lik_monoallelic1 + 0.499 * lik_monoallelic2+0.002/n)
    return loglik_monoallelic - loglik_biallelic

def compute_ase_score_from_llrs(llrs,prior_coef=2.0):
    return np.sum(llrs) / (len(llrs) + prior_coef)


def count_SNPs_gene_sample(ase_dir,sample,gene=None):
    if ase_dir is None: return 0
    chr = gene.chr
    start = gene.start
    end = gene.end
    df = pd.read_csv(os.path.join(ase_dir,sample+".tsv"),sep="\t",dtype={"contig":str})
    df["contig"] = [x.lstrip("chr") for x in df["contig"]]
    df = df.loc[df["contig"]==chr]
    df = df.loc[(df["position"]>=start) & (df["position"]<=end)]
    return df.shape[0]
