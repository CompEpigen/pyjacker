import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
import sys

from pyjacker.genes import read_genes_gtf

def read_CNAs(df_CNA,samples):
    CNAs={sample:{} for sample in samples}
    for i in df_CNA.index:
        sample = df_CNA.loc[i,"sample"]
        if not sample in samples: continue
        chr = df_CNA.loc[i,"chr"].lstrip("chr")
        #if (not chr in config.chromosomes) or not sample in samples: continue
        start = df_CNA.loc[i,"start"]
        end = df_CNA.loc[i,"end"]
        cn = df_CNA.loc[i,"cn"]
        if not chr in CNAs[sample]: CNAs[sample][chr]=[]
        CNAs[sample][chr].append((start,end,cn))
    return CNAs

def get_CN_gene(CNAs,sample,gene):
    chr = gene.chr
    cn_gene=2
    if sample in CNAs:
        if chr in CNAs[sample]:
            for (start,end,cn) in CNAs[sample][chr]:
                if gene.start > start and gene.end < end:
                    cn_gene = cn
    return cn_gene

def find_amplified_genes(CNAs,genes,sample,gene_set=None):
    """
    Returns a dictionary where the keys are the ids of the amplified genes in this sample and the values are the corresponding copy numbers.
    """
    amplified_genes = {}
    for x in genes:
        if gene_set is not None and (not genes[x].gene_id in gene_set): continue # Ignore genes for which we do not have expression anyway
        chr = genes[x].chr
        start = genes[x].start
        end = genes[x].end
        if sample in CNAs:
            if chr in CNAs[sample]:
                for (CNA_start,CNA_end,cn) in CNAs[sample][chr]:
                    if cn>2 and CNA_start<end and CNA_end > start: amplified_genes[genes[x].gene_id] = cn
    return amplified_genes

def compute_ploidy(sample,CNAs,chr_lengths):
    total_length=0
    sum_ploidy=0
    for chr in chr_lengths:
        length_diploid_chr = chr_lengths[chr]
        if sample in CNAs and chr in CNAs[sample]:
            for (start,end,cn) in CNAs[sample][chr]:
                sum_ploidy+= cn * (end-start) / 1000
                length_diploid_chr-=(end-start)
        sum_ploidy+=2*length_diploid_chr/1000
        total_length+=chr_lengths[chr] / 1000

    ploidy = sum_ploidy / total_length
    return ploidy
        

def correct_exp_cn(df_TPM,CNAs,genes=None,genome=None,gtf=None,chr_lengths=None):
    """
    When the copy number of a gene is >2, divide the expression by (copy number / 2).
    This (partially) filters out high expression simply due to amplification, as opposed to enhancer hijacking.
    This function needs a mapping from gene_id to coordinates. Either this mapping has already been computed,
    or it can be computed from a gtf file or using pyensembl, which requires to provide a genome version.
    """
    if genes is None:
        if gtf is not None:
            genes,_ = read_genes_gtf(gtf)
        elif genome is not None:
            pass
            #genes = read_genes_pyensembl(genome)
        else:
            sys.exit("Must provide a gtf file or a genome version.")
    
    set_genes_TPM = set(df_TPM.index)
    
    for sample in tqdm(df_TPM.columns,file=sys.stdout):
        if chr_lengths is not None: ploidy = compute_ploidy(sample,CNAs,chr_lengths)
        else: ploidy=2
        amplified_genes = find_amplified_genes(CNAs,genes,sample,set_genes_TPM)
        for gene_id in amplified_genes:
            #if not gene_id in set_genes_TPM: continue
            df_TPM.loc[gene_id,sample] =  df_TPM.loc[gene_id,sample] / (amplified_genes[gene_id]/ploidy)
    return df_TPM

def add_cna_breakpoints(CNAs,df_breakpoints = None,chr_lengths=None):
    """
    Identify breakpoints using copy number switches.
    This enables the use of pyjacker without WGS.
    Even if breakpoints from WGS are available, this can add breakpoints which were missed by the SV caller, possibly because of low mappability.
    If breakpoints are already provided, this will only add new breakpoints not already present.  
    Providing chromosome lengths will prevent adding breakpoints at the end of chromosomes.  
    """

    d={"sample":[],"chr1":[],"pos1":[],"orientation1":[],"chr2":[],"pos2":[],"orientation2":[]}
    for sample in CNAs:
        for chr in CNAs[sample]:
            for (start,end,cn) in CNAs[sample][chr]:
                # Only add breakpoints if they are not at the start or end of the chromosome.
                bps = []
                if start>10: bps.append((sample,chr,start))
                if chr_lengths is None or abs(chr_lengths[chr]-end)>10000: bps.append((sample,chr,end)) 
                
                # check that there is not already a breakpoint close by.
                for bp in bps:
                    bp_new = True
                    if df_breakpoints is not None:
                        df_breakpoints_sample = df_breakpoints.loc[df_breakpoints["sample"]==bp[0],:]
                        for x in df_breakpoints_sample.index:
                            if bp[1]==df_breakpoints_sample.loc[x,"chr1"] and abs(bp[2]-df_breakpoints_sample.loc[x,"pos1"])<50000:
                                bp_new=False
                                break
                    for i in range(len(d["sample"])):
                        if d["sample"][i]==bp[0] and d["chr1"][i]==bp[1] and abs(d["pos1"][i]-bp[2])<50000:
                            bp_new=False
                            break

                    if bp_new: # if the breakpoints is new
                        d["sample"].append(bp[0])
                        d["chr1"].append(bp[1])
                        d["pos1"].append(bp[2])
                        d["orientation1"].append("NA")
                        d["chr2"].append("NA")
                        d["pos2"].append("NA")
                        d["orientation2"].append("NA")
    df_breakpoints_cna = pd.DataFrame(d)
    if df_breakpoints is not None:
        df_breakpoints = pd.concat([df_breakpoints,df_breakpoints_cna])
        df_breakpoints.reset_index(inplace=True)
    else:
        df_breakpoints = df_breakpoints_cna
    return df_breakpoints

def compute_penalty_amplification(CNAs,sample,gene):
    if CNAs is None: return 0
    cn = get_CN_gene(CNAs,sample,gene)
    if cn>2:return  (cn-2)
    else: return 0
    
def compute_penalty_deletion(CNAs,sample,gene,chr_lengths):
    chr = gene.chr
    if sample in CNAs and chr in CNAs[sample]:
        for i in range(len(CNAs[sample][chr])):
            penalty=0
            if CNAs[sample][chr][i][1] > gene.start:
                # Get the copy number of the gene
                if CNAs[sample][chr][i][0]<gene.start: cn_gene= max(1,CNAs[sample][chr][i][2])
                else: cn_gene=2

                if cn_gene<2: return 3

                # Penalize if the gene has a lower copy number than the segment to its left.
                if i>0 and abs(CNAs[sample][chr][i-1][1]-gene.start)<2000000: # There is a CNA before the gene
                    if CNAs[sample][chr][i-1][2]>cn_gene: penalty += 0.5+CNAs[sample][chr][i-1][2]/cn_gene
                elif i==0 and abs(CNAs[sample][chr][i][1]-gene.start) < 2000000 and CNAs[sample][chr][i][0]>100: # There is no CNA before, but the breakpoint is close to the gene.
                    if 2>cn_gene: penalty += 2.5

                # Penalize if the gene has a lower copy number than the segment to its right.
                if i<len(CNAs[sample][chr])-1 and abs(CNAs[sample][chr][i+1][1]-gene.start)<3000000: # There is a CNA after the gene
                    if CNAs[sample][chr][i+1][2]>cn_gene: penalty += 0.5+CNAs[sample][chr][i+1][2]/cn_gene
                elif i==len(CNAs[sample][chr])-1 and abs(CNAs[sample][chr][i][1]-gene.end) < 3000000 and \
                    chr_lengths is not None and chr in chr_lengths and CNAs[sample][chr][i][1] < chr_lengths[chr]-1000: # There is no CNA after, but the breakpoint is close to the gene.
                    if 2>cn_gene: penalty += 2.5
                if penalty>0: print((sample,gene.gene_name,penalty))
                return penalty
    return 0

def gene_is_deleted(CNAs,sample,gene):
    if CNAs is None: return 0
    cn = get_CN_gene(CNAs,sample,gene)
    return cn<2