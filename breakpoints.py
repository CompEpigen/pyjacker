import gzip
import numpy as np
import os
import pandas as pd
import vcfpy
from collections import namedtuple
Breakpoint = namedtuple('Breakpoint', 'sample chr pos chr2 pos2')
CopyNumberVariant = namedtuple('CopyNumberVariant', 'chr start end cn')

import config



def read_breakpoints(df_breakpoints):
    breakpoints = {}
    for chr in config.chromosomes: breakpoints[chr] = []

    for i in df_breakpoints.index:
        # Filter out small SVs ??
        sample = df_breakpoints.loc[i,"sample"]
        chr1 = str(df_breakpoints.loc[i,"chr1"])
        pos1 = df_breakpoints.loc[i,"pos1"]
        if "chr2" in df_breakpoints.columns:
            chr2 = str(df_breakpoints.loc[i,"chr2"])
            pos2 = df_breakpoints.loc[i,"pos2"]
        else:
            chr2="NA"
            pos2="NA"
        
        if sample in config.samples:
            if chr1 in config.chromosomes and not Breakpoint(sample,chr1,pos1,chr2,pos2) in breakpoints[chr1]:
                breakpoints[chr1].append(Breakpoint(sample,chr1,pos1,chr2,pos2))
            if chr2 in config.chromosomes and not Breakpoint(sample,chr2,pos2,chr1,pos1) in breakpoints[chr2]:
                breakpoints[chr2].append(Breakpoint(sample,chr2,pos2,chr1,pos1))
    #for chr in config.chromosomes:
    #    breakpoints[chr] = sorted(breakpoints[chr])
    return breakpoints

def read_CNAs(df_CNA):
    CNAs={}
    for sample in config.samples:
        CNAs[sample] = {}
        for chr in config.chromosomes:
            CNAs[sample][chr]=[]
    for i in df_CNA.index:
        sample = df_CNA.loc[i,"sample"]
        chr = df_CNA.loc[i,"chr"]
        if (not chr in config.chromosomes) or not sample in config.samples: continue
        start = df_CNA.loc[i,"start"]
        end = df_CNA.loc[i,"end"]
        cn = df_CNA.loc[i,"cn"]
        CNAs[sample][chr].append(CopyNumberVariant(chr,start,end,cn))
    return CNAs


def find_TAD_boundary(chr,pos):
    dist_beyond_boundary = 100000 # 300000
    index_TAD=0
    while index_TAD +1 < len(config.TADs[chr]) and pos - dist_beyond_boundary > config.TADs[chr][index_TAD][1]:
        index_TAD+=1
    left_boundary = config.TADs[chr][index_TAD][0] - dist_beyond_boundary
    while index_TAD +2 < len(config.TADs[chr]) and pos+dist_beyond_boundary > config.TADs[chr][index_TAD+1][0]:
        index_TAD+=1
    right_boundary = config.TADs[chr][index_TAD][1] + dist_beyond_boundary
    return (left_boundary,right_boundary)

def find_samples_with_SV_in_same_TAD(gene,breakpoints,exclude_deletion=True,exclude_amplification=False):
    """Find samples which have breakpoints in the TAD where a gene is located. Also look 100kb beyond the TAD boundaries, since they may not be very precise. """
    dist_beyond_boundary = 100000 # 300000
    chr = gene.contig
    start = gene.start
    end = gene.end
    if gene.strand=="+":
        TSS_pos = start
    else:
        TSS_pos = end
    if not chr in config.chromosomes:
        return set()

    #index_TAD=0
    #while index_TAD +1 < len(config.TADs[chr]) and start > config.TADs[chr][index_TAD][1]:
    #    index_TAD+=1
    #left_boundary = config.TADs[chr][index_TAD][0] - dist_beyond_boundary
    #while index_TAD +2 < len(config.TADs[chr]) and end > config.TADs[chr][index_TAD+1][0]:
    #    index_TAD+=1
    #right_boundary = config.TADs[chr][index_TAD][1] + dist_beyond_boundary
    left_boundary,right_boundary = find_TAD_boundary(chr,TSS_pos)
    samples_SV = set()

    for breakpoint in breakpoints[chr]:
        if breakpoint.pos > left_boundary and breakpoint.pos < right_boundary:
            if (not exclude_deletion) or (not gene_is_deleted(breakpoint.sample,gene)): #or chr=="X"
                if (not exclude_amplification) or (not gene_is_amplified(breakpoint.sample,gene)):
                    samples_SV.add(breakpoint.sample)
    return samples_SV

def find_distance_breakpoint(sample,gene):
    chr = gene.contig
    start = gene.start
    end = gene.end
    min_dist=100000000
    for breakpoint in config.breakpoints[chr] :
        if breakpoint.sample == sample:
            if breakpoint.pos>=start and breakpoint.pos<=end: min_dist=0
            min_dist = min(min_dist,abs(breakpoint.pos-start))
            min_dist = min(min_dist,abs(breakpoint.pos-end))
    return min_dist





def find_fusion(sample,gene):
    if not "FusionTranscripts_file" in dir(config): 
        # Find gene fusions based on WGS, not RNAseq
        gene_partners=set()
        for breakpoint in config.breakpoints[gene.contig]:
            if breakpoint.pos2 == breakpoint.pos2 and breakpoint.pos2!="NA":
                if breakpoint.sample == sample and breakpoint.pos >=gene.start and breakpoint.pos <= gene.end and breakpoint.pos2==breakpoint.pos2:
                    for gene2 in config.data.genes_at_locus(breakpoint.chr2,int(breakpoint.pos2),int(breakpoint.pos2+1)):
                        gene_partners.add(gene2.gene_name)
        return ",".join(list(gene_partners))
    else:
        # Use gene fusions detected in RNAseq data
        gene_name = gene.gene_name
        fusion=""
        df_fusions = pd.read_csv(config.FusionTranscripts_file,sep="\t")
        for i in range(df_fusions.shape[0]):
            if df_fusions.loc[i,"sample"]==sample and (df_fusions.loc[i,"LeftGene"]==gene_name or df_fusions.loc[i,"RightGene"]==gene_name):
                fusion= df_fusions.loc[i,"LeftGene"]+"--"+df_fusions.loc[i,"RightGene"]
        return fusion


def get_CN_gene(sample,gene):
    chr = gene.contig
    #pos = (gene.start + gene.end) // 2
    cn=2
    if "df_CNAs" in dir(config):
        for CNA in config.CNAs[sample][chr]:
            if gene.start > CNA.start and gene.end < CNA.end:
                cn = CNA.cn
    return cn

def gene_is_deleted(sample,gene):
    # Here, only consider deletions where the whole gene is deleted. (In 15PB19457, only part of MECOM is deleted, but not EVI1).
    chr = gene.contig
    start = gene.start
    end = gene.end
    if "df_CNAs" in dir(config):
        for CNA in config.CNAs[sample][chr]:
            if start > CNA.start and end < CNA.end and CNA.cn <=1:
                return True
    return False

def gene_is_amplified(sample,gene):
    # More than 8 copies
    chr = gene.contig
    start = gene.start
    end = gene.end
    pos = (start+end)//2
    if "df_CNAs" in dir(config):
        for CNA in config.CNAs[sample][chr]:
            if pos > CNA.start and pos < CNA.end and CNA.cn >=8:
                return True
    return False

def correct_counts_CN(df):
    df = df.copy(deep=True)
    c=0
    for gene_id in df.index:
        c+=1
        if c%100==0: print(c)
        gene = config.data.gene_by_id(gene_id)
        for sample in df.columns:
            #print((sample,gene))
            df.loc[gene_id,sample] = df.loc[gene_id,sample] * 2 / max(2,get_CN_gene(sample,gene)) 
    return df

def correct_TPM_CN(df):
    df = df.copy(deep=True)
    c=0
    for gene_id in df.index:
        c+=1
        if c%100==0: print(c)
        gene = config.data.gene_by_id(gene_id)
        for sample in df.columns:
            #print((sample,gene))
            df.loc[gene_id,sample] = df.loc[gene_id,sample] * 2 / max(2,get_CN_gene(sample,gene)) 
    for c in df.columns:
        df[c] = df[c]*1000000 / np.sum(df[c])
    return df


def find_enhancers(sample,gene):
    """Look for enhancers on the other side of a breakpoint"""
    chr = gene.contig
    start = gene.start
    end = gene.end
    if gene.strand=="+":
        TSS_pos = start
    else:
        TSS_pos = end
    left_boundary,right_boundary = find_TAD_boundary(chr,TSS_pos)
    enhancers = []
    for bp in config.breakpoints[chr]:
        if bp.pos>=left_boundary and bp.pos<=right_boundary and bp.sample==sample and bp.chr in config.chromosomes and bp.chr2 in config.chromosomes:
            if bp.chr2!="NA" and bp.chr2==bp.chr2 and bp.chr2!="nan":
                left_boundary2,right_boundary2 = find_TAD_boundary(bp.chr2,bp.pos2)
                for i in config.df_enhancers.index:
                    if config.df_enhancers.loc[i,"chr"]==bp.chr2 and  config.df_enhancers.loc[i,"start"]>=left_boundary2 and config.df_enhancers.loc[i,"end"]<=right_boundary2:
                        reg = (bp.chr2,config.df_enhancers.loc[i,"start"],config.df_enhancers.loc[i,"end"])
                        if not reg in enhancers:
                            enhancers.append(reg)
    return enhancers
