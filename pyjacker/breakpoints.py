import gzip
import numpy as np
import os
import pandas as pd
import vcfpy
from collections import namedtuple
Breakpoint = namedtuple('Breakpoint', 'sample chr1 pos1 orientation1 chr2 pos2 orientation2')

#import config
from pyjacker.genes import genes_at_locus


def read_breakpoints(df_breakpoints):
    breakpoints = {}
    for i in df_breakpoints.index:
        # Filter out small SVs ??
        sample = df_breakpoints.loc[i,"sample"]
        chr1 = str(df_breakpoints.loc[i,"chr1"]).strip("chr")
        pos1 = df_breakpoints.loc[i,"pos1"]
        if "chr2" in df_breakpoints.columns:
            chr2 = str(df_breakpoints.loc[i,"chr2"]).lstrip("chr")
            pos2 = df_breakpoints.loc[i,"pos2"]
        else:
            chr2="NA"
            pos2="NA"
        if "orientation1" in df_breakpoints.columns:
            orientation1 = df_breakpoints.loc[i,"orientation1"]
            orientation2 = df_breakpoints.loc[i,"orientation2"]

        #if chr1==chr2 and abs(pos1-pos2)<5000: continue # Filter out small SVs ?

        bp1=Breakpoint(sample,chr1,pos1,orientation1,chr2,pos2,orientation2)
        bp2=Breakpoint(sample,chr2,pos2,orientation2,chr1,pos1,orientation1)
        
        if not chr1 in breakpoints: breakpoints[chr1] = []
        if not bp1 in breakpoints[chr1]: breakpoints[chr1].append(bp1)
        if not chr2 in breakpoints: breakpoints[chr2] = []
        if not bp2 in breakpoints[chr2]: breakpoints[chr2].append(bp2)
    return breakpoints

def find_samples_with_breakpoints_near_gene(samples,gene,breakpoints,TADs=None,max_dist_bp2tss=1500000):
    if not gene.chr in breakpoints: return set()
    if gene.strand=="+": TSS_pos = gene.start
    else: TSS_pos = gene.end
    if TADs is not None:
        left_boundary,right_boundary = find_TAD_boundary(gene.chr,TSS_pos,TADs,max_dist_bp2tss)
    else:
        left_boundary,right_boundary = TSS_pos - max_dist_bp2tss , TSS_pos + max_dist_bp2tss

    samples_SV = set()
    for bp in breakpoints[gene.chr]:
        if bp.pos1 > left_boundary and bp.pos1 < right_boundary:
            #if (not exclude_deletion) or (not gene_is_deleted(breakpoint.sample,gene)): #or chr=="X"
            #    if (not exclude_amplification) or (not gene_is_amplified(breakpoint.sample,gene)):
            if bp.sample in samples:
                samples_SV.add(bp.sample)

    return samples_SV

def read_TADs(TADs_file):
    TADs={}
    with open(TADs_file,"r") as infile:
        for line in infile:
            linesplit = line.rstrip("\n").split("\t")
            chr = linesplit[0].lstrip("chr")
            if not chr in TADs:
                TADs[chr] = []
            TADs[chr].append((int(linesplit[1]),int(linesplit[2])))
    return TADs

def find_TAD_boundary(chr,pos,TADs,max_dist_bp2tss=1500000):
    """For a given position, find the boundaries of its TAD. If close to a TAD border, also include the next TAD."""
    dist_beyond_boundary = 80000 
    index_TAD=0
    while index_TAD +1 < len(TADs[chr]) and pos - dist_beyond_boundary > TADs[chr][index_TAD][1]:
        index_TAD+=1
    left_boundary = TADs[chr][index_TAD][0] - dist_beyond_boundary

    # In case the position is not located in any TAD.
    if pos + dist_beyond_boundary < TADs[chr][index_TAD][0]: return (pos-max_dist_bp2tss,pos+max_dist_bp2tss)

    while index_TAD +2 < len(TADs[chr]) and pos+dist_beyond_boundary > TADs[chr][index_TAD+1][0]:
        index_TAD+=1
    right_boundary = TADs[chr][index_TAD][1] + dist_beyond_boundary
    return (left_boundary,right_boundary)


def find_distance_breakpoint(breakpoints,sample,gene):
    chr = gene.chr
    start = gene.start
    end = gene.end
    min_dist=100000000
    if chr in breakpoints:
        for breakpoint in breakpoints[chr] :
            if breakpoint.sample == sample:
                if breakpoint.pos1>=start and breakpoint.pos1<=end: min_dist=0
                min_dist = min(min_dist,abs(breakpoint.pos1-start))
                min_dist = min(min_dist,abs(breakpoint.pos1-end))
    return min_dist



def find_fusion(breakpoints,genes_index,sample,gene,df_fusions=None):
    if df_fusions is None: 
        # Find gene fusions based on WGS, not RNAseq
        gene_partners=set()
        if gene.chr in breakpoints:
            for breakpoint in breakpoints[gene.chr]:
                if breakpoint.pos2 == breakpoint.pos2 and breakpoint.pos2!="NA":
                    if breakpoint.sample == sample and breakpoint.pos1 >=gene.start and breakpoint.pos1 <= gene.end and breakpoint.pos2==breakpoint.pos2:
                        for gene2 in genes_at_locus(genes_index,breakpoint.chr2,int(breakpoint.pos2)):
                            gene_partners.add(gene2.gene_name)
        return ", ".join(list(gene_partners))
    else:
        # Use gene fusions detected in RNAseq data
        gene_name = gene.gene_name
        fusions=[]
        for i in range(df_fusions.shape[0]):
            if df_fusions.loc[i,"sample"]==sample and (df_fusions.loc[i,"LeftGene"]==gene_name or df_fusions.loc[i,"RightGene"]==gene_name):
                fusions.append(df_fusions.loc[i,"LeftGene"]+"--"+df_fusions.loc[i,"RightGene"])
        return ",".join(fusions)




def find_enhancers(breakpoints,TADs,df_enhancers,sample,gene):
    """Look for enhancers on the other side of a breakpoint"""
    chr = gene.chr
    start = gene.start
    end = gene.end
    if gene.strand=="+":
        TSS_pos = start
    else:
        TSS_pos = end
    left_boundary,right_boundary = find_TAD_boundary(chr,TSS_pos,TADs)
    enhancers = []
    for bp in breakpoints[chr]:
        if bp.pos1>=left_boundary and bp.pos1<=right_boundary and bp.sample==sample:
            if bp.chr2!="NA" and bp.chr2==bp.chr2 and bp.chr2!="nan":
                left_boundary2,right_boundary2 = find_TAD_boundary(bp.chr2,bp.pos2,TADs)
                for i in df_enhancers.index:
                    if df_enhancers.loc[i,"chr"]==bp.chr2 and df_enhancers.loc[i,"start"]>=left_boundary2 and df_enhancers.loc[i,"end"]<=right_boundary2:
                        reg = (bp.chr2,df_enhancers.loc[i,"start"],df_enhancers.loc[i,"end"])
                        if not reg in enhancers:
                            enhancers.append(reg)
    return enhancers



def find_enhancers_orientation(breakpoints,TADs,df_enhancers,sample,gene):
    """Look for enhancers on the other side of a breakpoint"""
    chr = gene.chr
    start = gene.start
    end = gene.end
    if gene.strand=="+":
        TSS_pos = start
    else:
        TSS_pos = end
    left_boundary,right_boundary = find_TAD_boundary(chr,TSS_pos,TADs)
    breakpoints_left=[]
    breakpoints_right=[]
    minpos_right=1e10
    maxpos_left=0
    for bp in breakpoints[chr]:
        if bp.pos1>=left_boundary and bp.pos1<=right_boundary and bp.sample==sample:
            if bp.pos1>start:
                minpos_right= min(minpos_right,bp.pos1)
                if bp.orientation1=="-":
                    breakpoints_right.append(bp)
            elif bp.pos1<end:
                maxpos_left= max(maxpos_left,bp.pos1)
                if bp.orientation1=="+":
                    breakpoints_left.append(bp)

    # Only keep breakpoints which are close to the first breakpoint encountered.
    breakpoints_filtered = []
    for bp in breakpoints_left:
        if abs(bp.pos1-maxpos_left)<=20000: breakpoints_filtered.append(bp)
    for bp in breakpoints_right:
        if abs(bp.pos1-minpos_right)<=20000: breakpoints_filtered.append(bp)

    enhancers = []
    super_enhancers = []
    scores=[]
    for bp in breakpoints_filtered:
        left_boundary2,right_boundary2 = find_TAD_boundary(bp.chr2,bp.pos2,TADs)
        if bp.orientation2=="+": left_boundary2=bp.pos2
        else: right_boundary2= bp.pos2
        for i in df_enhancers.index:
            if df_enhancers.loc[i,"CHROM"]==bp.chr2 and  df_enhancers.loc[i,"STOP"]>=left_boundary2 and df_enhancers.loc[i,"START"]<=right_boundary2:
                reg = bp.chr2+":"+str(df_enhancers.loc[i,"START"])+"-"+str(df_enhancers.loc[i,"STOP"])
                if not reg in enhancers:
                    enhancers.append(reg)
                    if "avg" in df_enhancers.columns:
                        scores.append(df_enhancers.loc[i,"avg"]/10000)
                    else:
                        if "isSuper" in df_enhancers.columns and df_enhancers.loc[i,"isSuper"]:
                            scores.append(2)
                        else: scores.append(0.5)
                    if "isSuper" in df_enhancers.columns and df_enhancers.loc[i,"isSuper"]: super_enhancers.append(reg)
    scores=sorted(scores,reverse=True)
    score=0
    for i in range(len(scores)):
        score+=scores[i] / (i+1)

    return ",".join(enhancers), ",".join(super_enhancers), score