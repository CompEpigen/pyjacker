import gzip
import numpy as np
import os
import pandas as pd
import pickle
import pyensembl
import sys
import yaml

from breakpoints import read_breakpoints, read_CNAs, correct_TPM_CN
from monoallelic_expression import read_SNPs_RNA, find_SNPs_DNA_RNA, compute_monoallelic_pvalues, compute_monoallelic_llrs

data = pyensembl.EnsemblRelease(75) 

def load_config(config_file):
    print("Loading config from "+ config_file)
    global chromosomes,samples, TADs_file, max_dist_bp2tss, TADs, FusionTranscripts_file, SNV_RNA_dir, BAM_DNA_template, BAM_RNA_template, df_TPM, df_TPM_other, df_TPM_normal, df_TPM_reference
    global df_CNAs, df_breakpoints, CNAs, breakpoints, df_enhancers, output_dir, cache_dir, SNPs, pvalues_monoallelic, llrs_monoallelic
    with open(config_file, 'r') as stream:
        data_yaml = yaml.safe_load(stream)

        # General data
        #if "chr_lengths_file" in data_yaml:
        #    chr_lengths_file = data_yaml["chr_lengths_file"]
        #else:
        #    sys.exit("Missing chr_lengths_file in the config file.")
        if "TADs_file" in data_yaml:
            TADs_file = data_yaml["TADs_file"]
        else:
            if "max_dist_bp2tss" in data_yaml:
                max_dist_bp2tss= data_yaml["max_dist_bp2tss"]
            else:
                max_dist_bp2tss=1500000
            print("No TAD file was provided the config file; will look for SVs located within "+str(max_dist_bp2tss)+"bp of the TSS.")

        chromosomes = [str(x) for x in range(1,23)] + ["X"]
        #chr_lengths={}
        #with gzip.open(chr_lengths_file, 'rt') as infile:
        #    for line in infile:
        #        linesplit = str(line).rstrip("\n").split("\t")
        #        chr_lengths[linesplit[0].lstrip("chr")] = int(linesplit[2])

        TADs={}
        with open(TADs_file,"r") as infile:
            for line in infile:
                linesplit = line.rstrip("\n").split("\t")
                chr = linesplit[0].lstrip("chr")
                if not chr in TADs:
                    TADs[chr] = []
                TADs[chr].append((int(linesplit[1]),int(linesplit[2])))


        
        #if "samples_file" in data_yaml:
        #    samples_file = data_yaml["samples_file"]
        #else:
        #    sys.exit("Missing samples_file in the config file.")
        

        #if "SV_file" in data_yaml:
        #    SV_file = data_yaml["SV_file"]
        #elif "SV_dir" in data_yaml:
        #    SV_dir = data_yaml["SV_dir"]
        #else:
        #    sys.exit("Missing SV_file or SV_dir in the config file.")

        #if "CNV_dir" in data_yaml:
        #    CNV_dir = data_yaml["CNV_dir"]
        #else:
        #    print("WARNING: no CNV file was provided in the config file.")

        if "FusionTranscripts_file" in data_yaml:
            FusionTranscripts_file = data_yaml["FusionTranscripts_file"]
        #else:
        #    sys.exit("Missing FusionTranscripts_file in the config file.")

        if "SNV_RNA_dir" in data_yaml:
            SNV_RNA_dir = data_yaml["SNV_RNA_dir"]
        else:
            SNV_RNA_dir = None

        if "BAM_DNA_template" in data_yaml and "BAM_RNA_template" in data_yaml:
            BAM_DNA_template = data_yaml["BAM_DNA_template"]
            BAM_RNA_template = data_yaml["BAM_RNA_template"]
        #else:
        #    sys.exit("Missing SNV_RNA_dir in the config file.")


        # Gene expression
        if "RNA_TPM_file" in data_yaml:
            df_TPM = pd.read_csv(data_yaml["RNA_TPM_file"],sep="\t",index_col=0)
            if "gene_name" in df_TPM.columns:
                df_TPM.drop("gene_name",axis=1,inplace=True)
            # Make sure that gene IDs are unique
            gene_ids_presents=set()
            selected_indices=[]
            for i in range(df_TPM.shape[0]):
                if not df_TPM.index[i] in gene_ids_presents:
                    gene_ids_presents.add(df_TPM.index[i])
                    selected_indices.append(i)
                else:
                    print("WARNING: " + df_TPM.index[i] + " has more than 1 entry. Only the first one was kept.")
            df_TPM = df_TPM.iloc[selected_indices,:]
            print(df_TPM)
        elif "RNA_counts_file" in data_yaml and "exonic_gene_lengths" in data_yaml:
            df_counts = pd.read_csv(data_yaml["RNA_counts_file"],sep="\t",index_col=0)
            if "gene_name" in df_counts.columns:
                df_counts.drop("gene_name",axis=1,inplace=True)
            # Make sure that gene IDs are unique
            gene_ids_presents=set()
            selected_indices=[]
            for i in range(df_counts.shape[0]):
                if not df_counts.index[i] in gene_ids_presents:
                    gene_ids_presents.add(df_counts.index[i])
                    selected_indices.append(i)
                else:
                    print("WARNING: " + df_counts.index[i] + " has more than 1 entry. Only the first one was kept.")
            df_counts = df_counts.iloc[selected_indices,:]

            # Compute TPM from read counts
            exonic_gene_length = pd.read_csv(data_yaml["exonic_gene_lengths"],sep="\t",index_col=0,header=None)
            df_TPM = df_counts
            for g in df_TPM.index:
                if g in exonic_gene_length.index:
                    df_TPM.loc[g,:] = df_TPM.loc[g,:] / exonic_gene_length.loc[g,1]
                else:
                    df_TPM.drop(g,axis=0,inplace=True)

            for c in df_TPM.columns:
                df_TPM[c] = df_TPM[c]*1000000 / np.sum(df_TPM[c])

        else:
            sys.exit("Missing RNA_TPM_file or RNA_counts_file in the config file.")

        selected_geneIDs=[]
        for gene_id in df_TPM.index:
            gene = data.gene_by_id(gene_id)
            # Exclude genes with low expression for all samples, non-protein coding genes, genes on Y chromosome and in the mitochondrial genome
            if np.sum(df_TPM.loc[gene_id,:])<1: continue
            if gene.biotype in ["snRNA","pseudogene","snoRNA","antisense","lincRNA","misc_RNA","sense_intronic","processed_transcript"]: continue
            if gene.contig in ["Y","MT"]: continue
            selected_geneIDs.append(gene_id)
        df_TPM = df_TPM.loc[selected_geneIDs,:]

        if "RNA_TPM_reference_file" in data_yaml:
            df_TPM_reference = pd.read_csv(data_yaml["RNA_TPM_reference_file"],sep="\t",index_col=0)
            # Select genes present in both the reference and studied samples
            selected_geneIDs=[]
            for x in df_TPM_reference.index:
                if x in df_TPM.index: selected_geneIDs.append(x)
            df_TPM = df_TPM.loc[selected_geneIDs,:]
            df_TPM_reference = df_TPM_reference.loc[selected_geneIDs,:]

        if "RNA_TPM_other_samples" in data_yaml:
            df_TPM_other = pd.read_csv(data_yaml["RNA_TPM_other_samples"],sep="\t",index_col=0)
        if "RNA_TPM_normal_samples" in data_yaml:
            df_TPM_normal = pd.read_csv(data_yaml["RNA_TPM_normal_samples"],sep="\t",index_col=0)

        # Breakpoints
        if "breakpoints" in data_yaml:
            df_breakpoints = pd.read_csv(data_yaml["breakpoints"],sep="\t")
        else:
            sys.exit("Missing breakpoints in the config file.")

        
            

        #samples = []
        
        #samples_RNA=set(df_TPM.columns)
        #samples_breakpoints=set(df_breakpoints["sample"])
        #for sample in samples_RNA:
        #    if sample in samples_breakpoints:
        #        samples.append(sample)
        #df_TPM = df_TPM[samples]
        samples = list(df_TPM.columns)
        breakpoints = read_breakpoints(df_breakpoints=df_breakpoints)

        # CNAs
        if "CNAs" in data_yaml:
            df_CNAs = pd.read_csv(data_yaml["CNAs"],sep="\t")
            CNAs = read_CNAs(df_CNAs)

        # Enhancers
        if "enhancers_bed" in data_yaml:
            df_enhancers = pd.read_csv(data_yaml["enhancers_bed"],sep="\t",header=None)
            columns = list(df_enhancers.columns)
            columns[0:3] = ["chr","start","end"]
            df_enhancers.columns = columns
            df_enhancers["chr"] = [x.lstrip("chr") for x in df_enhancers["chr"]]

        if "output_dir" in data_yaml:
            output_dir = data_yaml["output_dir"]
            cache_dir = os.path.join(output_dir,"cache")
            os.makedirs(cache_dir, exist_ok=True)
        else:
            sys.exit("Missing output_dir in the config file.")
        




    path_corrected_TPM = os.path.join(cache_dir,"RNA_TPM_corrected.tsv")
    if os.path.exists(path_corrected_TPM):
        df_TPM = pd.read_csv(path_corrected_TPM,sep="\t",index_col=0)
    else:
        if "df_CNAs" in dir():
            print("Correcting TPM with genomic copy number")
            df_TPM = correct_TPM_CN(df_TPM)
            df_TPM.to_csv(path_corrected_TPM,sep="\t")
        else:
            df_TPM = df_TPM.copy()
            df_TPM.to_csv(path_corrected_TPM,sep="\t")


    # Find SNPs in each gene/sample
    if os.path.exists(os.path.join(cache_dir,"SNPs.pickle")):
        filename_SNP = os.path.join(cache_dir,"SNPs.pickle")
        filename_pvalues= os.path.join(cache_dir,"SNP_pvals.pickle")
        filename_llrs= os.path.join(cache_dir,"SNP_llrs.pickle")
        print("Loaded cached SNP files")
        with open(filename_SNP, "rb") as f:
            SNPs = pickle.load(f)
        with open(filename_pvalues,"rb") as f:
            pvalues_monoallelic = pickle.load(f)
        with open(filename_llrs,"rb") as f:
            llrs_monoallelic = pickle.load(f)
    else:
        print("Looking for SNPs in the RNA data")
        #if "BAM_DNA_template" in dir():
        #    SNPs = find_SNPs_DNA_RNA(samples,df_TPM.index,BAM_DNA_template,BAM_RNA_template)
        #else:
        SNPs = read_SNPs_RNA(samples,df_TPM.index,SNV_RNA_dir)
        print("Computing pvalues for monoallelic expression")
        pvalues_monoallelic = compute_monoallelic_pvalues(SNPs)
        print("Computing loglikelihood ratios between monoallelic and biallelic expression")
        llrs_monoallelic = compute_monoallelic_llrs(SNPs)
        if "cache_dir" in dir():
            filename_SNP = os.path.join(cache_dir,"SNPs.pickle")
            filename_pvalues= os.path.join(cache_dir,"SNP_pvals.pickle")
            filename_llrs= os.path.join(cache_dir,"SNP_llrs.pickle")
            with open(filename_SNP, "wb") as f:
                pickle.dump(SNPs, f)
            
            with open(filename_pvalues, "wb") as f:
                pickle.dump(pvalues_monoallelic, f)
            
            with open(filename_llrs, "wb") as f:
                pickle.dump(llrs_monoallelic, f)









