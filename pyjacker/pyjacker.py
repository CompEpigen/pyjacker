import numpy as np
import pandas as pd
import random
from statsmodels.stats.multitest import multipletests
import os
import sys
import yaml
import multiprocessing
from tqdm import tqdm
import argparse

from pyjacker.genes import read_genes_gtf, index_genes_by_pos
from pyjacker.ase import compute_ase_matrix, count_SNPs_gene_sample
from pyjacker.outlier_expression import compute_ohe_score
from pyjacker.cna import read_CNAs, add_cna_breakpoints, correct_exp_cn, gene_is_deleted
from pyjacker.breakpoints import read_breakpoints, read_TADs, find_samples_with_breakpoints_near_gene, find_distance_breakpoint, find_fusion, find_enhancers_orientation
from pyjacker.html_reports import generate_main_report, generate_individual_reports


class pyjacker:
    def __init__(self,config_file):
        if not os.path.exists(config_file): sys.exit("The config file does not exist: "+str(config_file))
        with open(config_file, 'r') as stream:
            try:
                data_yaml = yaml.safe_load(stream)
            except: sys.exit("Could not load the config file: {}. Please make sure that it is in valide yaml format.".format(config_file))
            self.output_dir = data_yaml["output_dir"]
            if "n_reports" in data_yaml: self.n_reports= int(data_yaml["n_reports"])
            else: self.n_reports=-1
            if "image_format" in data_yaml: self.image_format=data_yaml["image_format"]
            else: self.image_format="svg"
            if "image_dpi" in data_yaml: self.image_dpi= float(data_yaml["image_dpi"])
            else: self.image_dpi=300
            self.cache_dir = os.path.join( self.output_dir,"cache")
            os.makedirs( self.cache_dir,exist_ok=True)

            if "n_threads" in data_yaml: self.n_threads=data_yaml["n_threads"]
            else: self.n_threads = 6

            # Params
            self.weights = {"OHE":4.0,"ASE":2.0,"enhancers":1.0,"deletion":5}
            if "weight_OHE" in data_yaml:  self.weights["OHE"] = float(data_yaml["weight_OHE"]) # outlier high expression
            if "weight_ASE" in data_yaml:  self.weights["ASE"] = float(data_yaml["weight_ASE"]) # allele-specific expression
            if "weight_enhancers" in data_yaml:  self.weights["enhancers"] = float(data_yaml["weight_enhancers"])
            if "weight_deletion" in data_yaml:  self.weights["deletion"] = float(data_yaml["weight_deletion"])

            # gtf
            if "gtf" in data_yaml:
                self.gtf_file = data_yaml["gtf"]
                self.genes = read_genes_gtf(self.gtf_file)
                self.genes={gene_id:gene for (gene_id,gene) in self.genes.items() 
                            if (gene.biotype in ["protein_coding",""] and gene.chr in [str(x) for x in range(1,23)] + ["X"])}
                self.genes_index = index_genes_by_pos(self.genes)
            else:
                sys.exit("Missing gtf file")

            # TADs
            self.TADs = None
            if "TADs_file" in data_yaml:
                TADs_file = data_yaml["TADs_file"]
                self.TADs = read_TADs(TADs_file)
            self.max_dist_bp2tss=1500000 # will only be used if no TADs are provided
            if "max_dist_bp2tss" in data_yaml: 
                self.max_dist_bp2tss= data_yaml["max_dist_bp2tss"]

            #cytobands
            if "cytobands" in data_yaml:
                self.cytobands = data_yaml["cytobands"]
                self.chr_lengths={}
                with open(self.cytobands,"r") as infile:
                    for line in infile:
                        linesplit = line.split("\t")
                        chr=linesplit[0].lstrip("chr")
                        pos=int(linesplit[2])
                        if chr in self.chr_lengths: self.chr_lengths[chr]=max(self.chr_lengths[chr],pos)
                        else: self.chr_lengths[chr]=pos
            else:
                self.cytobands = None
                self.chr_lengths=None
        
            # Gene expression
            if "RNA_TPM_file" in data_yaml:
                self.df_TPM = pd.read_csv(data_yaml["RNA_TPM_file"],sep="\t",index_col=0)
                if "gene_name" in self.df_TPM.columns:
                    self.df_TPM.drop("gene_name",axis=1,inplace=True)
                self.samples = list(self.df_TPM.columns)
                self.genes={gene_id:gene for (gene_id,gene) in self.genes.items() if gene_id in self.df_TPM.index}   
            else:
                sys.exit("Missing RNA_TPM_file from the config file.")

            self.df_TPM_normal= None
            if "RNA_TPM_normal_samples" in data_yaml:
                self.df_TPM_normal = pd.read_csv(data_yaml["RNA_TPM_normal_samples"],sep="\t",index_col=0)


            # Breakpoints
            if (not "breakpoints" in data_yaml) and (not "CNAs" in data_yaml): 
                sys.exit("Missing breakpoints (or at least CNAs) in the config file.")
            self.df_breakpoints = None
            if "breakpoints" in data_yaml:
                self.df_breakpoints = pd.read_csv(data_yaml["breakpoints"],sep="\t")
                # Remove breakpoints corresponding to samples for which we have no expression.
                selected_indices=[]
                for x in self.df_breakpoints.index:
                    if self.df_breakpoints.loc[x,"sample"] in self.samples: selected_indices.append(x)
                self.df_breakpoints = self.df_breakpoints.loc[selected_indices,:].reset_index(drop=True)

            if len(self.samples)<3: sys.exit("Error: Insufficient number of samples in the RNA TPM file. Pyjacker is designed to run with at least 10 samples, and cannot be run with fewer than 3.")
            print("Running pyjacker on {} samples.".format(len(self.samples)))
            if len(self.samples)<10:
                print("WARNING: pyjacker is designed to run with at least 10 samples.")

            # CNAs
            if "CNAs" in data_yaml:
                df_CNAs = pd.read_csv(data_yaml["CNAs"],sep="\t",dtype={"chr":str})
                self.CNAs = read_CNAs(df_CNAs,self.samples)
                self.df_breakpoints = add_cna_breakpoints(self.CNAs,self.df_breakpoints,self.chr_lengths)
                TPM_corrected_file = os.path.join(self.cache_dir,"TPM_corrected.tsv")
                if os.path.exists(TPM_corrected_file):
                    self.df_TPM = pd.read_csv(TPM_corrected_file,sep="\t",index_col=0)
                else:
                    print("Computing gene expression values corrected for copy number...")
                    self.df_TPM = correct_exp_cn(self.df_TPM,self.CNAs,self.genes,chr_lengths=self.chr_lengths)
                    self.df_TPM.to_csv(TPM_corrected_file,sep="\t")
                
            else: self.CNAs = None
            self.breakpoints = read_breakpoints(self.df_breakpoints)

            # Filter genes: the gene must be expressed in at least one sample.
            # Also require at least 30% of samples to have expression lower than max/10. 
            # This reduces the number of tests.
            self.genes={gene_id:gene for (gene_id,gene) in self.genes.items() 
                        if (np.sum(self.df_TPM.loc[gene_id,:])>=1
                            and np.percentile(self.df_TPM.loc[gene_id,:],30)<=np.max(self.df_TPM.loc[gene_id,:])*0.1)}   
            self.gene_list = [self.genes[x] for x in self.genes]

            # Allele specific expression
            if "ase_dir" in data_yaml:
                self.ase_dir = data_yaml["ase_dir"]
            else:
                self.ase_dir = None
                print("WARN: ase_dir was not provided, so allele-specific expression will not be used.")
            if "ase_dna_dir" in data_yaml:
                self.ase_dna_dir = data_yaml["ase_dna_dir"]
            else:
                self.ase_dna_dir = None
            cache_ase_file = os.path.join(self.cache_dir,"ase.tsv")
            if os.path.exists(cache_ase_file):
                self.df_ase = pd.read_csv((cache_ase_file),sep="\t",index_col=0)
            else:
                print("Computing allele-specific expression scores...")
                imprinted_genes_file = None
                if "imprinted_genes_file" in data_yaml:
                    imprinted_genes_file = data_yaml["imprinted_genes_file"]
                self.df_ase = compute_ase_matrix(self.samples,self.ase_dir,self.genes,self.genes_index,imprinted_genes_file=imprinted_genes_file,CNAs=self.CNAs)
                self.df_ase.to_csv(cache_ase_file,sep="\t")

            # Fusion transcripts
            if "fusions" in data_yaml:
                self.df_fusions = pd.read_csv(data_yaml["fusions"],sep="\t")
            else:
                self.df_fusions=None

            # Enhancers
            if "enhancers" in data_yaml:
                self.df_enhancers = pd.read_csv(data_yaml["enhancers"],sep="\t",comment="#",dtype={"CHROM":str})
                self.df_enhancers["CHROM"] = [x.lstrip("chr") for x in self.df_enhancers["CHROM"]]
            else:
                self.df_enhancers=None


            if "group_by_gene" in data_yaml:
                self.group_by_gene = data_yaml["group_by_gene"]
            else:
                self.group_by_gene = True

            # FDR
            if "n_iterations_FDR" in data_yaml:
                self.n_iterations_FDR = data_yaml["n_iterations_FDR"]
            else:
                self.n_iterations_FDR = 50
                


    def find_EH(self,random_candidates=False):
        """
        Find putative enhancer hijacking events, and score them.
        If random_candidates is true, samples with rearrangements will be randomly selected for each gene, in order to estimate a null distribution.
        This is parallelized over genes.
        """

        # Data that is copied for each process.
        data_dic={"breakpoints":self.breakpoints,"CNAs":self.CNAs,"genes":self.genes,"genes_index":self.genes_index,"chr_lengths":self.chr_lengths,
                          "df_TPM":self.df_TPM,"df_TPM_normal":self.df_TPM_normal,"df_fusions":self.df_fusions,"df_enhancers":self.df_enhancers,"samples":self.samples,
                          "df_ase":self.df_ase,"ase_dir":self.ase_dir,"TADs":self.TADs,"max_dist_bp2tss":self.max_dist_bp2tss,"weights":self.weights,
                          "random_candidates": random_candidates}

        #datas=[]
        #for gl in gene_lists:
        #    data_dic = {"gene_list":gl,"breakpoints":self.breakpoints,"CNAs":self.CNAs,"genes":self.genes,"genes_index":self.genes_index,"chr_lengths":self.chr_lengths,
        #                  "df_TPM":self.df_TPM,"df_TPM_normal":self.df_TPM_normal,"df_fusions":self.df_fusions,"df_enhancers":self.df_enhancers,"samples":self.samples,
        #                  "df_ase":self.df_ase,"ase_dir":self.ase_dir,"TADs":self.TADs,"max_dist_bp2tss":self.max_dist_bp2tss,"weights":self.weights,
        #                  "random_candidates": random_candidates}
        #
        #    datas.append(data_dic)

        # Multiprocessing
        def init_worker():
            global data
            data = data_dic
        with multiprocessing.Pool(self.n_threads, initializer=init_worker) as pool:
            results = list(tqdm(pool.imap(find_EH_gene,self.gene_list),total=len(self.gene_list),disable=random_candidates,file=sys.stdout))

        results_flattened = [x for y in results for x in y]

        df_result = pd.DataFrame(results_flattened)
        if not df_result.empty:
            df_result = df_result.sort_values("score",ascending=False)
        #print(df_result)
        return df_result
    
    
        
        

    def save_results(self,df_result,group_by_gene=True,null_distribution=None):

        if group_by_gene:
            df_result = compute_grouped_gene_scores(df_result)
            gene_ids,gene_scores = get_unique_gene_scores(df_result)
            

        if null_distribution is not None:
            def compute_empirical_pval(score,null_distribution):
                return np.sum([score<x for x in null_distribution]) / len(null_distribution)
            if group_by_gene:
                pvals = [compute_empirical_pval(score,null_distribution) for score in gene_scores]
                qvals = multipletests(pvals,alpha=0.05,method="fdr_bh")[1]
                gene2qval={}
                for i in range(len(gene_ids)):
                    gene2qval[gene_ids[i]] = qvals[i]
                qvals_all = [gene2qval[x] for x in df_result["gene_id"]]
                df_result["FDR"] = qvals_all
                df_result.sort_values(["FDR","gene_name","score"],ascending=[True,True,False],inplace=True)
            else:
                df_result["pval"] = [compute_empirical_pval(score,null_distribution) for score in df_result["score"]]
                df_result["FDR"] = multipletests(df_result["pval"],alpha=0.05,method="fdr_bh")[1]

        df_result = df_result.loc[df_result["score"]>0,:]
        df_result = df_result.reset_index(drop=True) 
        columns_output = ["FDR","score","gene_id","gene_name","chr","start","end","sample","distance_to_breakpoint","n_SNPs","fusion","OHE_score","n_std","ASE_score","enhancer_score","penalty_deletion","gene_score"]
        if "super_enhancers" in df_result.columns:
            columns_output+= ["super_enhancers"]
        if "enhancers" in df_result.columns:
            columns_output+= ["enhancers"]
        df_result = df_result[columns_output]
        os.makedirs(self.output_dir,exist_ok=True)
        df_result.to_csv(os.path.join(self.output_dir,"enhancer_hijacking.tsv"),index=False,sep="\t")
        generate_main_report(df_result,self.output_dir,300,filter_monoallelic=False)
        print("Finished searching for putative enhancer hijacking events. Found {} events with FDR<=20%.".format(np.sum(df_result["FDR"]<=0.20)))
        if self.n_reports!=0:
            print("Generating reports...")
            generate_individual_reports(df_result,self.df_TPM,self.breakpoints,self.CNAs,self.genes,self.ase_dir,self.ase_dna_dir,
                                        self.gtf_file,self.output_dir,self.cytobands,n_events=self.n_reports,
                                        image_format=self.image_format,image_dpi=self.image_dpi,n_threads=self.n_threads)
        print("Pyjacker completed successfully! The results are stored in "+self.output_dir+".")

    def estimate_null_distribution(self,seed=0):
        print("Estimating the null distribution of pyjacker scores in the absence of enhancer hijacking...")
        if self.n_iterations_FDR==0: return None
        random.seed(seed)
        scores=[]
        for i in tqdm(range(self.n_iterations_FDR),file=sys.stdout):
            df_result = self.find_EH(random_candidates=True)
            if self.group_by_gene:
                df_result = compute_grouped_gene_scores(df_result)
                scores+= get_unique_gene_scores(df_result)[1]
            else:
                scores+= list(df_result["score"])
        return sorted(scores)

def find_EH_gene(gene):

    l=[]

    # Identify candidate samples (near breakpoints) and reference samples (others)
    candidate_samples = find_samples_with_breakpoints_near_gene(data["samples"],data["genes"][gene.gene_id],data["breakpoints"],data["TADs"],data["max_dist_bp2tss"])
    reference_samples = []
    for sample in data["samples"]:
        if (not sample in candidate_samples):
            reference_samples.append(sample)

    # For estimating the null distribution: select some candidates among the reference samples, and remove them from the reference.
    # That way, all "candidate" samples should not be true enhancer hijacking events, so we can estimate the null distribution.
    if data["random_candidates"]:
        if len(data["samples"])>10:
            n_candidates = min(random.randint(1,3),len(reference_samples)-8)
        else:
            n_candidates = min(random.randint(1,3),len(reference_samples)-4)
        if n_candidates<=0: return []
        previous_candidates = candidate_samples
        candidate_samples = random.sample(reference_samples,n_candidates)
        for x in candidate_samples: reference_samples.remove(x)
        reference_samples = list(reference_samples) + list(previous_candidates) 
    
    for sample in candidate_samples:
        n_std,ohe_score = compute_ohe_score(data["df_TPM"],gene.gene_id,reference_samples,sample)
        d={}
        d["gene_id"]=gene.gene_id
        d["gene_name"]=gene.gene_name
        d["chr"]=gene.chr
        d["start"]=gene.start
        d["end"]=gene.end
        d["sample"]=sample
        d["distance_to_breakpoint"]=find_distance_breakpoint(data["breakpoints"],sample,gene)
        d["n_SNPs"]=count_SNPs_gene_sample(data["ase_dir"],sample,data["genes"][gene.gene_id])
        d["fusion"]=find_fusion(data["breakpoints"],data["genes_index"],sample,gene,df_fusions=data["df_fusions"])
        if data["df_enhancers"] is not None:
            d["enhancers"],d["super_enhancers"],enhancer_score = find_enhancers_orientation(data["breakpoints"],data["TADs"],data["df_enhancers"],sample,gene)
        else: d["enhancers"],d["super_enhancers"],enhancer_score="","",0


        # Score 
        OHE_score = data["weights"]["OHE"]*ohe_score
        d["OHE_score"]=OHE_score
        d["n_std"]=n_std
        ASE_score = data["weights"]["ASE"] * data["df_ase"].loc[gene.gene_id,sample]
        d["ASE_score"]=ASE_score
        enhancer_score=data["weights"]["enhancers"]* (enhancer_score)
        d["enhancer_score"] = enhancer_score
        penalty_deletion = - data["weights"]["deletion"] * gene_is_deleted(data["CNAs"],sample,data["genes"][gene.gene_id])
        d["penalty_deletion"]=penalty_deletion

        score = OHE_score + ASE_score + enhancer_score + penalty_deletion
        d["score"]=score
        l.append(d)
    return l

def compute_grouped_gene_scores(df_result):
        """Compute a score for each gene, eg the score is higher if the gene is overexpressed in several samples."""
        gene2index = {}
        for x in df_result.index:
            gene_id = df_result.loc[x,"gene_id"]
            if not gene_id in gene2index: gene2index[gene_id] = []
            gene2index[gene_id].append(x)
        gene2score={}
        gene2samples={}
        for gene_id in gene2index:
            scores=[]
            max_score = -100000
            for x in gene2index[gene_id]:
                if df_result.loc[x,"score"]>=1:
                    scores.append(df_result.loc[x,"score"])
                    if not gene_id in gene2samples: gene2samples[gene_id] = []
                    gene2samples[gene_id].append(df_result.loc[x,"score"])
                max_score = max(max_score,df_result.loc[x,"score"])
            if max_score <=1: gene2score[gene_id] = max_score
            else:
                gene2score[gene_id] = 5 * np.sum(scores) / (len(scores) +4) # Penalize by number of candidate samples.
               
        gene_scores=[]
        for x in df_result.index:
            gene_scores.append(gene2score[df_result.loc[x,"gene_id"]])
        df_result["gene_score"] = gene_scores
        return df_result

def get_unique_gene_scores(df_result):
    gene_scores=[]
    gene_ids=[]
    for x in df_result.index:
        if not df_result.loc[x,"gene_id"] in gene_ids:
            gene_scores.append(df_result.loc[x,"gene_score"])
            gene_ids.append(df_result.loc[x,"gene_id"])
    return gene_ids,gene_scores

def main():
    parser = argparse.ArgumentParser(prog='pyjacker', description='Tool for the detection of enhancer hijacking using WGS and RNAseq. Please check the documentation at https://github.com/CompEpigen/pyjacker for a description of the config file that must be provided as input.')
    parser.add_argument('config_file') 
    args = parser.parse_args()
    pyjack = pyjacker(args.config_file)
    print("Searching for enhancer hijacking events...")
    df_result=pyjack.find_EH(random_candidates=False)
    null_distribution = pyjack.estimate_null_distribution()
    pyjack.save_results(df_result,null_distribution=null_distribution)


if __name__=="__main__":
    main()