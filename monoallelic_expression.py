import os
import numpy as np
from scipy.stats import betabinom, combine_pvalues
from sklearn.feature_selection import GenericUnivariateSelect
import pysam
import vcfpy

from collections import namedtuple
SNP = namedtuple('SNP', 'chr pos ref alt RO AD popAF')

import config

def pval_betabinom(k,n,alpha=10,beta=10):
    if k<n/2:
        return 2*betabinom.cdf(k,n,alpha,beta)
    else:
        return 2*betabinom.cdf(n-k,n,alpha,beta)

def llr_betabinom(k,n,alpha=10,beta=10):
    loglik_biallelic = betabinom.logpmf(k,n,alpha,beta)

    lik_monoallelic1 = betabinom.pmf(k,n,1,9)
    lik_monoallelic2 = betabinom.pmf(k,n,9,1)
    loglik_monoallelic = np.log(0.5 * lik_monoallelic1 + 0.5 * lik_monoallelic2)
    return loglik_monoallelic - loglik_biallelic

def plot_betabinom_distributions():
    import matplotlib.pyplot as plt
    n=200
    x = np.linspace(0,1,40)
    y_biallelic = [betabinom.pmf(int(n*a),n,10,10) for a in x]
    y_monoallelic = [0.5*betabinom.pmf(int(n*a),n,1,9) + 0.5*betabinom.pmf(int(n*a),n,9,1) for a in x]
    plt.plot(x,y_biallelic,label="biallelic")
    plt.plot(x,y_monoallelic,label="monoallelic")
    plt.legend()
    plt.show()


'''
def read_SNPs_RNA(samples,gene_ids,SNV_RNA_dir,min_depth=6):
    d= {}
    for sample in samples:
        print(sample)
        d[sample] = {}
        filepath = os.path.join(SNV_RNA_dir,sample+".vcf.gz")
        if not os.path.exists(filepath):
            filepath = os.path.join(SNV_RNA_dir,sample,sample+".vcf.gz")
        reader = vcfpy.Reader.from_path(filepath)
        c=0
        for gene_id in gene_ids:
            c+=1
            if c%1000==0: print(c)
            gene = config.data.gene_by_id(gene_id)
            d[sample][gene_id] = []
            #try:
            for record in reader.fetch(gene.contig,begin=gene.start,end=gene.end):
                #if len(record.REF)>1 or len(record.ALT[0].value)>1:continue
                ref_counts = record.calls[0].data.get("AD")[0]
                alt_counts = record.calls[0].data.get("AD")[1]
                if ref_counts+alt_counts>=min_depth:
                    pval = pval_betabinom(alt_counts,alt_counts+ref_counts,10,10)
                    d[sample][gene_id].append(SNP(record.CHROM,record.POS,record.REF,str(record.ALT[0].value),ref_counts,alt_counts,pval))
            #except:
            #    pass
    return d
'''

def read_SNPs_RNA(samples,gene_ids,SNV_RNA_dir,min_depth=6):
    d= {}
    for sample in samples:
        d[sample] = {}
        for gene_id in gene_ids: d[sample][gene_id] = []
        if SNV_RNA_dir is not None:
            filepath = os.path.join(SNV_RNA_dir,sample+".vcf.gz")
            if not os.path.exists(filepath):
                filepath = os.path.join(SNV_RNA_dir,sample,sample+".vcf.gz")
                if not os.path.exists(filepath):
                    print("No RNA SNP file for sample " + sample)
                    continue
            reader = vcfpy.Reader.from_path(filepath)
            for record in reader:
                if len(record.REF)>1 or len(record.ALT[0].value)>1:continue
                genes = config.data.genes_at_locus(record.CHROM,record.POS)
                if len(record.calls[0].data.get("AD")) !=2: continue
                ref_counts = record.calls[0].data.get("AD")[0]
                alt_counts = record.calls[0].data.get("AD")[1]
                popAF=0
                if "popAF" in record.INFO: popAF = record.INFO["popAF"][0]
                for gene in genes:
                    gene_id = gene.gene_id
                    if gene_id in gene_ids and ref_counts+alt_counts>=min_depth:
                        d[sample][gene_id].append(SNP(record.CHROM,record.POS,record.REF,str(record.ALT[0].value),ref_counts,alt_counts,popAF))
    return d

def read_SNPs_RNA_merged(samples,gene_ids,SNV_RNA_vcf,min_depth=6):
    """Read a merged VCF file containing the SNPs detected in RNAseq data for all samples."""
    d= {}
    for sample in samples:
        d[sample] = {}
        for gene_id in gene_ids: d[sample][gene_id] = []

    reader = vcfpy.Reader.from_path(SNV_RNA_vcf)
    for record in reader:
        if len(record.REF)>1 or len(record.ALT[0].value)>1:continue
        genes = config.data.genes_at_locus(record.CHROM,record.POS)
        for call in record.calls:
            if len(call.data.get("AD")) ==2:
                sample = call.sample.split("_")[-1]
                ref_counts = call.data.get("AD")[0]
                alt_counts = call.data.get("AD")[1]
                for gene in genes:
                    gene_id = gene.gene_id
                    if gene_id in gene_ids and ref_counts+alt_counts>=min_depth:
                        d[sample][gene_id].append(SNP(record.CHROM,record.POS,record.REF,str(record.ALT[0].value),ref_counts,alt_counts))
    return d

def find_SNPs_DNA_RNA_sample_gene(sample,gene,min_depth=6):
    DNA_path = config.BAM_DNA_template.replace("{sample}",sample)
    RNA_path = config.BAM_RNA_template.replace("{sample}",sample)
    if not(os.path.exists(DNA_path) and os.path.exists(RNA_path)):
        print((DNA_path,RNA_path))
        print("Missing bam files for sample "+sample)
        return []
    samfile_DNA = pysam.AlignmentFile(DNA_path, "rb" )
    samfile_RNA = pysam.AlignmentFile(RNA_path, "rb" )
    variants_DNA=[]
    for pileupcolumn in samfile_DNA.pileup(gene.contig, gene.start,gene.end,truncate=True,max_depth=70): 
        bases={}
        for x in pileupcolumn.pileups:
            if not x.is_del and not x.is_refskip:
                b=x.alignment.query_sequence[x.query_position]
                if not b in bases: bases[b]=1
                else: bases[b]+=1
        base_nucleotides=[]
        base_counts=[]
        for x in bases:
            if bases[x]>5 and bases[x]/pileupcolumn.n>=0.30:
                base_counts.append(bases[x])
                base_nucleotides.append(x)
        if len(base_counts)==2: 
            variants_DNA.append(SNP(gene.contig,pileupcolumn.pos,base_nucleotides[0],base_nucleotides[1],base_counts[0],base_counts[1],0))
    variants_RNA=[]
    for variant in variants_DNA:
        for pileupcolumn in samfile_RNA.pileup(variant.chr, variant.pos,variant.pos+1,truncate=True,max_depth=300):
            bases={}
            total_counts=0
            for x in pileupcolumn.pileups:
                if not x.is_del and not x.is_refskip:
                    b=x.alignment.query_sequence[x.query_position]
                    if not b in bases: bases[b]=1
                    else: bases[b]+=1
                    total_counts+=1
            if total_counts>=min_depth:
                counts1 = bases[variant.ref] if variant.ref in bases else 0
                counts2 = bases[variant.alt] if variant.alt in bases else 0
                variants_RNA.append(SNP(gene.contig,pileupcolumn.pos,variant.ref,variant.alt,counts1,counts2,0))
    return variants_RNA


def find_SNPs_DNA_RNA(samples,gene_ids,min_depth=6):
    d= {}
    for sample in samples:
        d[sample] = {}
        for gene_id in gene_ids:
            d[sample][gene_id] = []
        for gene_id in gene_ids:
            gene = config.data.gene_by_id(gene_id)
            d[sample][gene_id] = find_SNPs_DNA_RNA_sample_gene(sample,gene,min_depth)
    return d



#d = read_SNPs_RNA(config.samples,config.df_counts.index,config.SNV_RNA_dir)
#d["C010-AML-15KM19129"][config.data.genes_by_name("BCL11B")[0].gene_id]

#d2 = read_SNPs_RNA_merged(config.samples,config.df_counts.index,"/home/e840r/Documents/AML/out/test_merged.vcf.gz")

#reader = vcfpy.Reader.from_path("/home/e840r/Documents/AML/out/test_merged.vcf.gz")
#record = next(reader)


def compute_monoallelic_pvalues(SNPs):
    pvalues={}
    for sample in SNPs:
        pvalues[sample]={}
        for gene in SNPs[sample]:
            pvals=[]
            for snp in SNPs[sample][gene]:
                pvals.append(pval_betabinom(snp.AD,snp.RO+snp.AD,10,10))
            if len(pvals)==0:
                pval=1
            else:
                pval = combine_pvalues(pvals)[1]
            pvalues[sample][gene] = pval
    return pvalues


def compute_monoallelic_llrs_sample_gene(SNPs):
    llrs=[]
    for snp in SNPs:
        llrs.append(llr_betabinom(snp.AD,snp.RO+snp.AD,10,10))
    return np.sum(llrs)

def compute_monoallelic_llrs(SNPs):
    loglikratios={}
    for sample in SNPs:
        loglikratios[sample]={}
        for gene in SNPs[sample]:
            loglikratios[sample][gene] = compute_monoallelic_llrs_sample_gene(SNPs[sample][gene])
    return loglikratios





def find_samples_monoallelic_expression(gene_id,threshold=4):
    if not "llrs_monoallelic" in dir(config):
        return []
    else:
        samples_monoallelic = []
        for sample in config.samples:
            if config.llrs_monoallelic[sample][gene_id]>=threshold:
                samples_monoallelic.append(sample)
        return samples_monoallelic

def monoallelic_expression2(sample,gene):
    if not "SNV_RNA_dir" in dir(config):
        return "?",0,1
    filepath = os.path.join(config.SNV_RNA_dir,sample+".vcf.gz")
    if not os.path.exists(filepath):
        filepath = os.path.join(config.SNV_RNA_dir,sample,sample+".vcf.gz")
    if os.path.exists(filepath) and os.path.exists(filepath+".tbi"):
        count_het=0
        count_hom=0
        reader = vcfpy.Reader.from_path(filepath)
        pval_biallelics=[]
        try:
            for record in reader.fetch(gene.contig,begin=gene.start,end=gene.end):
                #if len(record.REF)>1 or len(record.ALT[0].value)>1:continue
                #print(record.calls[0].data.get("AD"))
                ref_counts = record.calls[0].data.get("AD")[0]
                alt_counts = record.calls[0].data.get("AD")[1]
                if ref_counts <=1 or alt_counts <=1:
                    count_hom+=1
                else:
                    count_het+=1
                pval_biallelics.append(pval_betabinom(alt_counts,alt_counts+ref_counts,10,10))
            if len(pval_biallelics)==0:
                pval_biallelic=1
                n_SNPs=0
            else:
                pval_biallelic = combine_pvalues(pval_biallelics)[1]
                n_SNPs=len(pval_biallelics)
            if count_hom+count_het==0:
                return "?",n_SNPs,pval_biallelic
            elif count_hom>=2*count_het:
                return "monoallelic",n_SNPs,pval_biallelic
            else:
                return "biallelic",n_SNPs,pval_biallelic
        except:
            return "?",0,1
    else:
        return "?",0,1