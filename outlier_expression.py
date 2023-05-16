import numpy as np
from scipy.stats import nbinom, t
from scipy.special import loggamma
from scipy.optimize import minimize
#import mpmath
#mpmath.mp.dps = 60 # decimal digits of precision

import config

from breakpoints import *
from monoallelic_expression import find_samples_monoallelic_expression

#from scipy.stats import median_absolute_deviation

def estimate_mean_disp(X):
    """Maximum likelihood estimation of mean and dispersion parameters for a negative binomial distribution"""
    #TODO : add a prior ? ie penalty for very low mean ??
    X = np.array(X)
    def nbinom_loglik(X,mu,theta):
        logliks = loggamma(X+theta) - loggamma(X+1)- loggamma(theta) + theta * np.log(theta / (mu+theta))  + X * np.log(mu/(mu+theta))

        return -np.sum(logliks)

    res = minimize(lambda x: nbinom_loglik(X,x[0],x[1]),[np.mean(X),1],bounds=[(10,None),(0.01,1000)])
    return res.x
    #mean = np.median(X)
    #var = median_absolute_deviation(X) **2
    #theta = mean**2 / (var-mean)
    #if theta > 1000: theta = 1000
    #if theta < 0.01: theta = 0.01
    #return (mean,theta)
    


def nbinom_pval(k,mu,theta):
    p = mu / (mu + mu**2 / theta)
    n = theta 
    return nbinom.sf(k-1,n,p)
    #return 2*min(1-nbinom.cdf(k-1,n,p),nbinom.cdf(k-1,n,p)) # Two-sided

#def nbinom_pval_precise(k,mu,theta):
#    mu = mpmath.mpf(mu)
#    theta = mpmath.mpf(theta)
#
#    lik = lambda x : mpmath.gamma(x+theta) / mpmath.gamma(x+1) / mpmath.gamma(theta)  *  (theta /(mu+theta))**(theta) * (mu/(mu+theta))**x
#    pval = mpmath.mpf(1.0) - mpmath.nsum(lik,[0,int(k)-1])
#    return pval


def nbinom_pvals_gene(gene_id):
    """Compute p-values for outlier high expression based on the negative binomial distribution, for the gene of interest and all samples which have SVs close"""

    #TODO: Filter out genes with deletion ??

    counts_corrected = config.df_counts_corrected.loc[gene_id,:].copy()
    
    # Find the samples which have SVs close to the gene of interest, and use all other samples to estimate the normal counts distribution for this gene
    #samples_SV = find_SVs_close(gene,breakends)
    samples_SV = find_samples_with_SV_in_same_TAD(config.data.gene_by_id(gene_id),config.breakpoints,exclude_deletion=False)
    samples_monoallelic = find_samples_monoallelic_expression(gene_id)
    if len(samples_monoallelic) > max(5,2*len(samples_SV)): 
        samples_monoallelic = [] # in case too many samples are considered monoallelic.
        #print("Warning monoallelic...")
        #print(config.data.gene_by_id(gene_id).gene_name)
    counts_normal = []
    for sample in config.samples:
        if (not sample in samples_SV) and (not sample in samples_monoallelic):
            counts_normal.append(counts_corrected[sample])
    mu,theta = estimate_mean_disp(counts_normal)
    #print((mu,theta))
    #print(samples_monoallelic)
    #print(counts_normal)
    #print(counts_corrected)

    samples_candidate=[]
    samples_SV_nodeletion= find_samples_with_SV_in_same_TAD(config.data.gene_by_id(gene_id),config.breakpoints,exclude_deletion=True)
    for sample in config.samples:
        if sample in samples_SV_nodeletion:# or sample in samples_monoallelic:
            samples_candidate.append(sample)
    # Compute p-values for the samples which have SVs close to the gene of interest
    pvals={}
    for sample in samples_candidate:
        pval = nbinom_pval(counts_corrected[sample],mu,theta)
        #if pval<=0.00001:
        #    pval = nbinom_pval_precise(counts_corrected[sample],mu,theta)
        pvals[sample] = pval
    return pvals

def OHE_pvals_gene(gene_id):
    """Compute p-values for outlier high expression based on the student distribution distribution, for the gene of interest and all samples which have SVs close"""

    TPM = config.df_TPM.loc[gene_id,:].copy()
    # Log transform 
    TPM = np.log10(0.5+TPM)
    

    # Use as normal samples those that have no breakpoint close to the gene, and which do not have monoallelic expression
    samples_SV = find_samples_with_SV_in_same_TAD(config.data.gene_by_id(gene_id),config.breakpoints,exclude_deletion=False)
    samples_monoallelic = find_samples_monoallelic_expression(gene_id)
    if len(samples_monoallelic) > max(6,3*len(samples_SV)): 
        samples_monoallelic = [] # in case too many samples are considered monoallelic.
    TPM_normal = []
    for sample in config.samples:
        if (not sample in samples_SV) and (not sample in samples_monoallelic):
            TPM_normal.append(TPM[sample])
    if "df_TPM_reference" in dir(config):
        TPM_normal = TPM_normal + list(config.df_TPM_reference.loc[gene_id,:])
    mean_normals = np.mean(TPM_normal)
    std_normals = np.std(TPM_normal)+0.2

    # Select the candidate samples with a breakpoint close to the gene, excluding those where the whole gene is deleted.
    samples_candidate=[]
    samples_SV_nodeletion= find_samples_with_SV_in_same_TAD(config.data.gene_by_id(gene_id),config.breakpoints,exclude_deletion=True,exclude_amplification=True)
    for sample in config.samples:
        if sample in samples_SV_nodeletion:# or sample in samples_monoallelic:
            samples_candidate.append(sample)
    # Compute p-values for the samples which have SVs close to the gene of interest
    pvals={}
    for sample in samples_candidate:
        # Compute pvalue based on Student's t-test distribution, estimated using only the normal samples.
        pval = t.sf((TPM[sample]-mean_normals)/std_normals / np.sqrt(1+1/(len(TPM_normal)-1)),(len(TPM_normal)-1))
        #if pval<=0.00001:
        #    pval = nbinom_pval_precise(counts_corrected[sample],mu,theta)
        pvals[sample] = pval
    return pvals

def nbinom_pvals_gene_sample(gene_id,sampleSelected):
    """Compute p-values for outlier high expression based on the negative binomial distribution, for the gene of interest and all samples which have SVs close"""

    #TODO: Filter out genes with deletion ??
    counts_corrected = config.df_counts_corrected.loc[gene_id,:].copy()

    counts_normal = []
    for sample in config.samples:
        if sample !=sampleSelected:
            counts_normal.append(counts_corrected[sample])
    mu,theta = estimate_mean_disp(counts_normal)

    pval = nbinom_pval(counts_corrected[sampleSelected],mu,theta)
    return pval
    
def nbinom_meandisp_gene(gene_id):
    """Estimate mean and overdispersion parameter for a gene"""

    counts_corrected = config.df_counts_corrected.loc[gene_id,:].copy()
    
    # Find the samples which have SVs close to the gene of interest, and use all other samples to estimate the normal counts distribution for this gene
    #samples_SV = find_SVs_close(gene,breakends)
    samples_SV = find_samples_with_SV_in_same_TAD(config.data.gene_by_id(gene_id),config.breakpoints)
    samples_monoallelic = find_samples_monoallelic_expression(gene_id)
    if len(samples_monoallelic) > max(5,3*len(samples_SV)): samples_monoallelic = [] # in case too many samples are considered monoallelic.
    counts_normal = []
    for sample in config.samples:
        if (not sample in samples_SV) and (not sample in samples_monoallelic):
            counts_normal.append(counts_corrected[sample])
    mu,theta = estimate_mean_disp(counts_normal)

    return np.mean(counts_normal),theta