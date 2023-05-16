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

    #df_result = df_result.sort_values("pval",ascending=True)
    df_result = df_result.sort_values("score",ascending=False)
    print(df_result)

    df_result["padj"] = multipletests(df_result["pval"], alpha=0.05, method='fdr_bh',is_sorted=False)[1]
    df_result = df_result.loc[df_result["pval"]<0.1,:]
    df_result = df_result.reset_index(drop=True)
    
    df_result.to_csv(os.path.join(config.output_dir,"enhancer_hijacking.tsv"),index=False,sep="\t")
    return df_result


#df = find_genes_outlier_sample("C010-AML-15PB8708")

df_result = find_enhancer_hijacking()

#df_result = pd.read_csv("/home/e840r/Documents/AML/out_new/enhancer_hijacking.tsv",sep="\t")

#generate_main_report(df_result,"/home/e840r/Documents/AML/out_new/enhancer_hijacking2")
#generate_individual_reports(df_result,"/home/e840r/Documents/AML/out_new/enhancer_hijacking2")

generate_main_report(df_result,config.output_dir,100,filter_monoallelic=False)
generate_individual_reports(df_result,config.output_dir,n_events=200)


###########
"""
import pysam
samfile = pysam.AlignmentFile("/mnt/cluster/Projects/AML/WGS/BAM/15PB8708.bam", "rb" )
samfile_RNA = pysam.AlignmentFile("/mnt/cluster/RNAseq/pipeline/star_salmon/C010-AML-15PB8708.markdup.sorted.bam", "rb" )
variants_DNA=[]
for pileupcolumn in samfile.pileup("7", 156798548,156803357,truncate=True): #
    #print(pileupcolumn)
    #print(pileupcolumn.pileups)
    #print((pileupcolumn.pos,pileupcolumn.n))
    bases={}
    for x in pileupcolumn.pileups:
        if not x.is_del and not x.is_refskip:
            #print(x.query_position)
            b=x.alignment.query_sequence[x.query_position]
            #print("qq")
            if not b in bases: bases[b]=1
            else: bases[b]+=1
    #print(bases)
    count_alleles=0
    base_counts=[]
    for x in bases:
        if bases[x]>5 and bases[x]/pileupcolumn.n>=0.30:
            count_alleles+=1
            base_counts.append(bases[x])
    if count_alleles>1 : 
        print((pileupcolumn.pos,pileupcolumn.n))
        print(bases)
        variants_DNA.append(("7",pileupcolumn.pos,base_counts[0],base_counts[1]))


for variant in variants_DNA:
    for pileupcolumn in samfile_RNA.pileup(variant[0], variant[1],variant[1]+1,truncate=True): #
        bases={}
        for x in pileupcolumn.pileups:
            if not x.is_del and not x.is_refskip:
                #print(x.query_position)
                b=x.alignment.query_sequence[x.query_position]
                #print("qq")
                if not b in bases: bases[b]=1
                else: bases[b]+=1
        print(variant)
        print(bases)
        print("---")
    #print(bases)
    count_alleles=0
    base_counts=[]
    for x in bases:
        if bases[x]>5 and bases[x]/pileupcolumn.n>=0.30:
            count_alleles+=1
            base_counts.append(bases[x])

    #print(dir(pileupcolumn))

# Monoallelic expression for one gene

gene_id = config.data.genes_by_name("HOXB3")[0].gene_id
for sample in config.samples:
    llr = config.llrs_monoallelic[sample][gene_id]
    if llr>0:
        print((sample,llr))
        if llr>5:
            print(config.SNPs[sample][gene_id])


d={"gene":[],"count":[],"del5q_association":[]}
for gene_id in config.df_TPM.index:
    count=0
    counts=np.zeros((2,2))
    for sample in config.samples:
        llr = config.llrs_monoallelic[sample][gene_id]
        del5q = 1 if (get_CN_gene(sample,config.data.genes_by_name("KDM3B")[0])==1) else 0
        if llr>5 and get_CN_gene(sample,config.data.gene_by_id(gene_id))==2:
            count+=1
            counts[1,del5q]+=1
        elif llr<-5 and get_CN_gene(sample,config.data.gene_by_id(gene_id))==2:
            counts[0,del5q]+=1
    if count>2:
        d["gene"].append(config.data.gene_by_id(gene_id).gene_name)
        d["count"].append(count)
        d["del5q_association"].append(fisher_exact(counts))
df_monoallelic = pd.DataFrame(d).sort_values(by="count",ascending=False).reset_index()



nbinom_pvals_gene(config.data.genes_by_name("KLHL20")[0].gene_id)

nbinom_meandisp_gene(config.data.genes_by_name("MECOM")[0].gene_id)

for x in config.breakpoints:
    for b in config.breakpoints[x]:
        if b.sample=="16PB3075":
            print(b)


gene = config.data.genes_by_name("AGR2")[0]
config.SNPs["C010-AML-15PB8708"][gene.gene_id]

# Plot normal distribution
import matplotlib.pyplot as plt
import numpy as np
l= [3516.7406922293894, 2005.0675688560898, 2716.81426340231, 4986.60318970359, 4163.738270168329, 3672.06676469139, 1777.3972175672102, 1959.3029122594698, 2198.07402394354, 2177.4163995011, 2803.5656840773, 2386.01974947057, 2340.5100711526197, 2009.4370325485102, 2301.58610209839, 5611.57036987276, 2296.22674100357, 3580.5234812416497, 3114.38013784841, 3302.50461181745, 3772.4291667597204, 4329.71936832808, 4211.314327635349, 3038.8651803670195, 2203.92797074474, 2097.3560837315304, 2119.34618248918, 2684.9656443794897, 4190.13858181165, 6629.418729400661, 3393.6942887975106, 3185.1772301723104, 5159.56158090179, 3105.2365745009997, 4187.34708879195, 2047.68846685212, 3161.8273024735304, 3887.58799117201, 5226.64947474222, 2324.76110051368, 6253.27978742171, 4262.26866157043, 3434.56630153743, 3375.83614638904, 5885.443268546801, 3994.21811709685, 3341.3457954812297, 1547.53700431568, 3088.7144071760604, 3239.4881718290803, 5202.936117886729]
mu = 3402
theta = 8.244
p = mu / (mu + mu**2 / theta)
n = theta 
from scipy.stats import nbinom
x = np.linspace(0,np.max(l)*2,100)
y=[nbinom.pmf(int(a),n,p) for a in x]

plt.plot(x,y)
plt.hist(l,density=True)
plt.show()


# Plot distribution of pvalues
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("TkAgg")
pvalues = list(df_result['pval'])

plt.hist(pvalues,density=True)
plt.xlabel("pvalue")
plt.show()




import matplotlib.pyplot as plt
from outlier_expression import nbinom_meandisp_gene
from sklearn.linear_model import LinearRegression
means=[]
thetas=[]
c=0
for gene_id in config.df_counts_corrected.index:
    c+=1
    if c%100==0:
        print(c)
    if c==5000: break
    mean,theta = nbinom_meandisp_gene(gene_id)
    means.append(mean)
    thetas.append(theta)

log_means = np.log(means)
log_thetas = np.log(thetas)

reg = LinearRegression().fit(log_means.reshape(-1,1), log_thetas)
x = np.linspace(np.min(log_means),np.max(log_means),100)
y = np.exp(reg.predict(x.reshape(-1,1)))


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.scatter(means,thetas)
ax.plot(np.exp(x),y,color="red")
ax.set_yscale('log')
ax.set_xscale('log')
plt.show()




# ERG
from breakpoints import find_samples_with_SV_in_same_TAD
s  = find_samples_with_SV_in_same_TAD(config.data.genes_by_name("ERG")[0],config.breakpoints)
for sample in config.samples:
    print(sample)
    print(sample in s)
    print(config.df_TPM.loc["ENSG00000157554",sample])
    if (config.df_TPM.loc["ENSG00000157554",sample] >50) and (not sample in s): print("111111111111111111111111111111111111111111111")
    print("--")


from outlier_expression import nbinom_meandisp_gene
mean,theta = nbinom_meandisp_gene("ENSG00000157554")


def find_genes_outlier_sample(sample):
    d = {"gene_id":[],"gene_name":[],"chr":[],"start":[],"end":[],"sample":[],"pval":[],"n_SNPs":[],"pval_monoallelic":[],"llr_monoallelic":[]}
    for gene_id in config.df_counts_corrected.index:
        gene = config.data.gene_by_id(gene_id)
        gene_name = gene.gene_name

        # Exclude genes with low expression for all samples, non-protein coding genes, genes on Y chromosome and in the mitochondrial genome
        if np.sum(config.df_counts_corrected.loc[gene_id,:])<1000: continue
        if gene.biotype in ["snRNA","pseudogene","snoRNA","antisense","lincRNA"]: continue
        if gene.contig in ["Y","MT"]: continue
        
        pval = nbinom_pvals_gene_sample(gene_id,sample)
        
        d["gene_id"].append(gene_id)
        d["gene_name"].append(gene_name)
        d["chr"].append(gene.contig)
        d["start"].append(gene.start)
        d["end"].append(gene.end)
        d["sample"].append(sample)
        d["pval"].append(pval)
        #monoallelic,n_SNPs,pval_monoallelic=monoallelic_expression(sample,gene)
        #d["monoallelic"].append(monoallelic)
        d["n_SNPs"].append(len(config.SNPs[sample][gene_id]))
        d["pval_monoallelic"].append(config.pvalues_monoallelic[sample][gene_id])
        d["llr_monoallelic"].append(config.llrs_monoallelic[sample][gene_id])

    df_result = pd.DataFrame(d)
    

    df_result = df_result.sort_values("pval",ascending=True)

    df_result["pval_BH"] = multipletests(df_result["pval"], alpha=0.05, method='fdr_bh',is_sorted=True)[1]
    df_result = df_result.reset_index(drop=True)
    df_result.to_csv(os.path.join(config.output_dir,"outliers_"+sample+".tsv"),index=False,sep="\t")
    return df_result

"""