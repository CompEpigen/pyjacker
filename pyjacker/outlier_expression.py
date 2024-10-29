import numpy as np



def scoring_function(x):
    if x>2: return np.log(x-1)
    else: return -2*np.log(3-x)

def compute_ohe_score(df_TPM,gene_id,reference_samples,sample):
    exp = np.log(0.5+df_TPM.loc[gene_id,:]) #1.0
    mean_reference = np.mean(exp[reference_samples])
    std_reference = np.std(exp[reference_samples]) + 0.3 # 0.5
    v = (exp[sample]-mean_reference) / std_reference
    return v, scoring_function(v)


def find_outlier_samples(df_TPM,genes=None,TissGDB_file=None):
    """
    Detects samples for which many genes have an outlier expression (regardless of breakpoints).
    These samples are likely contaminated by other cell types.
    If TissGDB is provided, it will try to guess the most likely tissue from which the contamination originated.
    Returns a dictionary sample2tissue, whose keys are the outlier samples, and the values are the likely tissues (or just 'outlier')
    """

    d={}
    for sample in df_TPM.columns: d[sample] = set()
    df = np.log(0.5+df_TPM)
    for gene_id in df.index:
        exp = sorted(df.loc[gene_id,:])[:-1]
        mean = np.mean(exp)
        std=np.std(exp) +0.3
        for sample in df.columns:
            v = (df.loc[gene_id,sample]-mean) / std
            if v>5: d[sample].add(gene_id)
    
    sample2tissue={}
    for sample in d:
        if len(d[sample])>40:
            if (genes is None) or (TissGDB_file is None):
                sample2tissue[sample] = "outlier"
            else:
                genename2tissues = {}
                with open(TissGDB_file,"r") as infile:
                    for line in infile:
                        linesplit = line.rstrip("\n").split("\t")
                        if not linesplit[0] in genename2tissues: genename2tissues[linesplit[0]] = []
                        tissue = linesplit[2]
                        if tissue.startswith("\xa0"): tissue = tissue[1:]
                        if not tissue in genename2tissues[linesplit[0]]:
                            genename2tissues[linesplit[0]].append(tissue)
                count_tissue={}
                for gene_id in d[sample]:
                    gene_name = genes[gene_id].gene_name
                    if gene_name in genename2tissues:
                        for tissue in genename2tissues[gene_name]:
                            if not tissue in count_tissue: count_tissue[tissue]=0
                            count_tissue[tissue]+=1
                max_count=0
                for tissue in count_tissue:
                    if count_tissue[tissue]>max_count: max_count = count_tissue[tissue]
                tissues = []
                if max_count<8: sample2tissue[sample] = "Possible contamination"
                else:
                    for tissue in count_tissue:
                        if count_tissue[tissue]>0.6*max_count:
                            tissues.append(tissue)
                    sample2tissue[sample] = "Possible contamination from " +",".join(tissues)
    return sample2tissue

def penalty_gene_expressed_normal(df_TPM,df_TPM_normal,gene_id,sample):
    if (not df_TPM_normal is None) and gene_id in df_TPM_normal.index:
        exp_sample = df_TPM.loc[gene_id,sample]
        max_exp_normal= np.max(df_TPM_normal.loc[gene_id,:])
        penalty = 0
        if max_exp_normal > exp_sample/10:
            return np.log((exp_sample+0.00001)/10/max_exp_normal)
        return penalty
    else:
        return 0

