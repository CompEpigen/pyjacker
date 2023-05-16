# pyjacker

This is a tool to detect enhancer hijacking events from RNAseq + WGS. It does not require matched normals and can detect enhancer hijacking events which occur in a single sample. Briefly, it looks for outlier high and monoallelic expression of a gene in a sample which has a breakpoint close to the gene. It does not require a reference cohort, but will instead use the samples without breakpoints close to the gene to estimate the normal gene expression.

## Usage
`python pyjacker.py config.yaml`
A template config file is provided, and has to be modified

## Dependencies
Python 3 is required, and the dependencies can be installed with the following:
```
pip install `numpy pandas scipy statsmodels pyensembl
pyensembl install --release 75 --species human
```

## Inputs
### Gene expression table (required)
Rows are genes (ensembl IDs) and columns are samples. The expression data must be provided in TPM. 

### Breakpoints (required)
tsv file with columns: sample, chr1, pos1, chr2, pos2. 
The fields chr2 and pos2 are optional (for example if you only have copy number data).

### TADs (optional)
A bed file of topologically-associating domains can be used, in which case only the breakpoints in the same TAD as a gene are considered in the search for enhancer hijacking events. One TAD file is provided in the data directory. If not TAD file is provided, pyjacker will instead look for breakpoints within a fixed distance to the gene (1.5Mb by default).

### Allelic read counts at SNPs in RNAseq (optional)
This is used to detect monoallelic expression. 

### Copy number alterations (optional)
tsv file with the folowwing columns: sample, chr, start, end, cn
If provided, this will be used to:
- correct gene expression based on copy number (so high expression because of amplification will not be reported)
- filter out SNPs within deletions from the monoallelic expression detection


### Fusion transcripts
Fusion transcripts can also lead to aberrant high and monoallelic expression of a gene. If a list of fusion transcripts detected from RNAseq is provided, they will be used to annotate candidate enhancer hijacking events which are actually due to a fusion. A bed file of enhancers detected in HSPCs is provided in the data directory.

### Enhancers
A bed file of enhancers can be provided in order to annotate the results with putative enhancers.
