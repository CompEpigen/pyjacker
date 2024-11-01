# pyjacker

This is a tool to detect enhancer hijacking events in a cohort profiled with WGS and RNA-seq. It does not require matched normals and can detect enhancer hijacking events occurring in a single sample. Briefly, it looks for outlier high and monoallelic expression of a gene in a sample which has a breakpoint close to the gene.

## Usage
In an environment with python>=3.7: 
```
pip install pyjacker
python pyjacker.py config.yaml
```
The config file indicates all the parameters, including paths to the input files (see config_AML.yaml as an example). Alternatively, we provide a nextflow workflow that generates pyjacker's inputs from bam files, and run pyjacker: https://github.com/CompEpigen/wf_WGS.

## Inputs
### Gene expression table (required)
Rows are genes (ensembl IDs) and columns are samples. The expression data must be provided in TPM. See data/TPM_ckAML.tsv for an example.

### Breakpoints (required)
tsv file with columns: sample, chr1, pos1, chr2, pos2. 
The fields chr2 and pos2 are optional (for example if you only have copy number data). See data/breakpoints.tsv for an example.

### TADs (optional)
A bed file of topologically-associating domains can be used, in which case only the breakpoints in the same TAD as a gene are considered in the search for enhancer hijacking events. One TAD file is provided in the data directory. If not TAD file is provided, pyjacker will instead look for breakpoints within a fixed distance to the gene (1.5Mb by default). See data/HSPC_TADs.bed for an example.

### Allelic read counts at SNPs in RNAseq (optional)
This is used to detect monoallelic expression. This requires files generated by [fast_ase](https://github.com/e-sollier/fast_ase) or [GATK ASEReadCounter](https://gatk.broadinstitute.org/hc/en-us/articles/360037428291-ASEReadCounter).

### Copy number alterations (optional)
tsv file with the following columns: sample, chr, start, end, cn. See data/CNAs.tsv for an example.
If provided, this will be used to:
- correct gene expression based on copy number (so high expression because of amplification will not be reported)
- filter out SNPs within deletions from the monoallelic expression detection


### Fusion transcripts
Fusion transcripts can also lead to aberrant high and monoallelic expression of a gene. If a list of fusion transcripts detected from RNAseq is provided, they will be used to annotate candidate enhancer hijacking events which are actually due to a fusion. See data/fusions.tsv for an example file.

### Enhancers
A file of scored enhancers, generated by ROSE. See data/enhancers_myeloid.tsv for an example.
