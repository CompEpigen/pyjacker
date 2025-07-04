# pyjacker
[![DOI](https://zenodo.org/badge/641412900.svg)](https://doi.org/10.5281/zenodo.14227955)

This is a tool to detect enhancer hijacking events in a cohort of at least 10 samples profiled with WGS and RNA-seq. It does not require matched normals and can detect enhancer hijacking events occurring in only a single sample. Briefly, it looks for outlier high and monoallelic expression of a gene in a sample which has a breakpoint close to the gene.

## Usage
In an environment with python>=3.7: 
```
pip install pyjacker
pyjacker config.yaml
```
The config file indicates all the parameters, including paths to the input files (see config_AML.yaml as an example). Alternatively, we provide a nextflow workflow that generates pyjacker's inputs from bam files, and run pyjacker: https://github.com/CompEpigen/wf_WGS.

### Docker
Alternatively, pyjacker can be run using a docker image.
```
docker run -t -w `pwd` -v `pwd`:`pwd` esollier/pyjacker:latest pyjacker config.yaml
```

## Inputs
### Gene expression table (required)
Rows are genes (ensembl IDs) and columns are samples. The expression data must be provided in TPM. See data/TPM_ckAML.tsv for an example.

### Breakpoints (required)
tsv file with columns: sample, chr1, pos1, chr2, pos2. 
The fields chr2 and pos2 are optional (for example if you only have copy number data). See data/breakpoints.tsv for an example.

### Reference files: gtf and cytobands (required)
gtf file containing gene coordinates for your reference genome (see [data/Homo_sapiens.GRCh37.75.gtf.gz](https://github.com/CompEpigen/pyjacker/blob/main/data/Homo_sapiens.GRCh37.75.gtf.gz) or [data/Homo_sapiens.GRCh38.113.gtf.gz](https://github.com/CompEpigen/pyjacker/blob/main/data/Homo_sapiens.GRCh38.113.gtf.gz), and tsv file containing cytobands (see [data/cytobands_hg19](https://github.com/CompEpigen/pyjacker/blob/main/data/cytobands_hg19.tsv) or [data/cytobands_hg38](https://github.com/CompEpigen/pyjacker/blob/main/data/cytobands_hg38.tsv)).

### TADs (optional)
A bed file of topologically-associating domains can be used, in which case only the breakpoints in the same TAD as a gene are considered in the search for enhancer hijacking events. See  [data/TADs_Dixon_IMR90_hg19.bed](https://github.com/CompEpigen/pyjacker/blob/main/data/TADs_Dixon_IMR90_hg19.bed) or [data/TADs_Dixon_IMR90_hg38.bed](https://github.com/CompEpigen/pyjacker/blob/main/data/TADs_Dixon_IMR90_hg38.bed) for TADs derived from the data of [Dixon et al.](https://www.nature.com/articles/nature14222) or [data/TADs_HSPC_hg19.bed](https://github.com/CompEpigen/pyjacker/blob/main/data/TADs_HSPC_hg19.bed) for TADs derived from [HSPCs](https://doi.org/10.1182/bloodadvances.2023012161). If not TAD file is provided, pyjacker will instead look for breakpoints within a fixed distance to the gene (1.5Mb by default).

### Allelic read counts at SNPs in RNAseq (optional)
This is used to detect monoallelic expression. This requires files generated by [fast_ase](https://github.com/e-sollier/fast_ase) or [GATK ASEReadCounter](https://gatk.broadinstitute.org/hc/en-us/articles/360037428291-ASEReadCounter). See [data/ASE_ckAML](https://github.com/CompEpigen/pyjacker/tree/main/data/ASE_ckAML) for example files.

### Copy number alterations (optional)
tsv file with the following columns: sample, chr, start, end, cn. See [data/CNAs_ckAML.tsv](https://github.com/CompEpigen/pyjacker/blob/main/data/CNAs_ckAML.tsv) for an example.
If provided, this will be used to:
- correct gene expression based on copy number (so high expression because of amplification will not be reported)
- filter out SNPs within deletions from the monoallelic expression detection

### Enhancers (optional)
A file of scored enhancers, generated by ROSE. See [data/enhancers_myeloid_hg19.tsv](https://github.com/CompEpigen/pyjacker/blob/main/data/enhancers_myeloid_hg19.tsv) for an example.

### Fusion transcripts (optional)

Fusion transcripts can also lead to aberrant high and monoallelic expression of a gene. If a list of fusion transcripts detected from RNAseq is provided, they will be used to annotate candidate enhancer hijacking events which are actually due to a fusion. See [data/fusions_ckAML.tsv](https://github.com/CompEpigen/pyjacker/blob/main/data/fusions_ckAML.tsv) for an example file.

## Outputs

Pyjacker will write its outputs to the `output_dir` specified in the config file. The output consists of two main files:
* `report.html`: an html file containing the list of putative enhancer hijacking events, with links to detailed reports for each event.
* `enhancer_hijacking.tsv`: the same list of putative enhancer hijacking events, but in tsv format (easier to parse but less convenient to explore)

### report.html
The main html file contains a table with the ranked list of pairs (gene, sample) where the gene is predicted to be overexpressed in the sample due to a structural rearrangement. The fields are:
* Rank: the rank of the event (first events are the ones for which pyjacker is most confident that a structural rearrangement led to the gene overexpression)
* FDR: False discovery rate. If you have 100 events with an FDR of 20%, then approximately 80 events are true (the overexpression is caused by a structural rearrangement) and 20 are false (the overexpression is due to chance). Note that the FDR is computed per gene, not per pair (gene, sample): results are aggregated across samples for each gene, so that recurrent events present in multiple samples have a lower FDR.
* score: pyjacker score for the (gene, sample) pair. A high score means that the gene is likely overexpressed due to enhancer hijacking, but is somewhat arbitrary and difficult to interpret (see the [manuscript](https://doi.org/10.1158/2643-3230.BCD-24-0278) for details about the methods), so it's recommended to mainly look at the FDR for interpreting results. If the same gene is activated in multiple samples, the score can help determine in which samples the evidence for enhancer hijacking is stronger.
* gene: the gene that is identified to be overexpressed in this sample
* chr, start, end: the genomic coordinates of the gene
* sample: the sample in which the putative enhancer hijacking event was observed.
* distance to breakpoint: distance between the gene and the closest breakpoint (0 if the breakpoint is within the gene)
* fusion: if a fusion transcript involving the gene was detected, it is listed here. Fusion transcripts can also lead to monoallelic overexpression, so this helps identify which events are true enhancer hijacking events (without fusions) and which are fusions. Note that this requires a list of fusion transcripts to be provided as input (see above; this needs to be generated by a tool like arriba or STAR-Fusion). If no list of fusions is provided as input, pyjacker will try to find breakpoints where one end lies within the gene and the other end lies within another gene, but this is not as reliable as identifying the fusion from RNA-seq data.
* n_SNPs: number of SNPs that were identified in the gene and used to estimate whether the expression is monoallelic.
* OHE score: outlier high expression score. A score >0 indicates that the gene is overexpressed in the sample, compared to samples which don't have breakpoints near the gene, while a score <0 indicates that the gene is not overexpressed compared to other samples.
* ASE score: allele-specific expression score. A score >0 indicates that there is evidence for monoallelic expression, while a score <0 indicates that there is evidence for bi-allelic expression. A score of 0 indicates that there is no evidence (usually because there is no SNP covered).
* Enhancer score: if a list of scored enhancers was provided as inputs, this score indicates the strength of enhancers coming close to the gene

For each row of this main table, clicking on the sample name will lead you to a detailed report for this particular event, containing:
* an overexpression plot: all samples are ranked by the expression of the gene, and the sample of interest is colored in green (typically the expression should be much higher for this sample than for other samples).
* copy-number and sv plot: a plot indicating the copy number at each position of the chromosome where the gene lies, with arcs indicating structural variants. The position of the gene within the chromosome is highlighted.
* allele-specific expression plot: for each SNP present within the gene, shows the fraction of reads supporting the major allele (blue) and minor allele (red). If the expression is monoallelic, we should almost only see the major allele. At the bottom, lines show the position of the SNPs within the gene.
* Lists of structural variants and copy numbers present in the chromosome, and allelic read counts in the gene.

Note that these detailed reports require the "Data" directory to be in the same directory as `report.html`, so if you want to see them, you cannot just copy and paste the `report.html` file elsewhere, you also need to copy this directory.

## Runtime
Pyjacker takes approximately 5h to run on the ckAML dataset (39 samples) with default settings and 6 cores. The runtime is essentially proportional to the number of samples in the dataset and to the number of iterations used when estimating the null distribution of scores (used to compute the false discovery rate). This number of iterations is 50 by default, which ensures that accurate p-values are computed, but this can easily be reduced to 5-10 to reduce the runtime, without drastically altering the results.


## Citation
If you use pyjacker in your research, please consider citing:

Sollier E, Riedel A, Toprak UH, Wierzbinska JA, Weichenhan D, Schmid JP, Hakobyan M, Touzart A, Jahn E, Vick B, Brown-Burke F. Enhancer Hijacking Discovery in Acute Myeloid Leukemia by Pyjacker Identifies MNX1 Activation via Deletion 7q, Blood Cancer Discovery (2025). https://doi.org/10.1158/2643-3230.BCD-24-0278


