import gzip
from collections import namedtuple
Gene = namedtuple('Gene', 'gene_id gene_name chr start end strand biotype')


def read_genes_gtf(gtf_path):
    d={}
    if gtf_path.endswith(".gz") : infile = gzip.open(gtf_path,"rt")
    else: infile =open(gtf_path,"r")
    for line in infile:
        if line.startswith("#"): continue
        linesplit = line.split("\t")
        if linesplit[2]!="gene": continue
        chr = linesplit[0].lstrip("chr")
        start = int(linesplit[3])
        end = int(linesplit[4])
        strand = linesplit[6]
        gene_biotype=""
        for x in linesplit[8].split(";"):
            if "gene_id" in x:
                gene_id = x[x.find("\"")+1:-1]
            elif "gene_name" in x:
                gene_name = x[x.find("\"")+1:-1]
            elif "gene_biotype" in x:
                gene_biotype = x[x.find("\"")+1:-1]
        d[gene_id] = Gene(gene_id,gene_name,chr,start,end,strand,gene_biotype)
    infile.close()
    return d


def index_genes_by_pos(genes):
    genes_index={}
    for gene_id in genes:
        gene = genes[gene_id]
        if not gene.chr in genes_index: genes_index[gene.chr] ={}
        for rounded_pos in range(gene.start//1000000,gene.end//1000000 +1):
            if not rounded_pos in genes_index[gene.chr]: genes_index[gene.chr][rounded_pos] = []
            genes_index[gene.chr][rounded_pos].append(gene)
    return genes_index


def genes_at_locus(genes_index,chr,pos):
    l=[]
    if chr in genes_index and pos//1000000 in genes_index[chr]:
        for gene in genes_index[chr][pos//1000000]:
            if pos>=gene.start and pos<=gene.end: l.append(gene)
    return l

def find_gene_id(gene_name,gtf_file):
    gene_id=None
    if gtf_file.endswith(".gz") : infile = gzip.open(gtf_file,"rt")
    else: infile =open(gtf_file,"r")
    for line in infile:
        if line.startswith("#"): continue 
        linesplit = line.split("\t")
        if linesplit[2]=="gene" and "gene_name \""+gene_name+"\"" in linesplit[8]: 
            for x in linesplit[8].split(";"):
                if "gene_id" in x:
                    gene_id = x[x.find("\"")+1:-1]
    infile.close()
    return gene_id

def find_gene_name(gene_id,gtf_file):
    gene_id=None
    if gtf_file.endswith(".gz") : infile = gzip.open(gtf_file,"rt")
    else: infile =open(gtf_file,"r")
    for line in infile:
        if line.startswith("#"): continue 
        linesplit = line.split("\t")
        if linesplit[2]=="gene" and "gene_id\""+gene_id+"\"" in linesplit[8]: 
            for x in linesplit[8].split(";"):
                if "gene_name" in x:
                    gene_name = x[x.find("\"")+1:-1]
    infile.close()
    return gene_name