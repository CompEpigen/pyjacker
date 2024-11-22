import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os 
from tqdm import tqdm
import multiprocessing
import sys

#from plot import plot_WGS
from figeno import figeno_make
from pyjacker.genes import find_gene_id, find_gene_name

def format_float(x):
    x=float(x)
    if x<0.001:
        return "{:.2E}".format(x)
    else:
        return "{:.4f}".format(x)
    


def generate_main_report(df_result,outdir,n_events=400,filter_monoallelic=False):
    df_result = df_result.copy(deep=True)
    df_result = df_result.loc[df_result["score"]>=0,:].reset_index(drop=True)
    n=0
    if "FDR" in df_result:
        while n<df_result.shape[0] and df_result.loc[n,"FDR"]<=0.8:
            n+=1
    else:
        while n<df_result.shape[0] and df_result.loc[n,"score"]>=0:
            n+=1
    n_events = min(n_events,n)

    with open(os.path.join(outdir,"report.html"),"w") as outfile:
        tmp = outfile.write("<!DOCTYPE html>\n<html lang=\"en\"><head>\n<meta charset=\"UTF-8\">\n<meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n")
        tmp = outfile.write("<title>Enhancer hijacking candidates</title>\n<link rel=\"stylesheet\" href=\"https://cdn.datatables.net/1.10.22/css/jquery.dataTables.min.css\">\n")
        tmp = outfile.write("<script src=\"https://code.jquery.com/jquery-3.5.1.js\"></script>\n<script src=\"https://cdn.datatables.net/1.10.22/js/jquery.dataTables.min.js\"></script>\n")
        
        tmp = outfile.write("<style>\nbody {padding: 0px; margin: 0; font-family: Verdana, Geneva, Tahoma, sans-serif;}\n")
        tmp = outfile.write("tr {transition: all .2s ease-in; cursor: pointer; color: #2c3e50;}\nth, td {padding: 12px; text-align: left; border-bottom: 1px solid #ddd;}\n")
        tmp = outfile.write("#header {background-color: #2c3e50; color: #fff;}\n h1 {font-weight: 600; text-align: center; color: #2c3e50; padding: 10px 0px;}\n")
        #tmp = outfile.write("tr:hover {background-color: #f5f5f5; transform: scale(1.02); box-shadow: 2px 2px 12px rgba(0, 0, 0, 0.2), -1px -1px 8px rgba(0, 0, 0, 0.2);}\n ")
        tmp = outfile.write("@media only screen and (max-width: 768px) {table {width: 90%;}}\n.details tr { line-height: 15px; }\n")
        tmp = outfile.write(".details th, td {padding: 5px; text-align: left; border-bottom: 1px solid #ddd; font-family:'Courier New', Courier, monospace;}\n")
        tmp = outfile.write("img{height: 50%;width: 50%;object-fit: contain;}\n")
        tmp = outfile.write("a, a:visited{color: blue;}\n")
        tmp = outfile.write("</style>\n<script>\n$(document).ready(function() {\n$('#cosmic').DataTable();});\n</script>\n</head>\n")
        tmp = outfile.write("<h1>Candidate enhancer hijacking events</h1><table class=\"details\">\n")
        tmp = outfile.write("<p style=\"margin-top:2.5em\"> </p>\n")
        tmp = outfile.write("<p style=\"margin-top:2.5em\"> </p>\n")
        tmp = outfile.write("<p style=\"margin-top:2.5em\"> </p>\n")
        tmp = outfile.write("<div style=\"width:80%\">\n")
        tmp = outfile.write("<table id=\"cosmic\" class=\"display\" >\n")
        tmp = outfile.write("<thead><tr id=\"header\"><th>Rank</th><th>FDR</th><th>Score</th><th>Gene</th><th>Chr</th><th>Start</th><th>End</th><th>Sample</th><th>Distance to breakpoint</th><th>Fusion</th><th>n_SNPs</th><th>OHE score</th><th>ASE score</th><th>Enhancer score</th></tr></thead>\n")
        tmp = outfile.write("<tbody>\n")
        rank=1
        for i in range(min(df_result.shape[0],n_events)):
            if (not filter_monoallelic) or df_result.loc[i,"ASE_score"]>=0:
                if "FDR" in df_result.columns: FDR = format_float(df_result.loc[i,"FDR"])
                else: FDR = "."
                tmp = outfile.write("<tr><td>" +str(rank) + "</td><td>"+FDR+"</td><td>"+format_float(df_result.loc[i,"score"])+"</td><td>"+df_result.loc[i,"gene_name"]+"</td><td>"+ df_result.loc[i,"chr"]+ "</td><td>"\
                    +str(df_result.loc[i,"start"]) + "</td><td>" +str(df_result.loc[i,"end"]) + "</td><td>" \
                    + "<a target=\"_blank\" href=\"Data/html/"+df_result.loc[i,"gene_name"]+"_"+df_result.loc[i,"sample"]+".html\">" +df_result.loc[i,"sample"] + "</a>"+"</td><td>"
                    + f'{df_result.loc[i,"distance_to_breakpoint"]:,}'+ "</td><td>"+ df_result.loc[i,"fusion"] +"</td><td>"+ str(df_result.loc[i,"n_SNPs"])+ "</td><td>"+ format_float(df_result.loc[i,"OHE_score"])+"</td><td>"+ format_float(df_result.loc[i,"ASE_score"])+"</td><td>"+ format_float(df_result.loc[i,"enhancer_score"])+"</td></tr>\n")
                rank+=1
        tmp = outfile.write("</table>\n</tbody>\n</div>\n</body>\n</html>\n")
        tmp = outfile.write("")
        tmp = outfile.write("<script>$('#cosmic').dataTable({paging: false, info: false});</script>")


# Generate report for individual enhancer hijacking events

# make boxplot


def get_df_SNPs(sample,gene,ase_dir):
    df_ase = pd.read_csv(os.path.join(ase_dir,sample+".tsv"),sep="\t",dtype={"contig":str})
    df_ase["contig"] = [x.lstrip("chr") for x in df_ase["contig"]]
    df_ase = df_ase.loc[df_ase["contig"]==gene.chr]
    df_ase = df_ase.loc[(df_ase["position"]>=gene.start) & (df_ase["position"]<=gene.end)]
    return df_ase


def write_SNPs(sample,gene,ase_dir,outfile):
    df_ase = get_df_SNPs(sample,gene,ase_dir)
    for i in df_ase.index:
        tmp = outfile.write("<tr><td>" + df_ase.loc[i,"variantID"]+ "</td><td>" + df_ase.loc[i,"contig"] + "</td><td>"+str(df_ase.loc[i,"position"])+"</td><td>" \
            + df_ase.loc[i,"refAllele"]+ "</td><td>"+df_ase.loc[i,"altAllele"] + "</td><td>"+ str(df_ase.loc[i,"refCount"]) + "</td><td>" + str(df_ase.loc[i,"altCount"])+"</td></tr>\n")



def generate_individual_reports(df_result,df_TPM,breakpoints,CNAs,genes,ase_dir,ase_dna_dir,gtf_file,outdir,cytobands,df_TPM_normal=None,n_events=-1,image_format="png",image_dpi=200,n_threads=4):
    df_result = df_result.copy(deep=True)
    df_result = df_result.loc[df_result["score"]>=0,:].reset_index(drop=True)
    os.makedirs(os.path.join(outdir,"Data/html"), exist_ok=True)
    os.makedirs(os.path.join(outdir,"Data/Figures/Expression/"), exist_ok=True)
    os.makedirs(os.path.join(outdir,"Data/Figures/chr_plots"), exist_ok=True)
    os.makedirs(os.path.join(outdir,"Data/Figures/ASE/"), exist_ok=True)

    if n_events<0: n_events=df_result.shape[0]
    else: n_events=min(df_result.shape[0],n_events)

    # Data that is copied for each process.
    data_dic={"df_result":df_result,"breakpoints":breakpoints,"CNAs":CNAs,"genes":genes,
            "df_TPM":df_TPM,"df_TPM_normal":df_TPM_normal, "ase_dir":ase_dir,"ase_dna_dir":ase_dna_dir,
            "gtf_file":gtf_file,"cytobands":cytobands,
            "outdir":outdir,"image_format":image_format, "image_dpi":image_dpi}

    # Multiprocessing
    def init_worker():
        global data
        data = data_dic
    with multiprocessing.Pool(n_threads, initializer=init_worker) as pool:
        with tqdm(total=n_events, file=sys.stdout) as pbar:
            for _ in pool.imap(generate_report, range(n_events)):
                pbar.update() 

def generate_report(i):
    """Generate the report for the index i of df_result"""
    df_result=data["df_result"]
    genes=data["genes"]
    df_TPM=data["df_TPM"]
    breakpoints = data["breakpoints"]
    CNAs=data["CNAs"]

    sample = df_result.loc[i,"sample"]
    gene_name = df_result.loc[i,"gene_name"]
    gene_full = genes[df_result.loc[i,"gene_id"]]
    chr = df_result.loc[i,"chr"]
    plot_expression(os.path.join(data["outdir"],"Data/Figures/Expression/"+gene_name+"_"+sample+"."+data["image_format"]),
                    df_TPM,gene_id=gene_full.gene_id,gene_name=gene_full.gene_name,samples={sample:sample},colors={sample:"#27ae60"},image_dpi=data["image_dpi"])
    if CNAs is not None:
        d_SV={"chr1":[],"pos1":[],"chr2":[],"pos2":[],"color":[]}
        for bp in breakpoints[chr]:
            if bp.sample==sample and bp.pos2!="NA" and bp.chr1==chr:
                if bp.chr2==chr and bp.pos1>bp.pos2: continue
                d_SV["chr1"].append(bp.chr1)
                d_SV["pos1"].append(bp.pos1)
                d_SV["chr2"].append(bp.chr2)
                d_SV["pos2"].append(bp.pos2)
                if bp.chr2!=bp.chr1: d_SV["color"].append("#27ae60")
                elif bp.orientation1=="-" and bp.orientation2=="+": d_SV["color"].append("#4a69bd")
                elif bp.orientation1=="+" and bp.orientation2=="-": d_SV["color"].append("#e55039")
                elif bp.orientation1=="-" and bp.orientation2=="-": d_SV["color"].append("purple")
                else: d_SV["color"].append("#e58e26")
        df_SVs = pd.DataFrame(d_SV)
        if data["cytobands"] is not None:
            config_chr={"general":{"reference":"custom","layout":"horizontal","genes_file":data["gtf_file"],"cytobands_file":data["cytobands"]}}
            config_chr["output"] = {"file":os.path.join(data["outdir"],"Data/Figures/chr_plots/"+sample+"_"+gene_name+"."+data["image_format"]),"dpi":data["image_dpi"],"width":183.0}
            config_chr["regions"] = [{"chr":chr}]
            config_chr["tracks"] = [{"type":"sv","height": 10.0,"margin_above": 0.0,"bounding_box": True,"df_SVs":df_SVs},
                                    {"type":"copynumber","height": 30.0,"margin_above": 0.0,"bounding_box": True,"CNAs":CNAs[sample],"genes":gene_name,
                                    "min_cn":None,"max_cn":None},
                                    {"type":"chr_axis","height": 10.0,"unit":"Mb","margin_above":0.0}]
            figeno_make(config_chr)
    ase_plot_exists=False
    if (data["ase_dir"] is not None):
        ase_file=os.path.join(data["ase_dir"],sample+".tsv")
        vcf_DNA= None
        if data["ase_dna_dir"] is not None: vcf_DNA = os.path.join(data["ase_dna_dir"],sample+".vcf.gz")

        # ASE plot
        config_ase={"general":{"reference":"custom","layout":"horizontal","genes_file":data["gtf_file"]}}
        config_ase["output"] = {"file":os.path.join(data["outdir"],"Data/Figures/ASE/"+gene_name+"_"+sample+"."+data["image_format"]),"dpi":data["image_dpi"],"width":100.0}
        config_ase["regions"] = [{"chr":chr,"start":gene_full.start,"end":gene_full.end}]
        config_ase["tracks"] = [{"type":"ase","height":60,"margin_above":1.5,"file":ase_file,"vcf_DNA":vcf_DNA,"grid":False},
                                {"type":"genes","height":12.0,"margin_above":-2,"bounding_box": False,"genes":gene_name},
                                {"type":"chr_axis","height": 10.0}]
        figeno_make(config_ase)
        ase_plot_exists=True
    with open(os.path.join(data["outdir"],"Data/html/"+gene_name+"_"+sample+".html"),"w") as outfile:
        tmp = outfile.write("<!DOCTYPE html>\n<html lang=\"en\"><head>\n<meta charset=\"UTF-8\">")
        tmp = outfile.write("<div style=\"text-align:center\">\n")
        tmp = outfile.write("<meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n")
        tmp = outfile.write("<title>"+gene_name+" " + sample+"</title>\n")
        tmp = outfile.write("<link rel=\"stylesheet\" href=\"https://cdn.datatables.net/1.10.22/css/jquery.dataTables.min.css\">\n")
        tmp = outfile.write("<script src=\"https://code.jquery.com/jquery-3.5.1.js\"></script>\n<script src=\"https://cdn.datatables.net/1.10.22/js/jquery.dataTables.min.js\"></script>\n")
        tmp = outfile.write("<style>\nbody {padding: 0px; margin: 0; font-family: Verdana, Geneva, Tahoma, sans-serif;}\n")
        tmp = outfile.write("tr {transition: all .2s ease-in; cursor: pointer; color: #2c3e50;}\nth, td {padding: 12px; text-align: left; border-bottom: 1px solid #ddd;}\n")
        tmp = outfile.write("#header {background-color: #2c3e50; color: #fff;}\n")
        tmp = outfile.write("h1,h2,h3 {font-weight: 600; text-align: center; color: #2c3e50; padding: 10px 0px;}\n")
        tmp = outfile.write("tr:hover {background-color: #f5f5f5; transform: scale(1.02); box-shadow: 2px 2px 12px rgba(0, 0, 0, 0.2), -1px -1px 8px rgba(0, 0, 0, 0.2);}\n ")
        tmp = outfile.write("@media only screen and (max-width: 768px) {table {width: 90%;}}\n.details tr { line-height: 15px; }\n")
        tmp = outfile.write(".details th, td {padding: 5px; text-align: left; border-bottom: 1px solid #ddd; font-family:'Courier New', Courier, monospace;}\n")
        tmp = outfile.write("img{height: 50%;width: 50%;object-fit: contain;}")
        tmp = outfile.write("</style>\n<script>\n$(document).ready(function() {\n$('#cosmic').DataTable();});\n</script>\n</head>\n")
        #tmp = outfile.write("<style>\nimg{height: 20%;width: 80%;object-fit: contain;}\n</style>\n")
        #tmp = outfile.write("<style>img:hover {transform:scale(1.5);-ms-transform:scale(1.5); /* IE 9 */-moz-transform:scale(1.5); /* Firefox */-webkit-transform:scale(1.5); /* Safari and Chrome */-o-transform:scale(1.5); /* Opera */}</style>")
        tmp = outfile.write("<h1>Enhancer hijacking of " + gene_name+" in sample " + sample + "</h1>\n")
        tmp = outfile.write("<p>"+gene_name+" chr"+chr+":" +f'{gene_full.start:,}'+"-"+f'{gene_full.end:,}' + " (hg19) </p>")
        tmp = outfile.write("<img style=\"height:auto;width:30%;\" src=\"../Figures/Expression/"+gene_name+"_"+sample+"."+data["image_format"]+"\" alt=\"Expression plot\">\n")
        if (CNAs is not None) and (data["cytobands"] is not None):
            tmp = outfile.write("<h2> Copy number and SV plot</h2>\n")
            tmp = outfile.write("<img style=\"height:auto;width:70%;\" src=\"../Figures/chr_plots/"+sample+"_"+gene_name+"."+data["image_format"]+"\" alt=\"Copy number plot\">\n")
        tmp = outfile.write("")
        tmp = outfile.write("")
        tmp = outfile.write("")
        if data["ase_dir"] is not None and ase_plot_exists:
            tmp = outfile.write("<h2>Allele specific expression</h2>\n")
            tmp = outfile.write("<img style=\"height:auto;width:500px;\" src=\"../Figures/ASE/"+gene_name+"_"+sample+"."+data["image_format"]+"\" alt=\"ASE plot\">\n")

        # SVs
        tmp = outfile.write("<h2>Structural variants on chromosome "+chr+"</h2>\n")
        tmp = outfile.write("<p style=\"margin-top:2.5em\"> </p>\n")
        tmp = outfile.write("<table id=\"sv\" class=\"display\" style=\"width:60%\">\n")
        tmp = outfile.write("<thead><tr id=\"header\"><th>Chr1</th><th>Pos1</th><th>Chr2</th><th>Pos2</th></tr></thead>\n")
        tmp = outfile.write("<tbody>\n")
        for breakpoint in breakpoints[chr]:
            if breakpoint.sample == sample:
                if breakpoint.pos2!="NA":
                    tmp = outfile.write("<tr><td>" +str(breakpoint.chr1) + "</td><td>"+f'{breakpoint.pos1:,}'+"</td><td>"+ str(breakpoint.chr2)+ "</td><td>"+f'{breakpoint.pos2:,}'+"</td></tr>\n")
                else:
                    tmp = outfile.write("<tr><td>" +str(breakpoint.chr1) + "</td><td>"+f'{breakpoint.pos1:,}'+"</td><td>NA</td><td>NA</td></tr>\n")
        tmp = outfile.write("</table>\n</tbody>\n")
        tmp = outfile.write("<script>$('#sv').dataTable({searching: false, paging: false, info: false});</script>")

        # CNAs
        if CNAs is not None:
            tmp = outfile.write("<h2>Copy number alterations on chromosome "+chr+"</h2>\n")
            tmp = outfile.write("<p style=\"margin-top:2.5em\"> </p>\n")
            tmp = outfile.write("<table id=\"cnv\" class=\"display\" style=\"width:60%\">\n")
            tmp = outfile.write("<thead><tr id=\"header\"><th>Chr</th><th>Start</th><th>End</th><th>Copy Number</th></tr></thead>\n")
            tmp = outfile.write("<tbody>\n")
            if sample in CNAs and chr in CNAs[sample]:
                for CNA in CNAs[sample][chr]:
                    tmp = outfile.write("<tr><td>" +chr + "</td><td>"+f'{int(CNA[0]):,}'+"</td><td>"+ f'{int(CNA[1]):,}'+ "</td><td>"+str(CNA[2])+"</td></tr>\n")
                
            tmp = outfile.write("</table>\n</tbody>\n")
            tmp = outfile.write("<script>$('#cnv').dataTable({searching: false, paging: false, info: false});</script>")

        # SNPs
        if data["ase_dir"] is not None:
            tmp = outfile.write("<h2>SNPs in "+gene_name+"</h2>\n")
            tmp = outfile.write("<p style=\"margin-top:2.5em\"> </p>\n")
            tmp = outfile.write("<table id=\"snv\" class=\"display\" style=\"width:60%\">\n")
            tmp = outfile.write("<thead><tr id=\"header\"><th>ID</th><th>Chr</th><th>Pos</th><th>Ref</th><th>Alt</th><th>Ref reads</th><th>Alt reads</th></tr></thead>\n")
            tmp = outfile.write("<tbody>\n")
            write_SNPs(sample,gene_full,data["ase_dir"],outfile)
            tmp = outfile.write("</table>\n</tbody>\n")
            tmp = outfile.write("<script>$('#snv').dataTable({searching: false, paging: false, info: false});</script>")


        # Expression in normal samples
        if data["df_TPM_normal"] is not None:
            os.makedirs(os.path.join(data["outdir"],"Data/Figures/Expression_normal/"), exist_ok=True)
            other_plot_exists = plot_normal_samples(data["df_TPM_normal"],gene_full,os.path.join(data["outdir"],"Data/Figures/Expression_normal/"+gene_name+"_normal."+data["image_format"]),image_dpi=data["image_dpi"])
            if other_plot_exists:
                tmp = outfile.write("<h2>Expression of "+gene_name+" in normal samples</h2>\n")
                tmp = outfile.write("<img style=\"height:auto;width:40%;\" src=\"../Figures/Expression_normal/"+gene_name+"_normal."+data["image_format"]+"\" alt=\"Expression plot in normal samples\">\n")

            
        tmp = outfile.write("<p style=\"margin-top:2.5em\"> </p>\n")
        
        tmp = outfile.write("</div>\n")


def plot_normal_samples(df,gene,output_filename,image_dpi=200):
    if df is None: return False
    if not gene.gene_id in df.index: return False
    matplotlib.rcParams.update({'font.size': 20})
    plt.figure()
    

    plt.barh(list(df.columns),df.loc[gene.gene_id,:])
    plt.tick_params(axis='y', labelsize=10)
    plt.xlim(0,max(np.max(df.loc[gene.gene_id,:]),10))

    plt.savefig(output_filename,bbox_inches="tight",pad_inches=0.1,dpi=image_dpi)

    plt.cla() 
    plt.clf() 
    plt.close('all')
    return True

def plot_expression(output_filename,df_TPM,gene_name=None,gene_id=None,samples=None,colors=None,gtf_file=None,plot_legend=False,image_dpi=200):
    """
    samples: dictionary: group to list of samples in that group
    colors: dictionary: group to color
    """
    if (gene_name is None) and gene_id is None: return False
    if gene_id is None:
        if gtf_file is None: return False
        gene_id = find_gene_id(gene_name=gene_name,gtf_file=gtf_file)
        if gene_id is None: return False
    if gene_name is None: 
        if gtf_file is not None:
            gene_name = find_gene_name(gene_id,gtf_file=gtf_file)
        else: gene_name = ""

    matplotlib.rcParams.update({'font.size': 20})
    fig,ax = plt.subplots()
    values = df_TPM.loc[gene_id,:].sort_values()
    c = []
    for sample in values.index:
        group="other"
        if samples is not None:
            for g in samples:
                if sample in samples[g]: group = g
        color="#999999"
        if group in colors: color=colors[g]
        c.append(color)
    
    plt.scatter(range(len(values)),values,c=c,s=100)
    plt.ylabel(gene_name+" expression (TPM)")
    plt.xlabel("Samples")
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off

    #plt.text(i_sample-28,values[i_sample] - values[-1]*0.025,sample)

    if colors is not None and plot_legend:
        legend_elements=[]
        for g in colors:
            legend_elements.append(Line2D([0], [0], marker='o', color='w', label=g, markerfacecolor=colors[g], markersize=15))
        legend_elements.append(Line2D([0], [0], marker='o', color="w", label="other", markerfacecolor="black", markersize=15))
        ax.legend(handles=legend_elements, loc='upper left')
    
    plt.savefig(output_filename,bbox_inches="tight",pad_inches=0.1,dpi=image_dpi)
    plt.cla() 
    plt.clf() 
    plt.close('all')