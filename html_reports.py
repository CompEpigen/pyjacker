import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import os 
import vcfpy

import config


def format_float(x):
    x=float(x)
    if x<0.001:
        return "{:.2E}".format(x)
    else:
        return "{:.4f}".format(x)

def generate_main_report(df_result,outdir,n_events=400,filter_monoallelic=False):
    n=0
    while n<df_result.shape[0] and df_result.loc[n,"padj"]<=0.05:
        n+=1
    n_events = max(n_events,n)

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
        tmp = outfile.write("<thead><tr id=\"header\"><th>Rank</th><th>Score</th><th>Gene</th><th>Chr</th><th>Start</th><th>End</th><th>Sample</th><th>Distance to breakpoint</th><th>Fusion</th><th>pval outlier exp</th><th>n_SNPs</th><th>pval monoallelic</th><th>llr monoallelic</th><th>Enhancers</th></tr></thead>\n")
        tmp = outfile.write("<tbody>\n")
        rank=1
        for i in range(min(df_result.shape[0],n_events)):
            if (not filter_monoallelic) or df_result.loc[i,"llr_monoallelic"]>=0:
                tmp = outfile.write("<tr><td>" +str(rank) + "</td><td>"+format_float(df_result.loc[i,"score"])+"</td><td>"+df_result.loc[i,"gene_name"]+"</td><td>"+ df_result.loc[i,"chr"]+ "</td><td>"\
                    +str(df_result.loc[i,"start"]) + "</td><td>" +str(df_result.loc[i,"end"]) + "</td><td>" \
                    + "<a target=\"_blank\" href=\"Data/html/"+df_result.loc[i,"gene_name"]+"_"+df_result.loc[i,"sample"]+".html\">" +df_result.loc[i,"sample"] + "</a> </td><td>"\
                    + f'{df_result.loc[i,"distance_to_breakpoint"]:,}'+ "</td><td>"+ df_result.loc[i,"fusion"]+ "</td><td>"+format_float(df_result.loc[i,"pval"])+"</td><td>"+ str(df_result.loc[i,"n_SNPs"])+ "</td><td>"+ format_float(df_result.loc[i,"pval_monoallelic"])+ "</td><td>"+ format_float(df_result.loc[i,"llr_monoallelic"])+"</td><td>"+ df_result.loc[i,"enhancers"]+"</td></tr>\n")
                rank+=1
        tmp = outfile.write("</table>\n</tbody>\n</div>\n</body>\n</html>\n")
        tmp = outfile.write("")
        tmp = outfile.write("<script>$('#cosmic').dataTable({paging: false, info: false});</script>")


# Generate report for individual enhancer hijacking events

# make boxplot

def rank_order_plot(sample,gene,output_filename):
    matplotlib.rcParams.update({'font.size': 20})
    plt.figure()

    values = []
    for col in config.samples: #df_TPM.columns
        values.append(config.df_TPM.loc[gene.gene_id,col])
    values = sorted(values)

    colors = ["red" if x ==config.df_TPM.loc[gene.gene_id,sample] else "black" for x in values]

    i_sample = -1
    for i in range(len(colors)):
        if config.df_TPM.loc[gene.gene_id,sample]==values[i]:
            i_sample=i
    plt.scatter(range(len(values)),values,c=colors)
    plt.ylabel(gene.gene_name+" expression (TPM)")
    plt.xlabel("Samples")
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off

    #plt.text(i_sample-28,values[i_sample] - values[-1]*0.025,sample)
    plt.savefig(output_filename,bbox_inches="tight",pad_inches=0.1,dpi=300)
    plt.cla() 
    plt.clf() 
    plt.close('all')

def rank_order_plot_reference(gene,output_filename):
    if not "df_TPM_reference" in dir(config): return False
    matplotlib.rcParams.update({'font.size': 20})
    plt.figure()

    if not gene.gene_id in config.df_TPM_reference.index: return False

    values = []
    for col in config.df_TPM_reference.columns:
        values.append(config.df_TPM_reference.loc[gene.gene_id,col])
    values = sorted(values)

    colors = ["black" for x in values]

    plt.scatter(range(len(values)),values,c=colors)
    plt.ylabel(gene.gene_name+" expression (TPM)")
    plt.xlabel("Samples")
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off

    plt.savefig(output_filename,bbox_inches="tight",pad_inches=0.1,dpi=300)
    plt.ylim(0,max(np.max(values),10))

    plt.cla() 
    plt.clf() 
    plt.close('all')
    return True

def rank_order_plot_other(df,gene,output_filename):
    #if not "df_TPM_other" in dir(config): return False
    matplotlib.rcParams.update({'font.size': 20})
    plt.figure()

    if not gene.gene_id in df.index: return False

    values = []
    for col in df.columns:
        if df.loc[gene.gene_id,col]==df.loc[gene.gene_id,col] and df.loc[gene.gene_id,col]!="nan":
            values.append(df.loc[gene.gene_id,col])
    values = sorted(values)

    colors = ["black" for x in values]

    plt.scatter(range(len(values)),values,c=colors)
    plt.ylabel(gene.gene_name+" expression (TPM)")
    plt.xlabel("Samples")
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off

    plt.savefig(output_filename,bbox_inches="tight",pad_inches=0.1,dpi=300)
    if np.max(values)!=np.max(values):
        print(values)
    plt.ylim(0,max(np.max(values),10))

    plt.cla() 
    plt.clf() 
    plt.close('all')
    return True

def plot_normal_samples(df,gene,output_filename):
    #if not "df_TPM_other" in dir(config): return False
    matplotlib.rcParams.update({'font.size': 20})
    plt.figure()

    if not gene.gene_id in df.index: return False

    plt.barh(list(df.columns),df.loc[gene.gene_id,:])
    plt.tick_params(axis='y', labelsize=10)
    plt.xlim(0,max(np.max(df.loc[gene.gene_id,:]),10))

    plt.savefig(output_filename,bbox_inches="tight",pad_inches=0.1,dpi=300)

    plt.cla() 
    plt.clf() 
    plt.close('all')
    return True


def plot_SNP_allelefrac(SNPs,filename):
    y1=[]
    y2=[]
    x=[]
    for snp in SNPs:
        x.append(snp.chr+":"+str(snp.pos))
        y1.append(snp.RO)
        y2.append(snp.AD)
    plt.bar(x,y1,label="Reference counts")
    plt.bar(x,y2,label="Alternate counts",bottom=y1)
    plt.legend()
    plt.tick_params(axis='x',which='both',bottom=False,labelbottom=False) 
    plt.savefig(filename,bbox_inches="tight",pad_inches=0.1,dpi=300)
    plt.cla() 
    plt.clf() 
    plt.close('all')



def generate_individual_reports(df_result,outdir,n_events=600):
    n=0
    while  n<df_result.shape[0] and df_result.loc[n,"padj"]<=0.05:
        n+=1
    n_events = max(n_events,n)
    for i in range(min(df_result.shape[0],n_events)): #range(df_result.shape[0])
        sample = df_result.loc[i,"sample"]
        gene_name = df_result.loc[i,"gene_name"]
        gene_full = config.data.gene_by_id(df_result.loc[i,"gene_id"])
        chr = df_result.loc[i,"chr"]
        os.makedirs(os.path.join(outdir,"Data/Figures/Expression/"), exist_ok=True)
        rank_order_plot(sample,gene_full,os.path.join(outdir,"Data/Figures/Expression/"+gene_name+"_"+sample+".png"))
        os.makedirs(os.path.join(outdir,"Data/Figures/Expression_reference/"), exist_ok=True)
        reference_plot_exists = rank_order_plot_reference(gene_full,os.path.join(outdir,"Data/Figures/Expression_reference/"+gene_name+"_reference.png"))
        os.makedirs(os.path.join(outdir,"Data/Figures/SNPs/"), exist_ok=True)
        plot_SNP_allelefrac(config.SNPs[sample][gene_full.gene_id],os.path.join(outdir,"Data/Figures/SNPs/"+gene_name+"_"+sample+".png"))
        os.makedirs(os.path.join(outdir,"Data/html"), exist_ok=True)

        with open(os.path.join(outdir,"Data/html/"+gene_name+"_"+sample+".html"),"w") as outfile:
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
            tmp = outfile.write("<img style=\"height:auto;width:40%;\" src=\"../Figures/Expression/"+gene_name+"_"+sample+".png\" alt=\"Expression plot\">\n")
            tmp = outfile.write("<h2> Copy number plot</h2>\n")
            tmp = outfile.write("<img style=\"height:auto;width:1700px;\" src=\"../Figures/linearplots/"+sample+"/"+sample+"_chr"+chr+".png\" alt=\"Copy number plot\">\n")
            tmp = outfile.write("<h2>Circos plot</h2>\n")
            tmp = outfile.write("<img style=\"height:auto;width:1000px;\" src=\"../Figures/circos/circos_"+sample+".png\" alt=\"Circos plot\">\n")
            tmp = outfile.write("")
            tmp = outfile.write("")
            tmp = outfile.write("")

            # SVs
            tmp = outfile.write("<h2>Structural variants on chromosome "+chr+"</h2>\n")
            tmp = outfile.write("<p style=\"margin-top:2.5em\"> </p>\n")
            tmp = outfile.write("<table id=\"sv\" class=\"display\" style=\"width:60%\">\n")
            tmp = outfile.write("<thead><tr id=\"header\"><th>Chr1</th><th>Pos1</th><th>Chr2</th><th>Pos2</th></tr></thead>\n")
            tmp = outfile.write("<tbody>\n")
            for breakpoint in config.breakpoints[chr]:
                if breakpoint.sample == sample:
                    if breakpoint.pos2!="NA":
                        tmp = outfile.write("<tr><td>" +str(breakpoint.chr) + "</td><td>"+f'{breakpoint.pos:,}'+"</td><td>"+ str(breakpoint.chr2)+ "</td><td>"+f'{breakpoint.pos2:,}'+"</td></tr>\n")
                    else:
                        tmp = outfile.write("<tr><td>" +str(breakpoint.chr) + "</td><td>"+f'{breakpoint.pos:,}'+"</td><td>NA</td><td>NA</td></tr>\n")
            tmp = outfile.write("</table>\n</tbody>\n")
            tmp = outfile.write("<script>$('#sv').dataTable({searching: false, paging: false, info: false});</script>")

            # CNAs
            if "CNAs" in dir(config):
                tmp = outfile.write("<h2>Copy number alterations on chromosome "+chr+"</h2>\n")
                tmp = outfile.write("<p style=\"margin-top:2.5em\"> </p>\n")
                tmp = outfile.write("<table id=\"cnv\" class=\"display\" style=\"width:60%\">\n")
                tmp = outfile.write("<thead><tr id=\"header\"><th>Chr</th><th>Start</th><th>End</th><th>Copy Number</th></tr></thead>\n")
                tmp = outfile.write("<tbody>\n")
                for CNA in config.CNAs[sample][chr]:
                    tmp = outfile.write("<tr><td>" +CNA.chr + "</td><td>"+f'{int(CNA.start):,}'+"</td><td>"+ f'{int(CNA.end):,}'+ "</td><td>"+str(CNA.cn)+"</td></tr>\n")
                    
                tmp = outfile.write("</table>\n</tbody>\n")
                tmp = outfile.write("<script>$('#cnv').dataTable({searching: false, paging: false, info: false});</script>")

            # SNPs 
            if "SNV_RNA_dir" in dir(config):
                tmp = outfile.write("<h2>SNPs in "+gene_name+"</h2>\n")
                tmp = outfile.write("<p style=\"margin-top:2.5em\"> </p>\n")
                tmp = outfile.write("<table id=\"snv\" class=\"display\" style=\"width:60%\">\n")
                tmp = outfile.write("<thead><tr id=\"header\"><th>Chr</th><th>Pos</th><th>Ref</th><th>Alt</th><th>Ref reads</th><th>Alt reads</th><th>Gnomad Freq</th></tr></thead>\n")
                tmp = outfile.write("<tbody>\n")

                for snp in config.SNPs[sample][gene_full.gene_id]:
                    tmp = outfile.write("<tr><td>" + snp.chr + "</td><td>"+str(snp.pos)+"</td><td>"+ snp.ref+ "</td><td>"\
                                + snp.alt + "</td><td>"+ str(snp.RO) + "</td><td>" + str(snp.AD)\
                                +"</td><td>" + str(snp.popAF) +"</td></tr>\n")
                
                tmp = outfile.write("</table>\n</tbody>\n")
                tmp = outfile.write("<script>$('#snv').dataTable({searching: false, paging: false, info: false});</script>")

            #Expression in reference
            if reference_plot_exists:
                tmp = outfile.write("<h2>Expression of "+gene_name+" in the reference data</h2>\n")
                tmp = outfile.write("<img style=\"height:auto;width:40%;\" src=\"../Figures/Expression_reference/"+gene_name+"_reference.png\" alt=\"Expression plot in reference\">\n")

            # Expression in other samples
            os.makedirs(os.path.join(outdir,"Data/Figures/Expression_other/"), exist_ok=True)
            other_plot_exists = rank_order_plot_other(config.df_TPM_other,gene_full,os.path.join(outdir,"Data/Figures/Expression_other/"+gene_name+"_other.png"))
            if other_plot_exists:
                tmp = outfile.write("<h2>Expression of "+gene_name+" in other samples</h2>\n")
                tmp = outfile.write("<img style=\"height:auto;width:40%;\" src=\"../Figures/Expression_other/"+gene_name+"_other.png\" alt=\"Expression plot in other samples\">\n")

            # Expression in normal samples
            os.makedirs(os.path.join(outdir,"Data/Figures/Expression_normal/"), exist_ok=True)
            other_plot_exists = plot_normal_samples(config.df_TPM_normal,gene_full,os.path.join(outdir,"Data/Figures/Expression_normal/"+gene_name+"_normal.png"))
            if other_plot_exists:
                tmp = outfile.write("<h2>Expression of "+gene_name+" in normal samples</h2>\n")
                tmp = outfile.write("<img style=\"height:auto;width:40%;\" src=\"../Figures/Expression_normal/"+gene_name+"_normal.png\" alt=\"Expression plot in normal samples\">\n")

            
            tmp = outfile.write("<p style=\"margin-top:2.5em\"> </p>\n")
            
            tmp = outfile.write("</div>\n")


