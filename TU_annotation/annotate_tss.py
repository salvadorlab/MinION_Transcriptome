from email.policy import default
import os
import pandas as pd
import argparse
import math


parser=argparse.ArgumentParser(prog="tss_anot",usage='%(prog)s[options]',description='anot each tss as gTSS (pTSS/sTSS); iTSS; asTSS; oTSS:')
parser.add_argument('--input', "-i",type=str,help='TSS files')
parser.add_argument('--gff', "-g",type=str,help='reference gene annotation')
parser.add_argument('--sample', "-s",type=str,help='(optional) this is the combined and clustered file produced using combine_samples.r and cluster_tss.py, will show which sample each tss is present in', default="")
parser.add_argument('--output', "-o",type=str,help='output file')
parser.add_argument('--isCombined', "-c",type=bool,help='indicate if the input tss file is obtained from a single sample or combined. Default is False', default=False)
args = parser.parse_args()

filter_tss=args.input
gff=args.gff
output=args.output
sample=args.sample
isCombined=bool(args.isCombined)


# filter_tss="/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/cdna_output/tss/Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_PolyATail_rna_filtered.linear.TSS.merged.tab"
# gff="/scratch/rx32940/minION/polyA_cDNA/map/genome/reference/GCF_000007685.1_ASM768v1_genomic.gff"
# output="/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/cdna.filtered20.TSS.tab"
# sample=""
# isCombined=False

# put all reads in a dict 
# with structure {chromI: {+:[(start, end)], lags:[(start, end)]}, chromII: {+:[(start, end)], lags:[(start, end)]}}
def get_genes(gff):
    all_genes={}
    with open(gff) as gf:
        for line in gf:
            line_list=line.split('\t')
            if line[:1] == "#":
                continue
            elif line_list[2] == "region":
                chrom=line.split('\t')[0]
                all_genes[chrom] = {}
                all_genes[chrom]["+"] = []
                all_genes[chrom]["-"] = []
            elif line_list[2] != "CDS": 
                gene_info=line_list[8].split(';')
                gene_id=gene_info["ID=" in gene_info].split("-")[1]
                chrom=line_list[0]
                strand=line_list[6]
                all_genes[chrom][strand].append((int(line_list[3]), int(line_list[4]), gene_id))       
    return all_genes
 
def get_tss(all_genes, filter_tss):
    all_genes=all_genes
    all_tss={}
    for chrom in all_genes.keys():
        all_tss[chrom]={}
        for strand in all_genes[chrom].keys():
            all_tss[chrom][strand]=[]
    with open(filter_tss) as ft:
        for line in ft.readlines()[1:]:
            line_list=line.strip('\n').split('\t')
            if isCombined:
                all_tss[line_list[0]][line_list[5]].append((int(line_list[1]), int(line_list[2]), line_list[3],int(line_list[4]),int(line_list[6]), int(line_list[7])))
            else:
                all_tss[line_list[0]][line_list[5]].append((int(line_list[1]), int(line_list[2]), line_list[3],int(line_list[4])))
    return all_tss


def get_tss_anot(all_genes,chrom ,strand, item):
    tss=int(item[0])
     # each tss loop through all genes, both lead and lag
    current_gene=""
    current_anot=""
    for gene in all_genes[chrom][strand]:
        if strand == "+": # on lead strand
            if tss > gene[1]:
                continue
            elif gene[0] - tss > 300:
                break
            else:
                if gene[0] - tss <=300 and (gene[0] -tss >=0):
                    if current_gene == "":
                        current_gene = gene[2]
                    elif gene[2] != current_gene.split(";")[-1]:
                        current_gene = current_gene + ";" + gene[2]
                    else:
                        current_gene = current_gene
                    current_anot = current_anot + ";" +"gTSS"
                if tss >= gene[0] and tss <= gene[1]:
                    if current_gene == "":
                        current_gene = gene[2]
                    elif gene[2] != current_gene.split(";")[-1]:
                        current_gene = current_gene + ";" + gene[2]
                    else:
                        current_gene = current_gene
                    current_anot = current_anot + ";"+ "iTSS" 
        else: # on lag strand
            if (tss < gene[0]):
                break
            else:
                if (tss - gene[1] <= 300) and (tss - gene[1] >=0):
                    if current_gene == "":
                        current_gene = gene[2]
                    elif gene[2] != current_gene.split(";")[-1]:
                        current_gene = current_gene + ";" + gene[2]
                    else:
                        current_gene = current_gene
                    current_anot = current_anot + ";"+"gTSS"
                if tss>= gene[0] and tss <= gene[1]:
                    if current_gene == "":
                        current_gene = gene[2]
                    elif gene[2] != current_gene.split(";")[-1]:
                        current_gene = current_gene + ";" + gene[2]
                    else:
                        current_gene = current_gene
                    current_anot = current_anot + ";"+"iTSS"      
    if strand == "+":
        ostrand = "-"
    else:
        ostrand = "+"
    for gene in all_genes[chrom][ostrand]:
        if ostrand == "+":
            if tss > gene[1]:
                continue
            if (gene[0] - tss <= 100):
                if current_gene == "":
                    current_gene = gene[2]
                elif gene[2] != current_gene.split(";")[-1]:
                    current_gene = current_gene + ";" + gene[2]
                else:
                    current_gene = current_gene
                current_anot = current_anot + ";"+"asTSS"
            if tss>= gene[0] and tss <= gene[1]:
                if current_gene == "":
                    current_gene = gene[2]
                elif gene[2] != current_gene.split(";")[-1]:
                    current_gene = current_gene + ";" + gene[2]
                else:
                    current_gene = current_gene
                current_anot = current_anot + ";"+"asTSS"
        else:
            if (tss < gene[0]):
                break
            else:
                if (tss -gene[1] <= 100) and (tss -gene[1] >=0):
                    if current_gene == "":
                        current_gene = gene[2]
                    elif gene[2] != current_gene.split(";")[-1]:
                        current_gene = current_gene + ";" + gene[2]
                    else:
                        current_gene = current_gene
                    current_anot = current_anot + ";"+"asTSS"
                if tss>= gene[0] and tss <= gene[1]:
                    if current_gene == "":
                        current_gene = gene[2]
                    elif gene[2] != current_gene.split(";")[-1]:
                        current_gene = current_gene + ";" + gene[2]
                    else:
                        current_gene = current_gene
                    current_anot = current_anot + ";"+"asTSS"
    if current_anot == "":
        return (*item, "oTSS", "")
    else:
        current_anot = current_anot.strip(";")
        return (*item, current_anot, current_gene)

def annotate_tss(all_genes, all_tss):
    all_genes=all_genes
    all_tss=all_tss
    tss_anot={}
    # loop through tss to get annotated
    for chrom in all_tss.keys():
        tss_anot[chrom]={}
        for strand in all_tss[chrom].keys():
            tss_anot[chrom][strand]=[]
            for item in all_tss[chrom][strand]:
                tss_anot[chrom][strand].append(get_tss_anot(all_genes, chrom,strand, item))

    return tss_anot

    
def annotate_gTSS(tss_anot):
    annotated_df=pd.DataFrame(columns=["start","end", "geneFrom", "cov", "numSamples", "clusters", "tss_type","tssGene","chrom", "strand"])
    if isCombined:
        colnames=["start","end", "geneFrom", "cov", "numSamples", "clusters", "tss_type", "tssGene"]
    else:
        colnames=["start","end", "geneFrom", "cov", "tss_type", "tssGene"]
    for chrom in tss_anot.keys():
        for strand in tss_anot[chrom].keys():
            cur_df = pd.DataFrame(tss_anot[chrom][strand], columns=colnames)
            cur_df["gTSS"] = cur_df.apply(lambda x: "gTSS" if "gTSS" in x["tss_type"] else "",axis=1)
            cur_df["chrom"]=chrom
            cur_df["strand"]=strand
            annotated_df = pd.concat([annotated_df, cur_df], ignore_index=True)
    annotated_df.loc[annotated_df["gTSS"] != "", "gTSS_gene"]=annotated_df.loc[annotated_df["gTSS"] != ""].apply(lambda x: x['tssGene'].split(";")[x["tss_type"].split(";").index("gTSS")] if "gTSS" in x["tss_type"] else "", axis=1)
    annotated_df.loc[annotated_df["gTSS"] != "", "max_cov"]=annotated_df.loc[annotated_df["gTSS"] != ""].groupby(['gTSS_gene'])['cov'].transform(max)
    annotated_df.loc[annotated_df["gTSS"] == "gTSS", "gTSS_anot"]=annotated_df.loc[annotated_df["gTSS"] != ""].apply(lambda x: "pTSS" if x["cov"] == x["max_cov"] else "sTSS",axis=1)
    annotated_df.loc[annotated_df["gTSS_anot"] == "pTSS", "pTSS_count"]=annotated_df.loc[annotated_df["gTSS_anot"] == "pTSS"].groupby(['gTSS_gene'])['gTSS_anot'].transform('count')
    second_filter=""
    filter_name=""
    if isCombined:
        second_filter="numSamples"
        filter_name="max_sample"
        annotated_df.loc[annotated_df["gTSS_anot"] == "pTSS", filter_name]=annotated_df.loc[annotated_df["gTSS_anot"] == "pTSS"].groupby(['gTSS_gene'])[second_filter].transform(max)
        annotated_df.loc[annotated_df["gTSS"] == "gTSS", "gTSS_anot_final"]=annotated_df.loc[annotated_df["gTSS"] == "gTSS"].apply(lambda x: "pTSS" if x["pTSS_count"] == 1 else ("pTSS" if ((x[second_filter] == x[filter_name]) and (math.isnan(x["pTSS_count"]) == False)) else "sTSS"),axis=1)
        annotated_df = annotated_df[["chrom","start","end","tss_type","cov", "strand","gTSS_anot_final" ,"gTSS_gene"  ,"tssGene", second_filter, "geneFrom", "max_cov", 'clusters']]
    else: # if is tss obtained from a single sample, instead of break even by numSamples reporting the TSS, use the minimum TSS start site
        second_filter="start"
        filter_name="min_start"
        annotated_df.loc[annotated_df["gTSS_anot"] == "pTSS", "min_start"]=annotated_df.loc[annotated_df["gTSS_anot"] == "pTSS"].groupby(['gTSS_gene'])[second_filter].transform(min)
        annotated_df.loc[annotated_df["gTSS"] == "gTSS", "gTSS_anot_final"]=annotated_df.loc[annotated_df["gTSS"] == "gTSS"].apply(lambda x: "pTSS" if x["pTSS_count"] == 1 else ("pTSS" if ((x[second_filter] == x[filter_name]) and (math.isnan(x["pTSS_count"]) == False)) else "sTSS"),axis=1)
        annotated_df = annotated_df[["chrom","start","end","tss_type","cov", "strand","gTSS_anot_final" ,"gTSS_gene"  ,"tssGene",  "geneFrom", "max_cov"]]
    if sample == "":
        return annotated_df.sort_values(["chrom", "strand", "start"], ascending=[["NC_005823.1", "NC_005824.1"], ["+", "-"], True])
    else: 
        sample_df=pd.read_csv(sample, sep="\t")
        sample_df.update(sample_df.groupby(['chrom', 'strand','clusters']).ffill())
        present_df = sample_df.groupby(['chrom', 'strand','clusters']).tail(1)
        del present_df['start']
        combined_df = annotated_df.merge(present_df,  on=["chrom", "strand", "clusters"], how="left").sort_values(["chrom", "strand", "start"], ascending=[["NC_005823.1", "NC_005824.1"], ["+", "-"], True]).reset_index(drop=True)
        return combined_df

all_genes=get_genes(gff) 
all_tss=get_tss(all_genes, filter_tss)
tss_anot=annotate_tss(all_genes, all_tss)
anot_df = annotate_gTSS(tss_anot)

anot_df.to_csv(output, sep="\t", index=False)
            
            
            
    
    
                        
                        
                
                        

            
        

        
        