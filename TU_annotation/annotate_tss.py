import os
import pandas as pd
import argparse


parser=argparse.ArgumentParser(prog="tss_anot",usage='%(prog)s[options]',description='anot each tss as gTSS (pTSS/sTSS); iTSS; asTSS; oTSS:')
parser.add_argument('--input', "-i",type=str,help='TSS files')
parser.add_argument('--gff', "-g",type=str,help='reference gene annotation')
parser.add_argument('--output', "-o",type=str,help='output file')
args = parser.parse_args()

filter_tss=args.input
gff=args.gff
output=args.output

# filter_tss="/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/direct_output/tss/chromIlag.filtered20.TSS.tab"
# gff="/scratch/rx32940/minION/polyA_cDNA/map/genome/reference/GCF_000007685.1_ASM768v1_genomic.gff"
# output="/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/direct_output/tss/chromIlead.annotate20.TSS.tab"


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
            all_tss[line_list[0]][line_list[5]].append((int(line_list[1]), int(line_list[2]), line_list[3],int(line_list[4]),int(line_list[6])))
    return all_tss

def annotate_tss(all_genes, all_tss):
    all_genes=all_genes
    all_tss=all_tss
    tss_anot={}
    for chrom in all_tss.keys():
        tss_anot[chrom]={}
        for strand in all_tss[chrom].keys():
            tss_anot[chrom][strand]=[]
            if strand == "+":
                for item in all_tss[chrom][strand]:
                    tss=int(item[0])
                    for gene in all_genes[chrom][strand]:
                        if tss > gene[1]:
                            continue
                        elif gene[0] - tss <=300:
                            tss_anot[chrom][strand].append((*item, "gTSS", gene[2]))
                            break
                        elif tss >= gene[0]:
                            tss_anot[chrom][strand].append((*item, "iTSS", gene[2]))
                            break
                        else:
                            tss_anot[chrom][strand].append((*item, "oTSS", "."))
                            break
            else:
                for item in all_tss[chrom][strand]:
                    tss=int(item[0])
                    for gene in all_genes[chrom][strand]:
                        if (tss < gene[0]):
                            tss_anot[chrom][strand].append((*item, "oTSS", "."))
                            break
                        elif (tss - gene[1] <= 100) and (tss - gene[1] >=0):
                            tss_anot[chrom][strand].append((*item, "asTSS", gene[2]))
                            break
                        elif tss>= gene[0] and tss <= gene[1]:
                            tss_anot[chrom][strand].append((*item, "asTSS", gene[2]))
                            break
                        else:
                            continue

        return tss_anot


def annotate_gTSS(tss_anot):
    annotated_df=pd.DataFrame(columns=["start","end", "geneFrom", "cov", "numSamples", "tss_type", "tssGene", "chrom", "strand"])
    for chrom in tss_anot.keys():
        for strand in tss_anot[chrom].keys():
            cur_df = pd.DataFrame(tss_anot[chrom][strand], columns=["start","end", "geneFrom", "cov", "numSamples", "tss_type", "tssGene"])
            cur_df["chrom"]=chrom
            cur_df["strand"]=strand
            annotated_df = pd.concat([annotated_df, cur_df], ignore_index=True)
    annotated_df.loc[annotated_df["tss_type"] == "gTSS", "max_cov"]=annotated_df.loc[annotated_df["tss_type"] == "gTSS"].groupby(['tssGene'])['cov'].transform(max)
    annotated_df.loc[annotated_df["tss_type"] == "gTSS", "gTSS_anot"]=annotated_df.loc[annotated_df["tss_type"] == "gTSS"].apply(lambda x: "pTSS" if x["cov"] == x["max_cov"] else "sTSS",axis=1)
    return annotated_df


all_genes=get_genes(gff) 
all_tss=get_tss(all_genes, filter_tss)
tss_anot=annotate_tss(all_genes, all_tss)
anot_df = annotate_gTSS(tss_anot)[["chrom","start","end","tss_type","cov", "strand","gTSS_anot" ,"tssGene", "numSamples", "geneFrom", "max_cov"]]

anot_df.to_csv(output, sep="\t", index=False)
            
            
            
    
    
                        
                        
                
                        

            
        

        
        