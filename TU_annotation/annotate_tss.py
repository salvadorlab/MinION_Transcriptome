import os

filter_tss="/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/direct_output/tss/chromIlead.filtered20.TSS.tab"
gff="/scratch/rx32940/minION/polyA_cDNA/map/genome/reference/GCF_000007685.1_ASM768v1_genomic.gff"
output="/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/direct_output/tss/chromIlead.annotate20.TSS.tab"
col_index=1

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
all_genes=get_genes(gff) 

def get_tss(all_genes, filter_tss):
    all_genes=all_genes
    all_tss={}
    for chrom in all_genes.keys():
        all_tss[chrom]={}
        for strand in all_genes[chrom].keys():
            all_tss[chrom][strand]=[]
    with open(filter_tss) as ft:
        header= ft.readline()
        for line in ft:
            line_list=line.strip('\n').split('\t')
            all_tss[line_list[0]][line_list[4]].append((line_list[1], line_list[1], line_list[3],line_list[5],line_list[6]))
    return all_tss

all_tss=get_tss(filter_tss)

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
                        if tss > gene[0] and tss < gene[1]:
                            tss_anot[chrom][strand].append((*item, "iTSS"))
                        elif gene[0] - tss <=300:
                            tss_anot[chrom][strand].append((*item, "gTSS"))
                        else:
                            tss_anot[chrom][strand].append((*item, "oTSS"))
            else:
                for item in all_tss[chrom][strand]:
                    tss=int(item[0])
                    if tss - gene[1] <= 300:
                        tss_anot[chrom][strand].append((*item, "asTSS"))
                    else:
                        tss_anot[chrom][strand].append((*item, "oTSS"))
    return tss_anot

tss_anot=annotate_tss(all_genes, all_tss)

def annotate_gTSS(tss_anot):
    
    
                        
                        
                
                        

            
        

        
        