import os
import sys

bed=sys.argv[1]
gff_file=sys.argv[2]
output=sys.argv[3]

# bed="/scratch/rx32940/minION/homopolymer_cov_genome/bed/GCF_000007685.1_ASM768v1_homoA_5.bed"
# gff_file="/scratch/rx32940/minION/polyA_cDNA/map/genome/reference/GCF_000007685.1_ASM768v1_genomic.gff"
# output="/scratch/rx32940/minION/homopolymer_cov_genome/bed/GCF_000007685.1_ASM768v1_homoA_5.anot.bed"

gene_dict={}
with open(gff_file) as gff:
    for line in gff:
        line=line.strip("\n")
        line_ls=line.split("\t")
        if line[0] == "#":
            continue
        else:
            if line_ls[2] == "region":
                gene_dict[line_ls[0]]=[]
            elif line_ls[2] in ("CDS", "rRNA", "tRNA", "tmRNA", "riboswitch", "RNase_P_RNA"):
                anot_list=line_ls[8].split(";")
                locus=list(filter(lambda x: "locus_tag=" in x, anot_list))
                locus=locus[0].split("=")[1] if len(locus) > 0 else ""
                protein=list(filter(lambda x: "Name=" in x, anot_list))
                protein=protein[0].split("=")[1] if len(protein) > 0 else ""
                biotype=list(filter(lambda x: "gbkey=" in x, anot_list))
                biotype=biotype[0].split("=")[1] if len(biotype) > 0 else ""
                product=list(filter(lambda x: "product=" in x, anot_list))
                product=product[0].split("=")[1] if len(product) > 0 else ""
                gene_dict[line_ls[0]] =  gene_dict[line_ls[0]] + [(locus, line_ls[3], line_ls[4], line_ls[6], biotype,protein,product)]
            else:
                continue

with open(bed) as inputfile, open(output, "w+") as out:
    for line in inputfile:
        ln_ls = line.strip("\n").split("\t")
        chrom=ln_ls[0]
        homo_range = ln_ls[3].split(":")
        homo_start = homo_range[0]
        homo_end = homo_range[1]
        for gene in gene_dict[chrom]:
            if int(homo_start) > int(gene[2]):
                continue
            elif int(gene[1]) > int(homo_end):
                ln_ls.append("noncoding")
                a=out.write("\t".join(ln_ls) + "\n")
                break
            elif int(homo_start) <= int(gene[2]):
                gene_anot=gene[5] if gene[5] != "" else gene[0]
                score='.'
                strand=gene[3]
                # the start position of BED starts from 0-position, while GFF starts from 1, we will adjust GFF gene start into BED format as rest of our analysis is BED based.
                homo_region_start=ln_ls[1] if int(ln_ls[1]) > (int(gene[1]) -1) else str(int(gene[1])-1) # fix region 100 bps next to homo, if reset based on boundry of coding regions.
                homo_region_end=ln_ls[2] if int(ln_ls[2]) < int(gene[2]) else gene[2]
                anot = [ln_ls[0]] + [homo_region_start, homo_region_end ,gene_anot, score, strand] + ln_ls[3:] 
                a=out.write("\t".join(anot) + "\n")
                break
                
            
        
    
    