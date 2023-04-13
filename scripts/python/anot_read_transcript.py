import os
import sys

input_file=sys.argv[1]
output_file=sys.argv[2]
gff_file=sys.argv[3]

# polya_path="/scratch/rx32940/minION/polyA_cDNA/tailfindr/relative_map_position"
# map_path="/scratch/rx32940/minION/polyA_cDNA/map/genome/bam"
# input_file=os.path.join(map_path, "Patoc_Basecalled_Aug_16_2019_Direct-cDNA_NoPolyATail.bed")
# output_file=os.path.join(polya_path, "Patoc_Basecalled_Aug_16_2019_Direct-cDNA_NoPolyATail.anot.bed")
# gff_file="/scratch/rx32940/minION/polyA_cDNA/map/genome/reference/GCF_000017685.1_ASM1768v1_genomic.gff"

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
                
     
with open(input_file) as infile, open(output_file, "w+") as outfile:
    for line in infile:
        line=line.strip("\n")
        ln_ls=line.split("\t")
        map_chrom=ln_ls[0]
        map_start=int(ln_ls[1])
        map_end=int(ln_ls[2])
        find=0
        for gene in gene_dict[map_chrom]:
            if find == 1:
                break
            else:
                if map_start > int(gene[2]): # not overlap with gene CDS
                    # print("not this gene")
                    continue
                else:
                    if map_end >= int(gene[2]) and map_start <= int(gene[1]): # read complete covers gene CDS
                        outline="%s\t%s\t%s\n"%(line, "\t".join(gene), "complete")
                    elif map_end >= int(gene[2]) and map_start < int(gene[2]) and map_start > int(gene[1]) and gene[3] == "+": # read only has partial 3' overlap with gene CDS on lead strand
                        outline="%s\t%s\t%s\n"%(line, "\t".join(gene), "partial(3')")
                    elif map_end >= int(gene[2]) and map_start < int(gene[2]) and map_start > int(gene[1]) and gene[3] == "-": # read only has partial 3' overlap with gene CDS on lag strand
                        outline="%s\t%s\t%s\n"%(line, "\t".join(gene), "partial(5')")
                    elif map_start <= int(gene[1])  and map_end < int(gene[2]) and map_end > int(gene[1]) and gene[3] == "+": # read only has partial 5' overlap with gene CDS on lead strand
                        outline="%s\t%s\t%s\n"%(line, "\t".join(gene), "partial(5')")
                    elif map_start<= int(gene[1])  and map_end < int(gene[2]) and map_end > int(gene[1]) and gene[3] == "-": # read only has partial 5' overlap with gene CDS on lag strand
                        outline="%s\t%s\t%s\n"%(line, "\t".join(gene), "partial(3')")
                    elif map_start > int(gene[1]) and map_end < int(gene[2]): # read completely embeded in gene CDS
                        outline="%s\t%s\t%s\n"%(line, "\t".join(gene), "internal")
                    else:
                        outline="%s\t%s\t%s\n"%(line, "\t".join(gene), "noncoding")
                    a=outfile.write(outline)
                    find=1 
                
                        
                    
                    
                    
    