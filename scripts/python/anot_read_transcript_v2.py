import os
import sys
import tqdm

input_file=sys.argv[1]
output_file=sys.argv[2]
gff_file=sys.argv[3]

# polya_path="/scratch/rx32940/minION/polyA_cDNA/tailfindr/relative_map_position"
# map_path="/scratch/rx32940/minION/polyA_cDNA/map/genome/bed"
# input_file=os.path.join(map_path, "Patoc_Basecalled_Aug_16_2019_Direct-cDNA-NoPolyATail.bed")
# output_file="/scratch/rx32940/minION/polyA_cDNA/test.bed"
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
                
with tqdm.tqdm(total=os.path.getsize(input_file)) as pbar:     
    with open(input_file) as infile, open(output_file, "w+") as outfile:
        for line in infile: # loop through every read
            b=pbar.update(len(line))
            line=line.strip("\n")
            ln_ls=line.split("\t")
            map_chrom=ln_ls[0] # chrom current read mapped to
            map_start=int(ln_ls[1])
            map_end=int(ln_ls[2])
            find=0 #find read map position in genomic annotation
            for i in range(len(gene_dict[map_chrom])): # loop through all gene CDS anot in cur chrom
                gene = gene_dict[map_chrom][i] # cur gene
                if find == 1: # if find, to next read
                    break
                else:
                    if map_start > int(gene[2]): # not overlap with gene CDS, move to next gene
                        # print("not this gene")
                        continue
                    else:
                        if map_end >= int(gene[2]) and map_start <= int(gene[1]): # read complete covers gene CDS
                            gene_ls=[gene]
                            j = i+1
                            gene_count = 1
                            while j < len(gene_dict[map_chrom]) and map_end > int(gene_dict[map_chrom][j][1]): # check if the read is covering CDS of multiple genes (operon), loop through following CDS
                                if (map_end - int(gene_dict[map_chrom][j][1])) >= (int(gene_dict[map_chrom][j][2]) - int(gene_dict[map_chrom][j][1]))*0.9: # if mapped length is at least 90% of the current gene
                                    gene_ls=gene_ls + [gene_dict[map_chrom][j]]
                                    gene_count +=1
                                j +=1
                            if gene_count > 1:# if operon
                                gene_str = []
                                for z in range(len(gene)): # merge all covered genes anot
                                    gene_str =gene_str + ["|".join([item[z] for item in gene_ls])]
                                outline="%s\t%s\t%s\n"%(line, "\t".join(gene_str), "operon|" + str(gene_count))
                            else: # else anot as complet covering read
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

                    
                            
                        
                        
                        
        