import os

gff_files = [file for file in os.listdir("/scratch/rx32940/minION/polyA_cDNA/map/genome/reference/" )if ".gff" in file]
noncode_file = "/scratch/rx32940/minION/unannotated/noncoding_reads_genomic_positions.tab"
output = "/scratch/rx32940/minION/unannotated/noncoding_reads_genomic_positions_anot.tab"

gene_dict={}
for gff_file in gff_files:
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

with open(noncode_file) as noncode, open(output, "w+") as anot:
    header=noncode.readline()
    for region in noncode:
        ln_ls=region.strip('\n').split('\t')
        chrom=ln_ls[2]
        region_start = int(ln_ls[8])
        region_end = int(ln_ls[9])
        for i in range(len(gene_dict[chrom])):
            gene = gene_dict[chrom][i]
            gene_start = int(gene[1])
            gene_end = int(gene[2])
            if gene_end < region_start:
                continue
            else: #gene_end >= region_start
                before_gene = gene_dict[chrom][i-1]
                after_gene = gene
                break
        line = region.strip('\n') + "\t" + '\t'.join([before_gene[0], before_gene[3],before_gene[6],after_gene[0], after_gene[3],after_gene[6]]) + "\n"
        a = anot.write(line)
        
                    