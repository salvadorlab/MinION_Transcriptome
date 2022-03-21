import os
import sys
import pybedtools
import re
import subprocess
import argparse
from statistics import mode, multimode
import find_feature
from collections import Counter

class FindTU:
    def __init__(self, forceBed = True, bam="", gff="", output="", prefix=""):
        if bam == "" or gff == "":
            self.__parse_args()
        else:
            forceBed=bool(forceBed)
            self.bam=bam
            self.gff=gff
            self.sample=os.path.basename(self.bam).strip(".bam")
            self.output=output
            output_prefix=os.path.join(self.output,self.sample) if prefix == "" else os.path.join(self.output, prefix)
            self.tu_output=output_prefix  +  ".operon.tab"
            tss_obj=find_feature.FindFeature(bam=self.bam, gff=self.gff,feature="TSS")
            self.all_reads = tss_obj.all_reads
            self.all_genes = tss_obj.all_genes
            
        
    def __parse_args(self):
        parser=argparse.ArgumentParser(prog="findTU",usage='%(prog)s [options]',description='Arguments to identify transcriptional operons:')
        parser.add_argument( "-p",'--prefix',metavar='',type=str,help='prefix of the output file (only prefix), default is input file name', default="")
        parser.add_argument( "-o",'--output',metavar='',type=str,help='output file dir')
        parser.add_argument( "-b",'--bam',metavar='',type=str,help='read alignment in bam format')
        parser.add_argument( "-a",'--gff',metavar='',type=str,help='gff annotation file for the reference genome (.gff)', default="")
        parser.add_argument( "-F",'--forcebed',metavar='',type=str,help='new bed file will always generate from bam in the same dir', default=True)
        args = parser.parse_args()

        forceBed=bool(args.forcebed)
        self.bam=args.bam
        self.gff=args.gff
        self.sample=os.path.basename(self.bam).strip(".bam")
        self.output=args.output
        output_prefix=os.path.join(self.output,self.sample) if args.prefix == "" else os.path.join(self.output, args.prefix)
        self.tu_output=output_prefix  +  ".operon.tab"
        tss_obj=find_feature.FindFeature(bam=self.bam, gff=self.gff,feature="TSS")
        self.all_reads = tss_obj.all_reads
        self.all_genes = tss_obj.all_genes
    # bam="/scratch/rx32940/minION/polyA_directRNA/map/genome/bam/Q29_Copenhageni_Basecalled_May_22_2020_Direct-RNA_rna_filtered.linear.bam"
    # gff="/scratch/rx32940/minION/polyA_cDNA/map/genome/reference/GCF_000007685.1_ASM768v1_genomic.gff"
    # output="/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/input/"


    def __find_spread_reads(self):
        all_operons={}
        for chrom in self.all_reads.keys():
            all_operons[chrom]={}
            for strand in self.all_reads[chrom].keys():
                all_operons[chrom][strand] = []
                if strand == "-": # lag strand reads has second position as 5'
                    read_start=1
                    read_end=0
                else:
                    read_start=0
                    read_end=1
                for read in self.all_reads[chrom][strand]: # loop through each read
                    cur_operon_start = read[read_start]
                    cur_operon_end = read[read_end]
                    cur_operon_genes = ""
                    gene_in_operon = 0
                    i=0
                    while i < len(self.all_genes[chrom][strand]): # rotating through genes (read is stationary)
                        cur_gene = self.all_genes[chrom][strand][i]
                        if cur_gene[1] < cur_operon_start: # genes end is before read start
                            i += 1
                        elif cur_gene[0]> cur_operon_end : # if gene start after read ends, then no more future gene will overlap with this read if gene list sorted
                            cur_operon = (cur_operon_start, cur_operon_end, gene_in_operon, cur_operon_genes[1:])
                            i = len(self.all_genes[chrom][strand]) # stop looping through genes
                        else: # gene over lap with current read
                            gene_in_operon +=1 
                            cur_operon_genes = cur_operon_genes + ":"+ cur_gene[2] 
                            i += 1
                    all_operons[chrom][strand].append(cur_operon)
        return all_operons

    def find_operon(self):
        all_operons = self.__find_spread_reads()
        with open(self.tu_output, "w") as ow:
            ow.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%("chrom" ,"strand","start", "end","numGenes","numOfReads", "genesIncluded")) 
            operon_cov_dict={}
            for chrom in all_operons.keys():
                operon_cov_dict[chrom]={}
                for strand in all_operons[chrom].keys():
                    operon_cov_dict[chrom][strand]=[]
                    # get unique operons (with different genes included) and get number of reads supporting each operon
                    genes_in_operons=[x[3] for x in all_operons[chrom][strand]]
                    operons_cov = Counter(genes_in_operons)
                    for op in operons_cov.keys():
                        op_start= min([x[0] for x in all_operons[chrom][strand] if x[3] == op])
                        op_end = max([x[1] for x in all_operons[chrom][strand] if x[3] == op])
                        op_spread=[x[2] for x in all_operons[chrom][strand] if x[3] == op][0]
                        op_cov = operons_cov[op]
                        operon_cov_dict[chrom][strand].append((op_start, op_end,op_spread ,op_cov, op))
                        ow.write("%s\t%s\t%d\t%d\t%d\t%d\t%s\n"%(chrom, strand, op_start, op_end,op_spread ,op_cov, op)) 
            return operon_cov_dict
                    
def main():

    bam=""
    gff=""
    tu_obj = FindTU(bam, gff)
    tu_obj.find_operon()
    

if __name__ == "__main__":
    main()         

