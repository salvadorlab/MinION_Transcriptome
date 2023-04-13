from functools import partial
import os
import sys
import pybedtools
import re
import subprocess
import argparse
from statistics import mode, multimode
from collections import Counter


def find_spread_reads():
    all_operons={}
    for chrom in all_reads.keys():
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
                    if (cur_gene[1] < cur_operon_start) or (cur_gene[0] < cur_operon_start): # genes starts or ends is before read start, no partial overlap
                        i += 1
                    elif (cur_gene[0] > cur_operon_end) or (cur_gene[1]> cur_operon_end): # if gene start after read ends, then no more future gene will overlap with this read if gene list sorted
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

      