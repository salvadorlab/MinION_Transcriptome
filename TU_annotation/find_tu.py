import os
import sys
import pybedtools
import re
import subprocess
import argparse
from statistics import mode, multimode
import find_feature
from collections import Counter


parser=argparse.ArgumentParser(prog="findTU",usage='%(prog)s [options]',description='Arguments to identify Leptospira transcription units:')
parser.add_argument( "-p",'--prefix',metavar='',type=str,help='prefix of the output file (only prefix), default is input file name', default="")
parser.add_argument( "-o",'--output',metavar='',type=str,help='output file dir')
parser.add_argument( "-b",'--bam',metavar='',type=str,help='read alignment in bam format')
parser.add_argument( "-a",'--gff',metavar='',type=str,help='gff annotation file for the reference genome (.gff', default="")
parser.add_argument( "-F",'--forcebed',metavar='',type=str,help='new bed file will always generate from bam in the same dir', default=True)
parser.add_argument( "-t",'--threshold',metavar='',type=str,help='definition of read positions to classify in a cluster, default is within 10bps', default=10)
parser.add_argument( "-cf",'--genecovfilter',metavar='',type=str,help='proportion of total read coverage needed in a read cluster to keep the TSS, default is 0.3', default=0.3)
parser.add_argument( "-rf",'--genereadfilter',metavar='',type=str,help='number of reads needed in a read cluster to keep the TSS, default is 2 reads', default=2)
args = parser.parse_args()


forceBed=bool(args.forcebed)
bam=args.bam
gff=args.gff
sample=os.path.basename(bam).strip(".bam")
output=args.output
output_prefix=os.path.join(output,sample) if args.prefix == "" else os.path.join(output, args.prefix)
tu_output=output_prefix  +  ".operon.tab"

# bam="/scratch/rx32940/minION/polyA_directRNA/map/genome/bam/Q29_Copenhageni_Basecalled_May_22_2020_Direct-RNA_rna_filtered.linear.bam"
# gff="/scratch/rx32940/minION/polyA_cDNA/map/genome/reference/GCF_000007685.1_ASM768v1_genomic.gff"
# output="/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/input/"


def __find_spread_reads(all_reads, all_genes):
    all_operons={}
    for chrom in all_reads.keys():
        all_operons[chrom]={}
        for strand in all_reads[chrom].keys():
            all_operons[chrom][strand] = []
            if strand == "-": # lag strand reads has second position as 5'
                read_start=1
                read_end=0
            else:
                read_start=0
                read_end=1
            for read in all_reads[chrom][strand]: # loop through each read
                cur_operon_start = read[read_start]
                cur_operon_end = read[read_end]
                cur_operon_genes = ""
                gene_in_operon = 0
                i=0
                while i < len(all_genes[chrom][strand]): # rotating through genes (read is stationary)
                    cur_gene = all_genes[chrom][strand][i]
                    if cur_gene[1] < cur_operon_start: # genes end is before read start
                        i += 1
                    elif cur_gene[0]> cur_operon_end : # if gene start after read ends, then no more future gene will overlap with this read if gene list sorted
                        cur_operon = (cur_operon_start, cur_operon_end, gene_in_operon, cur_operon_genes[1:])
                        i = len(all_genes[chrom][strand]) # stop looping through genes
                    else: # gene over lap with current read
                        gene_in_operon +=1 
                        cur_operon_genes = cur_operon_genes + ":"+ cur_gene[2] 
                        i += 1
                all_operons[chrom][strand].append(cur_operon)
    return all_operons

def find_operon(all_reads, all_genes):
    all_operons = __find_spread_reads(all_reads, all_genes)
    with open(tu_output, "w") as ow:
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
    output=""
    tss_obj=find_feature.FindFeature(bam=bam, gff=gff, output=output, feature="TSS")
    # tss_dict = tss_obj.find_features()
    # tts_obj=find_feature.FindFeature(bam=bam, gff=gff, output=output, feature="TTS")
    # tts_dict = tts_obj.find_features() 

    all_reads = tss_obj.all_reads
    all_genes = tss_obj.all_genes
    find_operon(all_reads, all_genes)
    


if __name__ == "__main__":
    main()         

