import os
import sys
import pybedtools
import re
import subprocess
import argparse
from statistics import mode, multimode
import find_feature


# parser=argparse.ArgumentParser(prog="findTU",usage='%(prog)s [options]',description='Arguments to identify Leptospira transcription units:')
# parser.add_argument( "-p",'--prefix',metavar='',type=str,help='prefix of the output file (only prefix), default is input file name', default="")
# parser.add_argument( "-o",'--output',metavar='',type=str,help='output file dir')
# parser.add_argument( "-b",'--bam',metavar='',type=str,help='read alignment in bam format')
# parser.add_argument( "-a",'--gff',metavar='',type=str,help='gff annotation file for the reference genome (.gff', default="")
# parser.add_argument( "-F",'--forcebed',metavar='',type=str,help='new bed file will always generate from bam in the same dir', default=True)
# parser.add_argument( "-t",'--threshold',metavar='',type=str,help='definition of read positions to classify in a cluster, default is within 10bps', default=10)
# parser.add_argument( "-cf",'--genecovfilter',metavar='',type=str,help='proportion of total read coverage needed in a read cluster to keep the TSS, default is 0.3', default=0.3)
# parser.add_argument( "-rf",'--genereadfilter',metavar='',type=str,help='number of reads needed in a read cluster to keep the TSS, default is 2 reads', default=2)
# args = parser.parse_args()

# window=int(args.threshold)
# gene_aware_cov_filter=float(args.genecovfilter)
# gene_aware_count_filter=int(args.genereadfilter)
# forceBed=bool(args.forcebed)

# bam=args.bam
# gff=args.gff
# sample=os.path.basename(bam).strip(".bam")
# output=args.output
# output_prefix=os.path.join(output,sample) if args.prefix == "" else os.path.join(output, args.prefix)
# tu_output=output_prefix  +  ".tu.tab"

bam="/scratch/rx32940/minION/polyA_directRNA/map/genome/bam/Q29_Copenhageni_Basecalled_May_22_2020_Direct-RNA_rna_filtered.linear.bam"
gff="/scratch/rx32940/minION/polyA_cDNA/map/genome/reference/GCF_000007685.1_ASM768v1_genomic.gff"
output="/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/input/"

tss_obj=find_feature.FindFeature(bam=bam, gff=gff, output=output, feature="TSS")
tss_dict = tss_obj.find_features()
tts_obj=find_feature.FindFeature(bam=bam, gff=gff, output=output, feature="TTS")
tts_dict = tts_obj.find_features()

chrom="NC_005823.1"
strand="+"
# tss[chrom][strand][0:5]
# tts[chrom][strand][0:5]


all_orfs={}
for chrom in tts_dict.keys():
    all_orfs[chrom]={}
    for strand in tts_dict[chrom].keys():
        all_orfs[chrom][strand]=[]
        print(strand)
        if strand == "+":
            for tts in tts_dict[chrom][strand][0:5]:
                # print(tts)
                tts_start=tts[0] # start of the tts, from 3' -> 5'
                tts_end=tts[1]
                i=0
                orf_list=[]
                tss_start= tss_dict[chrom][strand][i][0]
                tss_end= tss_dict[chrom][strand][i][1]
                while i < len(tss_dict[chrom][strand]) and tts_start >= tss_end:
                    # print("i=" + str(i))
                    print("tss start:%s\ntss_end:%s\ntts_start:%s\ntts_end:%s\n"%(tss_start,tss_end,tts_start,tts_end))
                    if tss_start <= tts_start and tts_end <= tss_end:
                        print((tss_start, tts_start))
                        gene_id = tts_dict[chrom][strand][2] if tss_dict[chrom][strand][2] == "unaware" else tss_dict[chrom][strand][2]
                        orf_list.append((tss_start, tts_start))
                    tss_start= tss_dict[chrom][strand][i][0]
                    tss_end= tss_dict[chrom][strand][i][1]
                    i +=1

            all_orfs[chrom][strand] = all_orfs[chrom][strand] + orf_list
                
            
        
            
            
        

