import os
import sys
import pybedtools
import re
import subprocess
import argparse
from statistics import mode, multimode
import find_feature


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

window=int(args.threshold)
gene_aware_cov_filter=float(args.genecovfilter)
gene_aware_count_filter=int(args.genereadfilter)
forceBed=bool(args.forcebed)

bam=args.bam
gff=args.gff
sample=os.path.basename(bam).strip(".bam")
output=args.output
output_prefix=os.path.join(output,sample) if args.prefix == "" else os.path.join(output, args.prefix)
tu_output=output_prefix  +  ".tu.tab"

tss_obj=find_feature.FindFeature(bam=bam, gff=gff, output=output, feature="TSS")
tss = tss_obj.find_features()
tts_obj=find_feature.FindFeature(bam=bam, gff=gff, output=output, feature="TTS")
tts = tts_obj.find_features()


