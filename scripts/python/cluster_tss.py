import os
import argparse
import sys

####
# input: all samples's tss annotation combined using tu_annotate.Rmd script
#       - only TSS from the same chrom and same strand should be used
#       - TSS identified within 10 bps from each other is considered in the same cluster
# output: add TSS cluster id to the last column of each TSS row
######

parser=argparse.ArgumentParser(prog="tss_anot",usage='%(prog)s[options]',description='Arguments to identify Leptospira Transcript TSS:')
parser.add_argument('--input', "-i",type=str,help='TSS annotation from all samples (for a specific chrom or strand only from RMD script tu_annotate.rmd)')
parser.add_argument('--output', "-o",type=str,help='output dir')
parser.add_argument('--threshold', "-t",type=str,help='definition of reads in each cluster, default is within 10bps', default=10)
parser.add_argument('--colindex', "-c",type=str,help='index of TSS positition column, default is 2, start from 0', default=2)
args = parser.parse_args()

input=args.input
output=args.output
cluster_range=int(args.threshold)
col_index=int(args.colindex)
with open(input, 'r') as tss_file, open(output, "w+") as new_file:
    headers = tss_file.readline()
    new_file.write(headers)
    start=0
    cluster = []
    i=0
    cluster_dict={}
    for line in tss_file:
        tss = int(line.split(',')[col_index])
        print(start)
        if tss-start <= cluster_range:
            cluster.append(tss)
            start=int(tss)
            new_file.write(line.strip('\n') + "," + str(i) + "\n")
        else:
            cluster.append(tss)
            start=int(tss)
            cluster_dict[i] = cluster
            cluster=[]
            i += 1
            new_file.write(line.strip('\n') + "," + str(i) + "\n")

            
            
            
        
    