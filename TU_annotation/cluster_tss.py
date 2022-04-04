import os
import argparse
import sys
import subprocess

####
# input: all samples's tss annotation combined using combine_samples.r script
#       - only TSS from the same chrom and same strand should be used
#       - TSS identified within 10 bps from each other is considered in the same cluster
# output: add TSS cluster id to the last column of each TSS row
######

parser=argparse.ArgumentParser(prog="tss_cluster",usage='%(prog)s[options]',description='Arguments to identify Leptospira Transcript TSS:')
parser.add_argument('--input', "-i",type=str,help='TSS annotation from all samples combined, ex. cat *.merged.tab | grep chromI_ID | grep "+" > chromIlead.combined.tab')
parser.add_argument('--output', "-o",type=str,help='output file')
parser.add_argument('--threshold', "-t",type=str,help='definition of reads in each cluster, default is within 10bps', default=10)
parser.add_argument('--colindex', "-c",type=str,help='index of TSS positition column, default is 1, start from 0', default=1)
parser.add_argument('--chromIndex', "-ci",type=str,help='column index of the chromosome, default is 0', default=0)
parser.add_argument('--strand', "-s",type=str,help='column index of the chromosome, default is 3', default=3)
args = parser.parse_args()

input=args.input
output=args.output
cluster_range=int(args.threshold)
col_index=int(args.colindex)
chromIndex=int(args.chromIndex)
strandIndex=int(args.strand)

# input="/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/direct_output/tss/combined/filtered/combined_dRNA.filtered20.TSS.anot.tab"
# output="/home/rx32940/github/MinION_Transcriptome/TU_annotation/combined.clustered.TSS.tab"
# cluster_range=10
# col_index=1
# chromIndex=0
# strandIndex=5

# get chroms in the bed file
bash_str="cat " + input
p1 = subprocess.Popen(bash_str.split(), stdout=subprocess.PIPE, text=True)
p2 = subprocess.Popen(['awk', '{print $1}'], stdin=p1.stdout ,stdout=subprocess.PIPE)
p3 = subprocess.Popen("uniq", stdin=p2.stdout ,stdout=subprocess.PIPE, text=True)
bashout, error = p3.communicate()
chroms = [x for x in list(set(bashout.strip().split("\n"))) if x != "chrom"]

tss_dict={}
for chrom in chroms:
    tss_dict[chrom]={}
    tss_dict[chrom]["+"]=[]
    tss_dict[chrom]["-"]=[]

with open(input, 'r') as tss_file, open(output, "w") as new_file:
    headers = tss_file.readline()
    new_file.write(headers.strip('\n') + "\t" + "clusters\n")
    # separate tss based on chrom and strand
    for line in tss_file:
        line_list=line.split("\t")
        tss_dict[line_list[chromIndex]][line_list[strandIndex]].append((int(line_list[col_index]),line.strip("\n")))
    # cluster tss on each chrom and strand separately
    for chrom in tss_dict.keys():
        for strand in tss_dict[chrom].keys():
            start=0
            cluster = []
            i=0
            cluster_dict={}
            tss_dict[chrom][strand].sort()
            for line in tss_dict[chrom][strand]:
                tss = int(line[1].split('\t')[col_index])
                if tss-start <= cluster_range:
                    cluster.append(tss)
                    start=int(tss)
                    a=new_file.write(line[1] + "\t" + str(i) + "\n")
                else:
                    cluster.append(tss)
                    start=int(tss)
                    cluster_dict[i] = cluster
                    cluster=[]
                    i += 1
                    a=new_file.write(line[1] + "\t" + str(i) + "\n")



            
            
            
        
    