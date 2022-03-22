import os
import argparse
from statistics import mode, multimode


####
# input: all samples's tss annotation combined using combine_samples.r script and clustered using cluster tss (choose the cluster threshold you examined)
#       - only TSS from the same chrom and same strand should be used
#       - TSS identified within 10 bps from each other is considered in the same cluster
# output: Will extracted the min TSS from each cluster (if lagging strand, will get max), and return number of samples in each cluster
######

parser=argparse.ArgumentParser(prog="tss_filter",usage='%(prog)s[options]',description='Arguments to get one TSS per cluster:')
parser.add_argument('--input', "-i",type=str,help='combined and clustered TSS files (only take output of combine_samples.r followed by cluster_tss.py)')
parser.add_argument('--output', "-o",type=str,help='output file')
parser.add_argument('--strand', "-s",type=str,help='current strand tss was detected from')
args = parser.parse_args()

tss_file=args.input
filter_file=args.output
strand=args.strand

# tss_file="/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/direct_output/tss/chromIlead.clustered20.TSS.tab"
# filter_file="/home/rx32940/github/MinION_Transcriptome/TU_annotation/chromIlead.filtered20.TSS.tab"
# strand="+"

with open(tss_file) as tf, open(filter_file, "w") as ff:
    header=tf.readline().strip("\n")
    a=ff.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%("chrom", "start", "end", "gene", "strand", "numSamples","clusters"))
    cur_line=tf.readline().strip("\n").split("\t")
    cur_cluster = int(cur_line[11])
    cur_tss=int(cur_line[1])
    cur_tss_cluster=[cur_tss]
    cur_ends=[int(x) for x in cur_line[4:11] if x != "NA"]
    cur_gene=cur_line[2]
    num_in_cluster = 0
    for line in tf:
        cur_line=line.strip("\n").split("\t")
        if int(cur_line[11]) == cur_cluster:
            cur_cluster = int(cur_line[11])
            num_in_cluster += 1
            cur_tss=int(cur_line[1])
            cur_tss_cluster.append(int(cur_tss))
            cur_ends= cur_ends+ [int(x) for x in cur_line[4:11] if x != "NA"]
            # print(cur_ends)
            cur_gene= str(cur_gene + ":" + cur_line[2]) if cur_line[2] != cur_gene else cur_gene
        else:
            if strand == "+":
                min_tss = min(multimode(cur_tss_cluster))  if cur_tss_cluster.count(mode(cur_tss_cluster)) > 1 else min(cur_tss_cluster)
                max_end = max(multimode(cur_ends)) if cur_ends.count(mode(cur_ends)) > 1 else max(cur_ends) 
            else:
                min_tss = max(multimode(cur_tss_cluster))  if cur_tss_cluster.count(mode(cur_tss_cluster)) > 1 else max(cur_tss_cluster)
                max_end = min(multimode(cur_ends)) if cur_ends.count(mode(cur_ends)) > 1 else min(cur_ends) 
            a=ff.write("%s\t%d\t%d\t%s\t%s\t%d\t%d\n"%(cur_line[0], min_tss,max_end,cur_gene, strand, len(cur_ends),cur_cluster))
            # print("%s\t%d\t%d\t%s\t%s\t%d\t%d\n"%(cur_line[0], min_tss,max_end,cur_gene, strand, len(cur_ends),cur_cluster))
            cur_cluster = int(cur_line[11])
            cur_tss=int(cur_line[1])
            cur_tss_cluster = [cur_tss]
            cur_ends=[int(x) for x in cur_line[4:11] if x != "NA"]
            cur_gene=cur_line[2]
            num_in_cluster = 0