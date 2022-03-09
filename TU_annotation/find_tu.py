import os
import sys
import pybedtools
import re
import subprocess
import argparse
from statistics import mode, multimode


# parser=argparse.ArgumentParser(prog="findTU",usage='%(prog)s [TSS|TTS] [options]',description='Arguments to identify Leptospira Transcript TSS:')
# parser.add_argument('feature',metavar='',type=str,help='specify if TSS or TTS need to be annotated')
# parser.add_argument( "-p",'--prefix',metavar='',type=str,help='prefix of the output file (only prefix), default is input file name', default="")
# parser.add_argument( "-o",'--output',metavar='',type=str,help='output file dir')
# parser.add_argument( "-b",'--bam',metavar='',type=str,help='read alignment in bam format')
# parser.add_argument( "-a",'--gff',metavar='',type=str,help='gff annotation file for the reference genome (.gff', default="")
# parser.add_argument( "-F",'--forcebed',metavar='',type=str,help='new bed file will always generate from bam in the same dir', default=True)
# parser.add_argument( "-t",'--threshold',metavar='',type=str,help='definition of read positions to classify in a cluster, default is within 10bps', default=10)
# parser.add_argument( "-cf",'--genecovfilter',metavar='',type=str,help='proportion of total read coverage needed in a read cluster to keep the TSS, default is 0.3', default=0.3)
# parser.add_argument( "-rf",'--genereadfilter',metavar='',type=str,help='number of reads needed in a read cluster to keep the TSS, default is 2 reads', default=2)
# args = parser.parse_args()

# feature=args.feature
# if !(feature == "TSS" | feature == "TTS"):
#     print("feature to analyze must be either TSS or TTS\nCurrent feature is not recognized\nprogram exiting...\n")
#     exit()
# window=int(args.threshold)
# gene_aware_cov_filter=int(args.genecovfilter)
# gene_aware_count_filter=int(args.genereadfilter)
# forceBed=bool(args.forcebed)

# bam=args.bam
# gff=args.gff
# sample=os.path.basename(bam).strip(".bam")
# output_prefix=os.path.join(args.output,sample) if args.prefix == "" else os.path.join(args.output, args.prefix)
# gene_output=output_prefix + "."+ feature +  ".gene_aware.tab"
# unaware_output=output_prefix + "."+ feature +".unaware.tab" 
# merged_output=output_prefix + "."+ feature  +  ".merged.tab"

window=10
gene_aware_cov_filter=0.3
gene_aware_count_filter=2
feature="TTS"

bam="/scratch/rx32940/minION/polyA_directRNA/map/genome/bam/Q29_Copenhageni_Basecalled_May_22_2020_Direct-RNA_rna_filtered.linear.bam"
gff="/scratch/rx32940/minION/polyA_cDNA/map/genome/reference/GCF_000007685.1_ASM768v1_genomic.gff"
gene_output="/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/input/Q29_Copenhageni_Basecalled_May_22_2020_Direct-RNA_rna_filtered.linear.gene_aware.tts。tab"
unaware_output="/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/input/Q29_Copenhageni_Basecalled_May_22_2020_Direct-RNA_rna_filtered.linear.unaware.tab"
merged_output="/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/input/Q29_Copenhageni_Basecalled_May_22_2020_Direct-RNA_rna_filtered.linear.merged.tab"

# create bed from bam
forceBed=True # always create a new bed file
if(os.path.isfile(bam.strip(".bam") + ".bed") and forceBed != True):
    print("bed file exist, will use the current bed")
    bed_file=bam.strip(".bam") + ".bed"
else:
    bam_file=pybedtools.BedTool(bam)
    bam_file.bam_to_bed().saveas(bam.strip(".bam") + ".bed")
    bed_file=bam.strip(".bam") + ".bed"

# read in bed file
with open(bed_file) as bedfile:
    bedlist=bedfile.readlines()

# get chroms in the bed file
bash_str="cat " + bed_file
p1 = subprocess.Popen(bash_str.split(), stdout=subprocess.PIPE, text=True)
p2 = subprocess.Popen(['awk', '{print $1}'], stdin=p1.stdout ,stdout=subprocess.PIPE, text=True)
p3 = subprocess.Popen("uniq", stdin=p2.stdout ,stdout=subprocess.PIPE, text=True)
output, error = p3.communicate()

chroms = output.split("\n")
chroms_list = [chrom for chrom in chroms if chrom != ""]

# put all reads in a dict 
# with structure {chromI: {+:[(start, end)], lags:[(start, end)]}, chromII: {+:[(start, end)], lags:[(start, end)]}}
all_reads={}
all_genes={} #dict to keep all gene annotations (positions)
pos_cov={}
for chrom in chroms_list:
    all_reads[chrom]={}
    all_genes[chrom]={}
    pos_cov[chrom]={}
    if feature == "TSS": # if detecting TTS, read end for + strand will be clustered, read start for - strand will be clustered
        start=1
        end=2
    else: # if detecting tts, reads mapped to lead strand will be treated as reads mapped in the lagging strand in the tss analysis
        start=2
        end=1
    # return start and end of each read in a tuple
    all_reads[chrom]['+'] = [(int(read.split('\t')[start]), int(read.split('\t')[end])) for read in bedlist if ('\t+' in read) == True & (chrom in read) == True]
    # lagging strand has end site as start
    all_reads[chrom]['-'] = [(int(read.split('\t')[end]), int(read.split('\t')[start])) for read in bedlist if ('\t-' in read) == True & (chrom in read) == True]
    # initiate other dict
    all_genes[chrom]['+']=[]
    all_genes[chrom]['-']=[]
    pos_cov[chrom]['+']={}
    pos_cov[chrom]['-']={}

# reference genome annotation
genome_file=gff.strip(".gff") + ".genome.tab"
with open(gff) as g, open(genome_file, "w") as genome:
    for line in g:
        line_list = line.split("\t")
        if line[0] != "#" and line_list[2] == "region": # skip headers, write chrom region to genome file (for later calculating genomce cov)
            a=genome.write("%s\t%d\n"%(line_list[0], int(line_list[4])))
            # print("%s\t%d"%(line_list[0], int(line_list[4])))
        elif line[0] != "#" and line_list[2] != "CDS": 
            gene_info=line_list[8].split(";")
            gene_id=gene_info["ID=" in gene_info].split("-")[1]
            if line_list[6] == "+":
                strand = "+"
            else:
                strand = "-"
            all_genes[line_list[0]][strand].append((int(line_list[3]), int(line_list[4]), gene_id))

            

def get_cluster_pos(read_list, strand, feature):
    if feature == "TTS":
        cur_strand = "-" if strand == "+" else "+"
    else:
        cur_strand = strand
    if cur_strand == "+":
        start_list=[x[0] for x in read_list]
        cluster_start = (min(multimode(start_list)) if start_list.count(mode(start_list)) > 1 else min(start_list)) + 1
        cluster_end = max([x[1] for x in read_list])
    else:
        start_list=[x[0] for x in read_list]
        cluster_start = max(multimode(start_list)) if start_list.count(mode(start_list)) > 1 else max(start_list)
        cluster_end = min([x[1] for x in read_list])
    
    return cluster_start, cluster_end

chrom="NC_005823.1"
strand="+"
gene=all_genes['NC_005823.1']['+'][1]
def get_overlap_reads(gene,strand,cur_reads_list,feature):
    if feature == "TTS": # if tts, then gene end will be checked, if read start > gene end, then iteration will stop (for lead strand mapped read)
        cur_strand="-" if strand == "+" else "+"
    else:                    
        cur_strand=strand
    gene_start=gene[0] if cur_strand == "+" else gene[1] # antisense initiate transcription from gene's 3' end
    # get reads overlap with gene start site
    cur_read_index = 0
    overlap_read_index =[]
    # use to check if we can stop iterating the reads
    read_start=0 if cur_strand == "+" else 1  # always check the 5' end of the read 
    read_end =1  if cur_strand == "+" else 0
    for read in cur_reads_list:
        # if sorted read start site is already larger than gene_start site, then stop the iterating mapped reads
        if read[read_start] > gene_start: 
            break
        if read[read_start] <= gene_start and read[read_end] >= gene_start:
            overlap_read_index.append(cur_read_index)
        # count current index of read
        cur_read_index += 1
    return overlap_read_index

"""
function to combine tss produced with gene-aware and unaware approach into the same dict
"""
def combine_tss(tss_c1, tss_c2):
    i=0
    j=0
    combined_tss=[]
    
    while i < len(tss_c1) and j < len(tss_c2):
        if tss_c1[i][0] < tss_c2[j][0]:
            combined_tss.append(tss_c1[i])
            i += 1
        elif tss_c2[j][0] < tss_c1[i][0]:
            combined_tss.append(tss_c2[j])
            j += 1
        else:
            combined_tss.append(tss_c1[i])
            combined_tss.append(tss_c2[j])
            i += 1
            j += 1
    
    while i < len(tss_c1):
        combined_tss.append(tss_c1[i])
        i +=1
    while j < len(tss_c2):
        combined_tss.append(tss_c2[j])
        j +=1
    return combined_tss

"""
gene_aware TSS identification

"""

with open(gene_output, "w") as go:
    gene_clusters={}
    for chrom in all_genes.keys():
        gene_clusters[chrom]={}
        for strand in all_genes[chrom].keys():
            gene_clusters[chrom][strand]=[]
            all_reads[chrom][strand].sort() # sort all reads based on the start site (for lag strand, is end site)
            all_genes[chrom][strand].sort() # sort all genes based on the start site (for lag strand, is end site)
            for gene in all_genes[chrom][strand]:
                gene_id=gene[2]
                overlap_reads_index=get_overlap_reads(gene, strand, all_reads[chrom][strand], feature)
                total_overlapped_cov=len(overlap_reads_index)
                # print(overlap_reads_index)
                if total_overlapped_cov > 1:
                    pre_read=all_reads[chrom][strand][overlap_reads_index[0]]
                    pre_start=pre_read[0]
                    cur_cluster=[pre_read]
                    for i in overlap_reads_index[1:]: # loop through all reads overlapped with gene start/end site
                        read = all_reads[chrom][strand][i]
                        if read[0] - pre_start <= window:   
                            cur_cluster.append(read)
                            pre_start = read[0]
                        else:
                            cluster_start, cluster_end = get_cluster_pos(cur_cluster, strand, feature)
                            # print(len(cur_cluster))
                            if (len(cur_cluster) >= (total_overlapped_cov * gene_aware_cov_filter)) and len(cur_cluster) >= gene_aware_count_filter:
                                gene_clusters[chrom][strand].append((cluster_start, cluster_end, gene_id, len(cur_cluster)))
                                a=go.write("%s\t%d\t%d\t%s\t%d\t%s\n"%(chrom, cluster_start, cluster_end, gene_id, len(cur_cluster), strand))
                            pre_start = read[0]
                            cur_cluster = [read]
                    # the last tss cluster when rest of the reads didn't get to go into the else condition to write their cluster into the df
                    cluster_start, cluster_end = get_cluster_pos(cur_cluster, strand, feature)
                    if (len(cur_cluster) >= (total_overlapped_cov * gene_aware_cov_filter)) and len(cur_cluster) >= gene_aware_count_filter:
                        # print("%s\t%d\t%d\t%s\t%d\t%s\n"%(chrom, cluster_start, cluster_end, gene_id, len(cur_cluster), strand))
                        gene_clusters[chrom][strand].append((cluster_start, cluster_end, gene_id, len(cur_cluster)))
                        b=go.write("%s\t%d\t%d\t%s\t%d\t%s\n"%(chrom, cluster_start, cluster_end, gene_id, len(cur_cluster), strand))


""""
gene_unaware TSS identification:
    all chroms
    + and lag strand
    input:
        dict all_reads{chromI: {+:[(start, end)], lags:[(start, end)]}, chromII: {+:[(start, end)], lags:[(start, end)]}}
    output:
        dict all_unaware_clusters{chromI: {+:[(start, end, numReads)], lags:[(start, end, numReads)]}, chromII: {+:[(start, end, numReads)], lags:[(start, end, numReads)]}}
"""

# calculate coverage at each base for each strand (use bedtools genome_coverage)
lead_bam_cov = pybedtools.bedtool.BedTool(bed_file).genome_coverage(g=genome_file, d=True, strand="+")
lag_bam_cov = pybedtools.bedtool.BedTool(bed_file).genome_coverage(d=True,g=genome_file, strand="-")

# save cov at each position into a dict
for pos in open(lead_bam_cov.fn): # lead
    line=pos.split("\t")
    chrom=line[0]
    pos_cov[chrom]["+"][int(line[1])] = int(line[2])

for pos in open(lag_bam_cov.fn): # lag
    line=pos.split("\t")
    chrom=line[0]
    pos_cov[chrom]["-"][int(line[1])] = int(line[2])

with open(unaware_output, "w") as uo:
    # gene unaware clustering
    all_unaware_clusters={} # this is the dict to save all the read can be clustered together without annotation
    for chrom in all_reads.keys(): # loop through chroms
        # chrom="NC_005823.1"
        all_unaware_clusters[chrom]={}
        for strand in all_reads[chrom].keys(): # loop through strands # HOW TO HANDLE LAGGING STRAND???????
            # strand="-"
            all_unaware_clusters[chrom][strand]=[]
            all_reads[chrom][strand].sort() # sort all reads based on the start site (for lag strand, is end site)
            cur_cluster_start_index = 0 # current cluster starts from this read index, reset for every new cluster
            cur_read_index = 0 # will not reset for all reads in the same strand and chrom
            first_read=all_reads[chrom][strand][0] # save the first read
            pre_start=first_read[0] # first read start
            cur_cluster = [first_read] # keep all end to a list
            for read in all_reads[chrom][strand][1:]: # start looping from second read
                start=read[0]
                end=read[1]
                cur_read_index += 1
                if (start - pre_start <= window): # if the start sites of the two reads are less than 10
                    cur_cluster.append(read) # save the end sites to a list
                    pre_start = start # now the start of the current read become the pre-start for next iteration
                else: # if the start sites of the two reads are further than 10, save all the reads from previous loops into a cluster
                    # the cluster starts with the first read index recorded for this cluster, and end with the maximum end position for all the reads in this cluster, 
                    # also with size of the cluster
                    cluster_start, cluster_end = get_cluster_pos(cur_cluster, strand,feature)
                    if len(cur_cluster) > pos_cov[chrom][strand][cluster_start]:
                        a=uo.write("%s\t%d\t%d\t.\t%d\t%s\t%s\n"%(chrom,cluster_start, cluster_end,len(cur_cluster), strand, "filter=" + str(pos_cov[chrom][strand][cluster_start]))) 
                        all_unaware_clusters[chrom][strand].append((cluster_start, cluster_end,"unaware",len(cur_cluster)))
                    cur_cluster=[read] # end list reinitiated for the new cluster.
                    pre_start=start
                    cur_cluster_start_index = cur_read_index # the new clusters starts, from the current read 
            cluster_start, cluster_end = get_cluster_pos(cur_cluster, strand,feature)
            if len(cur_cluster) > pos_cov[chrom][strand][cluster_start]:
                all_unaware_clusters[chrom][strand].append((cluster_start, cluster_end,"unaware",len(cur_cluster)))
                b=uo.write("%s\t%d\t%d\t.\t%d\t%s\t%s\n"%(chrom,cluster_start, cluster_end,len(cur_cluster), strand, "filter=" + str(pos_cov[chrom][strand][cluster_start]))) 




"""
combine gene-aware and unaware TSS
    - TSS within 10bps from two approaches are merged in one

"""

with open(merged_output, "w") as mo:
    merged_clusters={}
    for chrom in gene_clusters.keys():
        merged_clusters[chrom]={}
        for strand in gene_clusters[chrom].keys():
            combined_list = combine_tss(gene_clusters[chrom][strand], all_unaware_clusters[chrom][strand])
            merged_tss=[]
            pre_tss = combined_list[0]
            for tss in combined_list[1:]:
                if tss[0] - pre_tss[0] <= window:
                    merged_start=pre_tss[0]
                    merged_end=max([pre_tss[1], tss[1]])
                    merged_gene=tss[2] if pre_tss[2] == "unaware" else pre_tss[2]
                    merged_cov=(tss[3]+pre_tss[3])/2
                    pre_tss=(merged_start, merged_end, merged_gene, merged_cov)
                else:
                    merged_tss.append(pre_tss)
                    a=mo.write("%s\t%d\t%d\t%s\t%d\t%s\n"%(chrom,pre_tss[0], pre_tss[1],pre_tss[2],pre_tss[3], strand))
                    pre_tss=tss 
            merged_clusters[chrom][strand]=merged_tss
            

""""
TU identification:

"""
for chrom in all_unaware_clusters.keys():
    for strand in all_unaware_clusters[chrom].keys():
        first_cluster=all_unaware_clusters[chrom][strand]
        pre_start = first_cluster[0]
        pre_end = first_cluster[1]
        for cls in all_unaware_clusters[chrom][strand][1:]:
            cur_start = cls[0]
            if cur_start <= pre_end:
                overlap_start = cur_start
                overlap_end = pre_end


                                





    
    




    



