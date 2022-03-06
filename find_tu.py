import os
import sys
import pybedtools
import re
import subprocess

window=10

bam="/mnt/d/Dropbox/downloads/Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_PolyATail_rna_filtered.linear.bam"
bed="/mnt/d/Dropbox/downloads/Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_PolyATail_rna_filtered.linear.bed"
gff="/mnt/d/Dropbox/5.Rachel-projects/Transcriptomics_PolyATail/Map_statistics_10042021/reference/GCF_000007685.1_ASM768v1_genomic.gff"


# create bed from bam
forceBed=True # always create a new bed file
if bed != "":
    print('bed file provided, will use input instead of bam file')
    bed_file=bed
elif(os.path.isfile(bam.strip(".bam") + ".bed") & forceBed != True):
    print("bed file exist, will use the current bed")
    bed_file=bam.strip(".bam") + ".bed"
else:
    bam_file=pybedtools.BedTool(bam)
    bam_file.bam_to_bed().saveas(bam.strip(".bam") + ".bed")
    bed_file=bam.strip(".bam") + ".bed"

# read in bed file
with open(bed) as bedfile:
    bedlist=bedfile.readlines()

# get chroms in the bed file
bash_str="cat " + bed
p1 = subprocess.Popen(bash_str.split(), stdout=subprocess.PIPE, text=True)
p2 = subprocess.Popen(['awk', '{print $1}'], stdin=p1.stdout ,stdout=subprocess.PIPE, text=True)
p3 = subprocess.Popen("uniq", stdin=p2.stdout ,stdout=subprocess.PIPE, text=True)
output, error = p3.communicate()

chroms = output.split("\n")
chroms_list = [chrom for chrom in chroms if chrom != ""]

# put all reads in a dict 
# with structure {chromI: {lead:[(start, end)], lags:[(start, end)]}, chromII: {lead:[(start, end)], lags:[(start, end)]}}
all_reads={}
all_genes={} #dict to keep all gene annotations (positions)
all_gene_clusters={} # this is the dict to save all the read can be clustered together based on gene annotations
all_unaware_clusters={} # this is the dict to save all the read can be clustered together without annotation
for chrom in chroms_list:
    all_reads[chrom]={}
    all_genes[chrom]={}
    all_gene_clusters[chrom]={}
    all_unaware_clusters[chrom]={}
    # return start and end of each read in a tuple
    all_reads[chrom]["lead"] = [(int(read.split('\t')[1]), int(read.split('\t')[2])) for read in bedlist if ("\t+" in read) == True & (chrom in read) == True]
    # lagging strand has end site as start
    all_reads[chrom]["lag"] = [(int(read.split('\t')[2]), int(read.split('\t')[1])) for read in bedlist if ("\t-" in read) == True & (chrom in read) == True]
    # initiate other dict
    all_genes[chrom]["lead"]=[]
    all_genes[chrom]["lag"]=[]
    all_gene_clusters[chrom]["lead"]=[]
    all_gene_clusters[chrom]["lag"]=[]
    all_unaware_clusters[chrom]["lead"]=[]
    all_unaware_clusters[chrom]["lag"]=[]

# reference genome annotation
with open(gff) as g:
    for line in g:
        line_list = line.split("\t")
        if line[0] != "#" and line_list[2] != "region" and line_list[1] != "Protein Homology": # skip headers
            gene_info=line_list[8].split(";")
            gene_id=gene_info["ID=" in gene_info].split("-")[1]
            if line_list[6] == "+":
                strand = "lead"
            else:
                strand = "lag"
            all_genes[line_list[0]][strand].append((int(line_list[3]), int(line_list[4]), gene_id))

def get_cluster_pos(read_list, strand):
    if strand == "lead":
        cluster_start = min([x[0] for x in read_list])
        cluster_end = max([x[1] for x in read_list])
    else:
        cluster_start = max([x[0] for x in read_list])
        cluster_end = min([x[1] for x in read_list])
    
    return cluster_start, cluster_end


"""
gene_aware TSS identification

"""

for chrom in all_genes.keys():
    for strand in all_genes[chrom].keys():
        all_reads[chrom][strand].sort() # sort all reads based on the start site (for lag strand, is end site)
        all_genes[chrom][strand].sort() # sort all genes based on the start site (for lag strand, is end site)
        for gene in all_genes[chrom][strand]:
            gene_start=gene[0]
            gene_end=gene[1]
            gene_id=gene[2]
            # get reads overlap with gene start site
            is_overlap_started = 0
            cur_read_index = 0
            for read in all_reads[chrom][strand]:
                # if sorted read start site is already larger than gene_start site, then stop the iterating mapped reads
                if read[0] > gene_start: 
                    cur_cluster_reads = [all_reads[chrom][strand][i] for i in cur_cluster_index]
                    cluster_start, cluster_end = get_cluster_pos(cur_cluster_reads, strand)
                    all_gene_clusters[chrom][strand].append((cluster_start, cluster_end, gene[2],len(cur_cluster_reads)))
                    pre_start = read[0]
                    cur_cluster_index = [cur_read_index]
                    break
                if read[0] <= gene_start and read[1] >= gene_start:
                    if is_overlap_started == 0: # is this the first overlapped read
                        is_overlap_started = 1
                        pre_start = read[0]
                        cur_cluster_index = [cur_read_index]
                    elif read[0] - pre_start <= window:   
                        cur_cluster_index.append(cur_read_index)
                        pre_start = read[0]
                    else:
                        cur_cluster_reads = [all_reads[chrom][strand][i] for i in cur_cluster_index]
                        cluster_start, cluster_end = get_cluster_pos(cur_cluster_reads, strand)
                        all_gene_clusters[chrom][strand].append((cluster_start, cluster_end,gene[2] ,len(cur_cluster_reads)))
                        pre_start = read[0]
                        cur_cluster_index = [cur_read_index]
                # count current index of read
                cur_read_index += 1


""""
gene_unaware TSS identification:
    all chroms
    lead and lag strand
    filter: ???
    input:
        dict all_reads{chromI: {lead:[(start, end)], lags:[(start, end)]}, chromII: {lead:[(start, end)], lags:[(start, end)]}}
    output:
        dict all_unaware_clusters{chromI: {lead:[(start, end, numReads)], lags:[(start, end, numReads)]}, chromII: {lead:[(start, end, numReads)], lags:[(start, end, numReads)]}}
"""

for chrom in all_reads.keys(): # loop through chroms
    print(chrom)
    # chrom="NC_005823.1"
    for strand in all_reads[chrom].keys(): # loop through strands # HOW TO HANDLE LAGGING STRAND???????
        print(strand)
        # strand="lag"
        all_reads[chrom][strand].sort() # sort all reads based on the start site (for lag strand, is end site)
        cur_cluster_start_index = 0 # current cluster starts from this read index, reset for every new cluster
        cur_read_index = 0 # will not reset for all reads in the same strand and chrom
        first_read=all_reads[chrom][strand][0] # save the first read
        pre_start=first_read[0] # first read start
        end_list = [first_read[1]] # keep all end to a list
        for read in all_reads[chrom][strand][1:]: # start looping from second read
            start=read[0]
            end=read[1]
            if (start - pre_start <= window): # if the start sites of the two reads are less than 10
                end_list.append(end) # save the end sites to a list
                cur_read_index += 1
                pre_start = start # now the start of the current read become the pre-start for next iteration
            else: # if the start sites of the two reads are further than 10, save all the reads from previous loops into a cluster
                # the cluster starts with the first read index recorded for this cluster, and end with the maximum end position for all the reads in this cluster, 
                # also with size of the cluster
                if strand == "lead":
                    cluster_start = min([x[0] for x in all_reads[chrom][strand][cur_cluster_start_index: cur_read_index+1]])
                    cluster_end = max(end_list)
                else:
                    cluster_start = max([x[0] for x in all_reads[chrom][strand][cur_cluster_start_index: cur_read_index+1]])
                    cluster_end = min(end_list)
                all_unaware_clusters[chrom][strand].append((cluster_start, cluster_end, len(end_list))) 
                end_list=[end] # end list reinitiated for the new cluster.
                pre_start=start
                cur_read_index += 1
                cur_cluster_start_index = cur_read_index # the new clusters starts, from the current read 

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


                                





    
    




    



