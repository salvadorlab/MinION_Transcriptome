import os
import sys
import pybedtools
import re
import subprocess

window=10

bam="/mnt/d/Dropbox/downloads/Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_PolyATail_rna_filtered.linear.bam"
bed="/mnt/d/Dropbox/downloads/Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_PolyATail_rna_filtered.linear.bed"
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
all_clusters={} # this is the dict to save all the read can be clustered together
for chrom in chroms_list:
    all_reads[chrom]={}
    all_clusters[chrom]={}
    # return start and end of each read in a tuple
    all_reads[chrom]["lead"] = [(int(read.split('\t')[1]), int(read.split('\t')[2])) for read in bedlist if ("\t+" in read) == True & (chrom in read) == True]
    # lagging strand has end site as start
    all_reads[chrom]["lag"] = [(int(read.split('\t')[2]), int(read.split('\t')[1])) for read in bedlist if ("\t-" in read) == True & (chrom in read) == True]
    all_clusters[chrom]["lead"]=[]
    all_clusters[chrom]["lag"]=[]

for chrom in all_reads.keys(): # loop through chroms
    print(chrom)
    # chrom="NC_005823.1"
    for strand in all_reads[chrom].keys(): # loop through strands # HOW TO HANDLE LAGGING STRAND???????
        print(strand)
        # strand="lead"
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
                all_clusters[chrom][strand].append((all_reads[chrom][strand][cur_cluster_start_index][0], max(end_list), len(end_list))) 
                end_list=[end] # end list reinitiated for the new cluster.
                pre_start=start
                cur_read_index += 1
                cur_cluster_start_index = cur_read_index # the new clusters starts, from the current read 

# Next step: cluster overlapping clusters into ORFs, how to filter clusters
                                





    
    




    



