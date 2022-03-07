import pybedtools.bedtool as pybed
import os
import pandas as pd
from tqdm import tqdm
import math
import sys
import argparse


parser=argparse.ArgumentParser(prog="tss_anot",usage='%(prog)s [TSS|TTS] [options]',description='Arguments to identify Leptospira Transcript TSS:')
parser.add_argument('feature',metavar='',type=str,help='specify if TSS or TTS need to be annotated')
parser.add_argument( "-g",'--gene',metavar='',type=str,help='gene annotation file')
parser.add_argument( "-o",'--output',metavar='',type=str,help='output file dir')
parser.add_argument( "-b",'--bam',metavar='',type=str,help='read alignment in bam format')
parser.add_argument( "-c",'--cov',metavar='',type=str,help='genome coverage file at each position of the bam alignment')
parser.add_argument( "-t",'--threshold',metavar='',type=str,help='definition of read positions to classify in a cluster, default is within 10bps', default=10)
parser.add_argument( "-cf",'--genecovfilter',metavar='',type=str,help='proportion of total read coverage needed in a read cluster to keep the TSS, default is 0.3', default=0.3)
parser.add_argument( "-rf",'--genereadfilter',metavar='',type=str,help='number of reads needed in a read cluster to keep the TSS, default is 2 reads', default=2)
args = parser.parse_args()

workdir="/scratch/rx32940/minION/polyA_directRNA/TU_Annotation"

FEATURE=args.feature

if FEATURE is "TSS":
    lead_start="start"
    lead_end="end"
    lag_start="end"
    lag_end="start"
else:
    lead_start="end"
    lead_end="start"
    lag_start="start"
    lag_end="end"
    

output=os.path.join(args.output,"temp")
if not os.path.exists(output):
    os.makedirs(output)

gene_file=args.gene # gene annotation file

bam=args.bam # bam
sample_name=bam.split(".")[0].split("/")[-1]
bed_file=os.path.join(output, str(sample_name) + ".bed") # bed

bam = pybed.BedTool(bam)
bed=bam.bam_to_bed().to_dataframe()
bed.to_csv(bed_file,sep="\t", index = False)

# filter for gene-aware approach
genecovfilter=int(args.genecovfilter)
genereadfilter=int(args.genereadfilter)

    
bed_filtered = bed#.loc[bed["score"] >= 30] # filter read quality

read_map_lead=bed_filtered.loc[bed["strand"] == "+"]
read_map_lag = bed_filtered.loc[bed["strand"] == "-"]

gene_anot=pd.read_csv(gene_file, sep="\t")

lead_strand_genes = gene_anot.loc[gene_anot["strand"] == "+"]
lag_strand_genes = gene_anot.loc[gene_anot["strand"] == "-"]

threshold = int(args.threshold) # the position of reads with start codon less than threshold are put into the same cluster 

chroms = gene_anot["genomic_accession"].unique() # chromosomes in the genome

# function to get reads overlap at either "start" or "end" (start in lagging strand) of a gene
def get_overlap_reads(read_mapped_cur_chrom, cur_chrom_gene ,gene, position):
    
    gene_start = int(cur_chrom_gene.loc[cur_chrom_gene["locus_tag"]==gene,position]) # start codon of the current gene at 5' end
    reads_before_start = read_mapped_cur_chrom.loc[read_mapped_cur_chrom["start"] <= gene_start] # reads with start sites mapped before gene start codon
    reads_overlap_start = reads_before_start.loc[reads_before_start["end"] >= gene_start] # from reads start before start codon, get reads overlapped with the start codon (ends after start codon)
    reads_overlap_start.reset_index(drop=True, inplace=True)
    return reads_overlap_start

def get_tss_pos(reads_overlap, cluster_start, cluster_end,position):
    tss_start= min(reads_overlap.loc[cluster_start:cluster_end, "start"])
    tss_end = max(reads_overlap.loc[cluster_start:cluster_end, "end"]) # find the end sites (read with furtherest mapping) of the current cluster
    if position == "start":
        return tss_start+1, tss_end
    else:
        return tss_end, tss_start+1
            
def cluster_tss(reads_overlap,position, gene_name, chrom):
    cur_gene_tss_clusters = pd.DataFrame() # tss clusters at for current gene's start codon
    
    
    reads_overlap = reads_overlap.sort_values(by=position) # sort all reads overlapped with current gene's start codon based on their start sites
    reads_overlap.reset_index(drop=True, inplace=True)
    

    cluster_start = 0 # start position of the cluster  
    cluster_end = 0 # start position of the cluster  (not included)
    for index in range(1,len(reads_overlap.iloc[:,1])):
        current_read = reads_overlap.iloc[index]
        if (current_read[position] - reads_overlap.loc[index-1,position]) <= threshold: # if the start sites of two reads are close enough 
            cluster_end += 1 # assign to the same cluster
        else:# if the start site of reads are too distant from each other to be in the same cluster
            tss_pos = get_tss_pos(reads_overlap,cluster_start, cluster_end , position)
            cur_gene_tss_clusters = cur_gene_tss_clusters.append({"chrom" : cur_chrom,"start" : tss_pos[0], "end":tss_pos[1], "gene" : cur_gene,"num_reads":int(cluster_end - cluster_start  + 1)},ignore_index=True)
            cluster_start = cluster_end= index 
            index += 1

    # the last tss cluster when rest of the reads didn't get to go into the else condition to write their cluster into the df
    tss_pos = get_tss_pos(reads_overlap, cluster_start, cluster_end, position)
    cur_gene_tss_clusters = cur_gene_tss_clusters.append({"chrom" : chrom,"start" : tss_pos[0], "end":tss_pos[1], "gene" : gene_name,"num_reads":int(cluster_end - cluster_start  + 1)},ignore_index=True)

    return cur_gene_tss_clusters

lead_gene_tss_clusters = pd.DataFrame()

for cur_chrom in chroms:
    
    cur_chrom_gene = lead_strand_genes.loc[lead_strand_genes["genomic_accession"] == cur_chrom]
    cur_chrom_gene.reset_index(drop=True, inplace=True)
    
    read_mapped_cur_chrom = read_map_lead.loc[read_map_lead["chrom"]== cur_chrom]
    
    for gene in tqdm(range(len(cur_chrom_gene.iloc[:,1]))): # iterate through all genes
        cur_gene = cur_chrom_gene.loc[gene,"locus_tag"] # current gene name
        reads_overlap_start = get_overlap_reads(read_mapped_cur_chrom, cur_chrom_gene ,cur_gene, "start")
        cov_at_pos = len(reads_overlap_start.iloc[:,1]) # find the coverage of the current gene's start codon, (all the reads overlap at gene start codon)
        
        if cov_at_pos >= 1: # if the gene was mapped
            cur_gene_tss_clusters = cluster_tss(reads_overlap_start, "start", cur_gene, cur_chrom)
            lead_gene_tss_clusters  = lead_gene_tss_clusters.append(cur_gene_tss_clusters[(cur_gene_tss_clusters["num_reads"] > (cov_at_pos * genecovfilter)) | (cur_gene_tss_clusters["num_reads"] >= genereadfilter)])

# write lead strand tss annotation with start codon awareness into file
lead_gene_tss_clusters["strand"] = "+" # write lead strand gene's tss read clusters 

lag_gene_tss_clusters = pd.DataFrame()

for cur_chrom in chroms:
    
    cur_chrom_gene = lag_strand_genes.loc[lag_strand_genes["genomic_accession"] == cur_chrom]
    cur_chrom_gene.reset_index(drop=True, inplace=True)
    
    read_mapped_cur_chrom = read_map_lag.loc[read_map_lag["chrom"]== cur_chrom]
    
    for gene in tqdm(range(len(cur_chrom_gene.iloc[:,1]))): # iterate through all genes
        cur_gene = cur_chrom_gene.loc[gene,"locus_tag"] # current gene name
        reads_overlap_end = get_overlap_reads(read_mapped_cur_chrom, cur_chrom_gene ,cur_gene, "end")
        cov_at_pos = len(reads_overlap_end.iloc[:,1]) # find the coverage of the current gene's start codon, (all the reads overlap at gene start codon)
        
        if cov_at_pos >= 1: # if the gene was mapped
            cur_gene_tss_clusters = cluster_tss(reads_overlap_end, "end", cur_gene, cur_chrom)
            lag_gene_tss_clusters  = lag_gene_tss_clusters.append(cur_gene_tss_clusters[cur_gene_tss_clusters["num_reads"] > (cov_at_pos * genecovfilter) |  (cur_gene_tss_clusters["num_reads"] >= genereadfilter)])
            
lag_gene_tss_clusters["strand"]="-"

gene_aware_tss = pd.concat([lead_gene_tss_clusters, lag_gene_tss_clusters], ignore_index = True)

gene_aware_tss = gene_aware_tss.dropna() # remove nan in the dataframe
# gene_aware_tss= gene_aware_tss.loc[gene_aware_tss["num_reads"] > 2]

gene_aware_tss= gene_aware_tss.astype({'start':int, 'end':int, "num_reads":int})
gene_aware_tss.to_csv(os.path.join(output,str(sample_name)+ ".gene_aware_tss.tab"), sep="\t", index=False)

#Identification of TSS with no annotation

read_map_lead = bed_filtered.loc[bed_filtered["strand"] == "+"]
read_map_lag = bed_filtered.loc[bed_filtered["strand"] == "-"]

def cluster_tss_unaware(cur_chrom_reads,position, chrom):
    cur_tss_clusters = pd.DataFrame()
    cluster_start = 0 # start position of the cluster  
    cluster_end = 0 # start position of the cluster  (not included)
    for index in tqdm(range(1,len(cur_chrom_reads.iloc[:,1]))): # loop through reads
        current_read = cur_chrom_reads.iloc[index]
        if (current_read[position] - cur_chrom_reads.loc[index-1,position]) <= threshold: # if the start sites of two reads are close enough 
            cluster_end += 1 # assign to the same cluster
        else:# if the start site of reads are too distant from each other to be in the same cluster
            tss_pos = get_tss_pos(cur_chrom_reads,cluster_start, cluster_end , position)
            cur_tss_clusters = cur_tss_clusters.append({"chrom" : chrom,"start" : tss_pos[0], "end":tss_pos[1], "name":"gene_unaware","num_reads":int(cluster_end - cluster_start  + 1)},ignore_index=True)
            cluster_start = cluster_end = index 
            index += 1

    # the last tss cluster when rest of the reads didn't get to go into the else condition to write their cluster into the df
    tss_pos = get_tss_pos(cur_chrom_reads, cluster_start, cluster_end, position)
    cur_tss_clusters = cur_tss_clusters.append({"chrom" : cur_chrom,"start" : tss_pos[0], "end":tss_pos[1], "name":"gene_unaware","num_reads":int(cluster_end - cluster_start  + 1)},ignore_index=True)

    return cur_tss_clusters

tss_clusters = pd.DataFrame() # tss clusters at for all reads' start codon


for cur_chrom in chroms:
    
    # read mapped to one chromosome
    read_map_cur_chrom_lead = read_map_lead.loc[read_map_lead["chrom"] == cur_chrom]

    read_map_cur_chrom_lead = read_map_cur_chrom_lead.sort_values(by="start") # sort all reads overlapped with current gene's start codon based on their start sites
    read_map_cur_chrom_lead.reset_index(drop=True, inplace = True)
    cur_tss_clusters = cluster_tss_unaware(read_map_cur_chrom_lead, "start", cur_chrom)

    tss_clusters =tss_clusters.append(cur_tss_clusters ,ignore_index=True)

lead_tss_clusters = tss_clusters
lead_tss_clusters.reset_index(drop=True, inplace = True)
lead_tss_clusters["strand"] = "+"
lead_tss_clusters= lead_tss_clusters.astype({'start':int, 'end':int, "num_reads":int})
# lead_tss_clusters.to_csv(os.path.join(output, str(sample_name) + ".lead.tss.unaware.bed"), sep="\t", index = False)

tss_clusters = pd.DataFrame() # tss clusters at for all reads' start codon (lagging strand at 3' end)


for cur_chrom in chroms:
    
    # read mapped to one chromosome
    read_map_cur_chrom_lag = read_map_lag.loc[read_map_lag["chrom"] == cur_chrom]

    read_map_cur_chrom_lag = read_map_cur_chrom_lag.sort_values(by="end") # sort all reads overlapped with current gene's start codon based on their start sites
    read_map_cur_chrom_lag.reset_index(drop=True, inplace = True)
    cur_tss_clusters = cluster_tss_unaware(read_map_cur_chrom_lag, "end", cur_chrom)

    tss_clusters =tss_clusters.append(cur_tss_clusters ,ignore_index=True)
    

lag_tss_clusters = tss_clusters
lag_tss_clusters.reset_index(drop=True, inplace = True)
lag_tss_clusters["strand"] = "-"
lag_tss_clusters= lag_tss_clusters.astype({'start':int, 'end':int, "num_reads":int})
# lag_tss_clusters.to_csv(os.path.join(output, str(sample_name) +".lag.tss.unaware.bed"), sep="\t", index = False)

gene_unaware_tss= pd.concat([lead_tss_clusters,lag_tss_clusters], ignore_index = True)
gene_unaware_tss.to_csv(os.path.join(output, str(sample_name) +".tss.unaware.tab"), sep="\t", index = False)

# filter gene unaware tss
# separate reads mapped to leading and lagging strand in the input bam file 
lead_bed = bed_filtered.loc[bed["strand"] == "+"]
lag_bed = bed_filtered.loc[bed["strand"] == "-"]

# write reads mapped to each bam file in bed (separate file for each stand)
lead_bed.to_csv(os.path.join(output, str(sample_name) + ".bamtobed_lead.bed"), header=False,sep="\t", index=False)
lag_bed.to_csv(os.path.join(output, str(sample_name) +".bamtobed_lag.bed"), header=False,sep="\t", index = False)

# use bed tools to read in each strand's bed file
lead_bam_bed = pybed.BedTool(os.path.join(output, str(sample_name) + ".bamtobed_lead.bed"))
lag_bam_bed = pybed.BedTool(os.path.join(output, str(sample_name) +".bamtobed_lag.bed"))

# calculate coverage at each base for each strand (use bedtools genome_coverage)
lead_bam_cov = lead_bam_bed.genome_coverage(g=args.cov, d=True)
lag_bam_cov = lag_bam_bed.genome_coverage(d=True,g=args.cov)

lead_bam_cov_df = pd.read_table(lead_bam_cov.fn, names =[ "chrom", "start", "cov"] )
lag_bam_cov_df = pd.read_table(lag_bam_cov.fn, names=["chrom", "start", "cov"])

# write per-base cov to bed file
# lead_bam_cov_df.to_csv(os.path.join(output, str(sample_name) + ".bamtobed_lead_cov.bed"), header=False,sep="\t", index=False)
# lag_bam_cov_df.to_csv(os.path.join(output, str(sample_name) + ".bamtobed_lag_cov.bed"), header=False,sep="\t", index=False)

# get cov of the site of each detected tss
# lead
lead_unaware_tss_start_cov = lead_tss_clusters.merge(lead_bam_cov_df,how="left", on=["chrom", "start"]) 
# lead_unaware_tss_start_cov.to_csv(os.path.join(output,str(sample_name) + ".unaware_tss_filtered.tab")
lead_unaware_tss_filtered= lead_unaware_tss_start_cov.loc[lead_unaware_tss_start_cov["num_reads"] >= lead_unaware_tss_start_cov["cov"]] # filter Tss based on cov at the site
# lag
lag_unaware_tss_start_cov = lag_tss_clusters.merge(lag_bam_cov_df,how="left", on=["chrom", "start"]) 
lag_unaware_tss_filtered= lag_unaware_tss_start_cov.loc[lag_unaware_tss_start_cov["num_reads"] >= lag_unaware_tss_start_cov["cov"]]

# combine lead and lag strand, write filtered tss to file
lead_unaware_tss_filtered.reset_index(drop=True, inplace=True)
lag_unaware_tss_filtered.reset_index(drop=True, inplace=True)
gene_unaware_tss_filtered= pd.concat([lead_unaware_tss_filtered, lag_unaware_tss_filtered])
gene_unaware_tss_filtered.to_csv(os.path.join(output,str(sample_name) + ".unaware_tss_filtered.tab"), sep = "\t", index=False)


# combine TSS from both approach
combine_approaches = gene_unaware_tss_filtered[["chrom","start", "end", "strand"]].merge(gene_aware_tss, how="outer", on=["chrom","start","strand"])
combine_approaches["end"] = combine_approaches.apply(lambda row: row["end_y"] if math.isnan(row["end_x"]) else row["end_x"],axis=1)
combine_approaches_edited = combine_approaches.drop(["end_x", "end_y","num_reads"],axis=1)
combine_approaches_edited = combine_approaches_edited.astype({"start":int, "end":int})
combine_approaches_edited.to_csv(os.path.join(output, "../" + str(sample_name) + ".combined_tss.tab"), index=False, sep="\t")

