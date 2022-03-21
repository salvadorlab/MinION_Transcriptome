import pysam
import os
import pandas as pd

# ml Pysam/0.16.0.1-GCC-8.3.0
# ml pandas/0.25.3-intel-2019b-Python-3.7.4
map_path="/scratch/rx32940/minION/polyA_cDNA/map"
map_to="transcriptome"

"""
create tabular index for reference annotation
 1) download gff3 file from NCBI
 2) sort gff based on contig positions: 
    2.1) download tool: GFF3sort (https://github.com/billzt/gff3sort.git)  
    2.2) need to load Bioperl to run: ml BioPerl/1.7.2-GCCcore-8.3.0
    2.3) sort gff: perl gff3sort/gff3sort.pl GCF_000007685.1_ASM768v1_genomic.gff | bgzip > GCF_000007685.1_ASM768v1_genomic_sorted.gff.gz
 3) now we can use tabix to index the gff file (if use cmd line tabix, also need use bgzip to re-compress the gff file)
"""

# tabix = pysam.tabix_index(os.path.join(map_path, map_to, "GCF_000007685.1_ASM768v1_genomic_sorted.gff.gz"), preset="gff")
# gff = pysam.TabixFile(os.path.join(map_path, map_to, "GCF_000007685.1_ASM768v1_genomic_sorted.gff.gz"))

# cur_sample
# for i in gff.fetch("NC_005823.1", parser=pysam.asGTF()):
#     print(str(i.attributes))

"""
Each read's mapping statistics to reference transcriptome

parse use pysam
"""
bam_path=os.path.join(map_path,map_to,"bam")

bam_files=[file for file in os.listdir(bam_path) if file.split(".")[-1] == "bam"]

all_sample=[]
trans_map_stat= pd.DataFrame()
for file in bam_files:
    
    sample_name=file.split(".")[0] # sample name from path
    abs_filepath=os.path.join(bam_path,file) # absolute path to file
    
    # read in the current bam file (rb, read bam)
    bamfile = pysam.AlignmentFile(abs_filepath, "rb")
    
    # number of reads mapped to each transcripts
    cur_sample_map_stat = pd.DataFrame(bamfile.get_index_statistics())
    cur_sample_map_stat["Sample"] = sample_name
    trans_map_stat = pd.concat([trans_map_stat,cur_sample_map_stat ])

    for transcript in bamfile.references: # loop through all transcripts

        transcript_name = transcript # current transcript
        transcript_len = bamfile.get_reference_length(transcript_name) # reference length of the transcript
        reads_map_transcript = bamfile.fetch(transcript_name) # fetch reads mapped to the current transcript from the current bam file

        for read in reads_map_transcript: # for every read in the reads mapped to the current transcript
            all_sample.append(
                {
                "Sample": sample_name,
                "Transcript Name": transcript_name,
                "Transcript Len": transcript_len,
                "ReadID": read.query_name,
                "Read Length": read.query_length, # include soft-clipped nucleotides but not hard-clipped
                "Aligned Length": read.reference_length, # aligned length of the read on the reference genome
                "Map Start": read.reference_start,
                "Map End":read.reference_end,
                "CIGAR string": read.cigarstring, 
                "is_secondary": read.is_secondary, 
                "is_supplementary": read.is_supplementary,
                "coverage fraction": round(read.reference_length/transcript_len,4), # get the coverage fraction of the current read to the reference transcript it mapped to
                "Mapping quality": read.mapping_quality
                }
                )


read_cov_frac = pd.DataFrame(all_sample) # turn list to dataframe

read_cov_frac.to_csv(os.path.join(map_path,map_to,"cDNA_read_map_cov_frac_to_transcript.csv"))

trans_map_stat.to_csv(os.path.join(map_path,map_to,"cDNA_transcript_read_map_stats.csv"))

