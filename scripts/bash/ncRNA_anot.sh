#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=anot_blastX
#SBATCH --ntasks=1                    	
#SBATCH --cpus-per-task=12            
#SBATCH --time=100:00:00
#SBATCH --mem=10G
#SBATCH --output=%x.%j.out       
#SBATCH --error=%x.%j.out        
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL 


# TABLE="/Users/rx32940/Dropbox/5.Rachel-projects/Transcriptomics_PolyATail/manuscript/cdna/Lepto_transcriptome_ONT_cDNA/5.tables"

# # # ##### make BED file from the genomic regions identified
# cat $TABLE/noTransMap_genomic_regions_cDNA_region_anot.csv | \
# awk -F',' '{print $5"\t"$7"\t"$8"\t"$14":"$15"\t"$13"\t.\t"$11"\t"$12}' | \
# tail -n +2 > ../5.tables/noTransMap_genomic_regions.bed 

###################################################################################
#
# Move to Sapelo2 from this step
#
###################################################################################
# source activate tu_annotation

WKDIR="/scratch/rx32940/minION/ncRNA"
REF="/scratch/rx32940/minION/map/genome/reference"

# ml BEDTools/2.30.0-GCC-8.3.0
# ml GenomeTools/1.6.1-GCC-8.3.0-Python-2.7.16
# ml seqkit/0.16.1
# ml Infernal/1.1.3-foss-2019b

# # # #### extract sequences all regions from reference genomes only mapped during genome mapping
# for ref in $REF/*v1_genomic.fna;
# do
# echo $ref
# bedtools getfasta -bed $WKDIR/noTransMap_genomic_regions.bed -fi $ref >> $WKDIR/fasta/noTransMap_genomic_regions.fasta
# done

# # # # # ### only keep the unique sequences in the upstream fasta
# gt sequniq $WKDIR/fasta/noTransMap_genomic_regions.fasta > $WKDIR/fasta/noTransMap_genomic_regions.uniq.fasta

# rm $WKDIR/fasta/*.fasta.*

####### not use for now
# ####seqkit head -n 1000 $BED/all_oTSS_q29q36.uniq.m50.fasta > $BED/all_oTSS_q29q36.uniq.m50.fasta_1
# ####prodigal -p meta -i $BED/all_transcripts_q29q36.uniq.m50.fasta -o $BED/all_transcripts_q29q36.uniq.m50.orf.gff -f gff 
# ### bedtools getfasta -name -bed $BED/all_transcripts_q29q36.uniq.m50.orf.gff -fi $BED/all_transcripts_q29q36.uniq.m50.fasta -fo $BED/all_transcripts_q29q36.uniq.m50.orf.fasta
# ###seqkit head -n 1000 $BED/all_transcripts_q29q36.uniq.m50.orf.fasta > $BED/all_transcripts_q29q36.uniq.m50.orf.fasta_1
# ####seqkit range -r 1001:2000 $BED/all_transcripts_q29q36.uniq.m50.orf.fasta > $BED/all_transcripts_q29q36.uniq.m50.orf.fasta_2
# ####seqkit range -r 2001:3000 $BED/all_transcripts_q29q36.uniq.m50.orf.fasta > $BED/all_transcripts_q29q36.uniq.m50.orf.fasta_3

###################################################################################
#
# Rfam
#
###################################################################################

# BED="/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/direct_output/tss/sRNA"
# # z score calculation number of nucleotide (cat $WKDIR/fasta/noTransMap_genomic_regions.uniq.fasta | wc)*2/1000000
# # https://docs.rfam.org/en/latest/genome-annotation.html
# cmscan -Z 0.029984 --cut_ga --rfam --anytrunc --nohmmonly --cpu 8 --tblout $WKDIR/noTransMap_genomic_regions.rfam.tblout --fmt 2 --clanin $BED/Rfam/Rfam.clanin $BED/Rfam/Rfam.cm $WKDIR/fasta/noTransMap_genomic_regions.uniq.fasta > $WKDIR/Rfam/noTransMap_genomic_regions.Rfam.out

###################################################################################
#
# BLASTN
#
###################################################################################

# conda activate ncbi
# DBNAME="/scratch/rx32940/metagenomics/blastn/db/nt"

# time blastn -query $WKDIR/fasta/noTransMap_genomic_regions.uniq.fasta \
# -db $DBNAME -outfmt 0 -num_alignments 1 -num_descriptions 1 \
# -out $WKDIR/blast/noTransMap_genomic_regions.blast.txt -num_threads 12

##############################################################
# diamond blastx
##############################################################

source activate megan
DBNAME="/scratch/rx32940/metagenomics/megan/db"

# time diamond blastx -d $DBNAME/nr -q $WKDIR/fasta/noTransMap_genomic_regions.uniq.fasta -o $WKDIR/blast/noTransMap_genomic_regions.blastx -f 100 -p 12
time diamond view -a $WKDIR/blast/noTransMap_genomic_regions.blastx.daa > $WKDIR/blast/noTransMap_genomic_regions.blastx.txt
