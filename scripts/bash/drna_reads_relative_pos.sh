#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=anotQ29
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=15:00:00
#SBATCH --mem=50G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL

##########################################################################################
#
# In this script, we want to determine where in the genome did the all dRNA reads were mapped to
#
##########################################################################################


##########################################################################################
#
# Turn Bam to Bed
#
##########################################################################################

# BAM="/scratch/rx32940/minION/polyA_directRNA/map/genome/bam"
# ml BEDTools/2.30.0-GCC-8.3.0
# for bam in $BAM/*.sorted.bam;
# do
# sample=$(basename $bam .sorted.bam)

# bedtools bamtobed -i $bam -cigar > $BAM/$sample.bed

# done


##########################################################################################
#
# anot all dRNA reads (where did these reads mapped to)
# use python script anot_FL_polyA_reads.py
# arg 1: bed file from the previous step
# args 2: anot file, output file 
# args 3: GFF of the sample
#
##########################################################################################
BAM="/scratch/rx32940/minION/polyA_directRNA/map/genome/bam"
WORK="/scratch/rx32940/minION/polyA_directRNA/relative_map_pos"
REFP="/scratch/rx32940/minION/polyA_cDNA/map/genome/reference"
REF="$REFP/GCF_000007685.1_ASM768v1_genomic.gff"



# all dRNA reads
##########################################################################################
for bed in $BAM/Q29*RNA.bed;
do
sample=$(basename $bed '.bed')

echo $bed
time python $WORK/anot_FL_polyA_reads.py $bed $WORK/$sample.anot.bed $REF

done

##########################################################################################
#
# grep list of PASS (nanopolish) reads from bed file 
#
##########################################################################################
# BAM="/scratch/rx32940/minION/polyA_cDNA/map/genome/bam"
# WORK="/scratch/rx32940/minION/polyA_directRNA/relative_map_pos"
# NANO="/scratch/rx32940/minION/polyA_directRNA/nanopolish07252022/nanopolish_output"

# cat $NANO/Q29_Copenhageni_Direct-RNA_FAST5_Pass.nanopolish.pass.tsv | awk '{print $1}' > $WORK/Q29_pass_list.txt
# cat $NANO/Q29_Copenhageni_Direct-RNA_Repeat_FAST5_Pass.nanopolish.pass.tsv | awk '{print $1}' >> $WORK/Q29_pass_list.txt

# # # dRNA reads with a PASS status only (evaulated by nanopolish)
# ##########################################################################################
# BAM="/scratch/rx32940/minION/polyA_directRNA/map/genome/bam"
# for bed in $BAM/Q29*RNA.bed;
# do
# sample=$(basename $bed '.bed')
# grep -wFf $WORK/Q29_pass_list.txt $bed > $WORK/$sample.nanoPASS.bed
# done


##########################################################################################
#
# anot dRNA PASS reads (where did these reads mapped to)
# use python script anot_FL_polyA_reads.py
# arg 1: bed file from the previous step
# args 2: anot file, output file 
# args 3: GFF of the sample
#
##########################################################################################
# WORK="/scratch/rx32940/minION/polyA_directRNA/relative_map_pos"
# REFP="/scratch/rx32940/minION/polyA_cDNA/map/genome/reference"
# REF="$REFP/GCF_000007685.1_ASM768v1_genomic.gff"


# # # all dRNA reads
# ##########################################################################################
# for bed in $WORK/*.nanoPASS.bed;
# do
# sample=$(basename $bed '.bed')

# python $WORK/anot_FL_polyA_reads.py $bed $WORK/$sample.anot.bed $REF

# done

