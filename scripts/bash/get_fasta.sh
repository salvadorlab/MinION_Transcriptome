#!/bin/bash
#SBATCH --partition=batch_30d
#SBATCH --job-name=get_30bps_around_3end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=100:00:00
#SBATCH --mem=50G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL

ml BEDTools/2.30.0-GCC-8.3.0
ml seqkit/0.16.1

REF="/scratch/rx32940/minION/polyA_cDNA/map/genome/reference"
WORK="/scratch/rx32940/minION/unannotated"

cat $REF/*fna > $WORK/reference_genomic_combined.fna

cat $WORK/noncoding_reads_genomic_positions_anot.tab | awk -F'\t' '{print $3"\t"$9"\t"$10"\t"$11","$14"\t"$4}' > $WORK/noncoding_regions.bed

sort $WORK/noncoding_regions.bed > $WORK/noncoding_regions.sorted.bed

split -l 1000 $WORK/noncoding_regions.sorted.bed $WORK/split_beds/noncoding_regions

for bed in $WORK/split_beds/*;
do

sample=$(basename $bed)

bedtools getfasta -fi $WORK/reference_genomic_combined.fna -bed $bed | seqkit seq -g -M 200000 | seqkit rmdup -s > $WORK/split_fna/$sample.noDup.fna 

done

