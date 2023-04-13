#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=pychopper_rescue
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=100:00:00
#SBATCH --mem=30G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL

source activate pychopper

conda install -c epi2melabs -c conda-forge -c bioconda "epi2melabs::pychopper"

WORK="/scratch/rx32940/minION/polyA_cDNA"
OUT="/scratch/rx32940/minION/polyA_cDNA/pychopper/output"

for fastq in $WORK/fastq/*fastq;
do

sample=$(basename $fastq ".fastq")
mkdir -p $OUT/$sample 

echo $sample

pychopper \
-r $OUT/$sample/PCB109_report.pdf \
-A $OUT/$sample/aln_hits_PCB109.bed \
-S $OUT/$sample/PCB109_statistics.tsv \
-u $OUT/$sample/unclassified_PCB109.fastq \
-t 8 \
-w $OUT/$sample/rescue_PCB109.fastq \
$fastq \
$OUT/$sample/PCB109.FL.fastq

echo "############################################################################################################################################################"
done


##############################################################################
#
# rescue unclassified reads 
#
##############################################################################
echo "### DCS rescue ############################################################################################################################################################"

for dir in $OUT/*/;
do
sample=$(basename $dir)

echo $sample 

pychopper \
-S $OUT/$sample/DCS109_statistics.tsv \
-A $OUT/$sample/aln_hits_DCS109.bed \
-r $OUT/$sample/DCS109_report.pdf \
-x rescue \
-w $OUT/$sample/rescue_DCS109.fastq \
-u $OUT/$sample/unclassified_DCS109.fastq \
-t 8 \
$OUT/$sample/unclassified_PCB109.fastq \
$OUT/$sample/DCS109.FL.fastq

cat $OUT/$sample/DCS109.FL.fastq $OUT/$sample/PCB109.FL.fastq $OUT/$sample/rescued_PCB109.fastq $OUT/$sample/rescue_DCS109.fastq > $OUT/$sample/$sample.final.fastq

echo "############################################################################################################################################################"
done

