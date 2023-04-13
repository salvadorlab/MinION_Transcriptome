#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=cutadapt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=100:00:00
#SBATCH --mem=30G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL

PYCHOPPER="/scratch/rx32940/minION/polyA_cDNA/pychopper/output"
OUT="/scratch/rx32940/minION/polyA_cDNA/cutadapt/out"

ml cutadapt/3.4-GCCcore-8.3.0-Python-3.7.4

for file in $PYCHOPPER/*/*.final.fastq;
do

sample=$(basename $file '.final.fastq')

echo $sample

cutadapt \
-g "TTTCTGTTGGTGCTGATATTGCTGGG" \
-e 1 \
-j 4 \
-o $OUT/$sample.cutadapt.fastq \
$file

done