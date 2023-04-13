#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=fastqc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=100:00:00
#SBATCH --mem=30G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL

ml FastQC/0.11.9-Java-11
PYCHOPPER="/scratch/rx32940/minION/polyA_cDNA/pychopper/output"
OUT="/scratch/rx32940/minION/polyA_cDNA/pychopper/fastq_QC"

for fastq in $PYCHOPPER/*/*.final.fastq;
do
sample=$(basename $fastq '.final.fastq')

mkdir -p $OUT/$sample

fastqc $fastq -o $OUT/$sample

done

ml multiqc/1.11-GCCcore-8.3.0-Python-3.8.2