#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=minIONQC
#SBATCH --ntasks=1                      
#SBATCH --cpus-per-task=3            
#SBATCH --time=100:00:00
#SBATCH --mem=100G
#SBATCH --output=%x.%j.out       
#SBATCH --error=%x.%j.out        
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL 

ml R/4.0.0-foss-2019b
seq="/scratch/rx32940/minION/dRNA_Comparison/DATA/guppy/rebasecalled_fast5"
minIONQC="/scratch/rx32940/minION/dRNA_Comparison/DATA/QC"
minqcOUT="/scratch/rx32940/minION/dRNA_Comparison/DATA/QC/minionQC"

# for file in $seq/*;
# do
# Rscript $minIONQC/MinIONQC.R -i $file -o $minqcOUT -p 3 -f pdf
# done

Rscript $minIONQC/MinIONQC.R -i $seq/fast5_pass-LIC-NO-POLYA -o $minqcOUT -p 3 -f pdf
#########################################################################
#
# fastQC
#
##########################################################################
# ml FastQC/0.11.9-Java-11
# ml MultiQC/1.8-foss-2019b-Python-3.7.4

# fastqc="/scratch/rx32940/minION/dRNA_Comparison/DATA/QC/fastqc"


# all_fastq="$(ls "$seq"/*.fastq)"

# fastqc -t 3 -o $fastqc -f fastq $all_fastq


# multiqc $fastqc/*_fastqc.zip -o $fastqc/multiQC
