#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=minIONQC
#SBATCH --ntasks=1                      
#SBATCH --cpus-per-task=3            
#SBATCH --time=100:00:00
#SBATCH --mem=50G
#SBATCH --output=%x.%j.out       
#SBATCH --error=%x.%j.out        
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL 

ml R/4.0.0-foss-2019b

QC="/scratch/rx32940/minION/QC"
cDNA="/scratch/rx32940/minION/Guppy/demultiplex/rebasecall"


for file in $cDNA/b*;
do
Rscript $QC/MinIONQC.R -i $file -o $cDNA/../QC -p 3 -f pdf
done
