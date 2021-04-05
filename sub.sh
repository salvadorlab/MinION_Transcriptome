#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=batch
#SBATCH --ntasks=1                    	
#SBATCH --cpus-per-task=3            
#SBATCH --time=100:00:00
#SBATCH --mem=100G
#SBATCH --output=%x.%j.out       
#SBATCH --error=%x.%j.out        
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL 

ml R/4.0.0-foss-2019b

cd /scratch/rx32940/minION

for file in data/*;
do
Rscript MinIONQC.R -i $file -o /scratch/rx32940/minION/minIONQC -p 3 -f pdf
done