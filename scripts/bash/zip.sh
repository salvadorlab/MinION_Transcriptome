#!/bin/sh
#SBATCH --partition=batch
#SBATCH --job-name=tar.gz.cDNA
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=100:00:00
#SBATCH --mem=100G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL

cd /scratch/rx32940/minION/Guppy

tar -zvcf PolyA_cDNA_FAST5_Pass.tar.gz PolyA_cDNA_FAST5_Pass/