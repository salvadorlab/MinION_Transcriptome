#!/bin/sh
#SBATCH --partition=batch
#SBATCH --job-name=lortiaTest
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=100:00:00
#SBATCH --mem=100G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL

ml Pysam/0.16.0.1-GCC-8.3.0
ml Biopython/1.75-intel-2019b-Python-3.7.4
ml BEDTools/2.30.0-GCC-8.3.0
ml scipy/1.4.1-intel-2019b-Python-3.7.4

cd /scratch/rx32940/minION/polyA_directRNA/lortia

python LoRTIA/LoRTIA -f True -s poisson ../map/genome/bam/Q29_Copenhageni_Basecalled_May_22_2020_Direct-RNA.sorted.bam ./ ../map/reference/GCF_000007685.1_ASM768v1_genomic.fna
