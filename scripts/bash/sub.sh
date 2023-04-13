#!/bin/sh
#SBATCH --partition=batch
#SBATCH --job-name=Patoc_Basecalled_Aug_16_2019_Direct-cDNA_PolyATail_tailfindr
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=168:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL

source activate tailfindr
Rscript run_tailfindr.r /scratch/rx32940/minION/Guppy/fast5/Patoc_Basecalled_Aug_16_2019_Direct-cDNA_PolyATail /scratch/rx32940/minION/polyA_cDNA/tailfindr/tailfindr_2021_basecall Patoc_Basecalled_Aug_16_2019_Direct-cDNA_PolyATail.csv 12
conda deactivate
