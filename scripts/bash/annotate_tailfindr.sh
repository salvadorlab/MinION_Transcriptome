#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=Anntailfindr
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=168:00:00
#SBATCH --mem-per-cpu=30G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL


#################################################################
#
# 1) This file was made to match Tailfindr output for each reads 
# to genome and transcriptome alignment (sam file) of the reads 
# for annotation

# 2) check the sample name, all sample name for each file need to be captialized 
# for this script to work

# 3) use tailfindr_stats.sh file to obtain statistics of the output files
#
#####################################################################

source activate tailfindr

ALLSAMPLES="Copenhageni Icterohaemorrhagiae Patoc Mankarso"


GENOME="genome transcriptome"
POLYA="noPolyA PolyA"
MAP="/scratch/rx32940/minION/polyA_cDNA/map"
TAILFINDR="/scratch/rx32940/minION/polyA_cDNA/tailfindr/tailfindr_multiplex_output"

for SAMPLE in $ALLSAMPLES;
do

    for polya in $POLYA;
    do 

        if [ $polya = "noPolyA" ]
        then
        NAME="NoPolyATail"
        else
        NAME="PolyATail" 
        fi

        for file in $GENOME;
        do
            ABB="."
            if [ $file = "transcriptome" ]
            then
            ABB=".trans." 
            fi

        echo $MAP/$file/sam/${SAMPLE}_Basecalled_Aug_16_2019_Direct-cDNA_${NAME}.sam
        
        Rscript annotate_tailfindr.r $MAP/$file/sam/${SAMPLE}_Basecalled_Aug_16_2019_Direct-cDNA_${NAME}.sam \
        $TAILFINDR/cDNA_PolyATail.csv \
        $TAILFINDR/annotated/${SAMPLE}_$polya.annotated${ABB}csv

        done

    done

done

conda deactivate
