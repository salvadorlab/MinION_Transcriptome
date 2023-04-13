#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=htseq_count
#SBATCH --ntasks=1                      
#SBATCH --cpus-per-task=1           
#SBATCH --time=100:00:00
#SBATCH --mem=50G
#SBATCH --output=%x.%j.out       
#SBATCH --error=%x.%j.out        
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL 

# ml HTSeq/0.9.1-intel-2019b-Python-2.7.16
source activate tu_annotation


BAM="/scratch/rx32940/minION/polyA_cDNA/map/genome/bam"
ref="/scratch/rx32940/minION/polyA_cDNA/map/genome/reference"
OUT="/scratch/rx32940/minION/polyA_cDNA/map/htseq/out"


for file in $BAM/*sorted.bam;
do

sample=$(basename $file ".sorted.bam")
echo $sample

    if [[ $sample == Copen* ]]
    then
        REF="GCF_000007685.1_ASM768v1_genomic";
    elif [[ $sample == Man* ]] || [[ $sample == Icter* ]]
    then
        REF="GCF_014858895.1_ASM1485889v1_genomic"
    elif [[ $sample == Pa* ]] 
    then
        REF="GCF_000017685.1_ASM1768v1_genomic"
    else
        REF=""
    fi

# no strandness for reads in this data
htseq-count --nonunique all --type="CDS" --idattr="ID" -f bam -s no --additional-attr="product" $file $ref/$REF.gff > $OUT/$sample.txt
done