#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=get_30bps_around_3end
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=100:00:00
#SBATCH --mem=50G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL


WORK="/scratch/rx32940/minION/polyA_cDNA/map/after_map_end/output"
BED="/scratch/rx32940/minION/polyA_cDNA/map/genome/stats"
ref="/scratch/rx32940/minION/polyA_cDNA/map/genome/reference"

ml BEDTools/2.30.0-GCC-8.3.0

for dir in $BED/*.anot.raw.stats;
do
sample=$(basename $dir ".anot.raw.stats")
mkdir -p $WORK/$sample

cat $BED/$sample.anot.raw.stats | grep -v rRNA > $WORK/$sample/CDS_mapped_only.bed

# create BED file for positions with 30 bps before/after 3' end
# python $WORK/../get_position_after_3_prime.py $BED/$sample.anot.raw.stats $WORK/$sample/position_after_3_prime.bed $WORK/$sample/position_before_3_prime.bed
python $WORK/../get_position_after_3_prime.py $WORK/$sample/CDS_mapped_only.bed $WORK/$sample/position_after_3_prime.CDS.bed $WORK/$sample/position_before_3_prime.CDS.bed


# get sequence from 3' end of all mapping read
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

bedtools getfasta -fi $ref/$REF.fna -bed $WORK/$sample/position_after_3_prime.CDS.bed -fo $WORK/$sample/$sample.after.CDS.fasta
bedtools getfasta -fi $ref/$REF.fna -bed $WORK/$sample/position_before_3_prime.CDS.bed -fo $WORK/$sample/$sample.before.CDS.fasta

done

