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

PYCHOPPER="/scratch/rx32940/minION/polyA_cDNA/pychopper/output"
WORK="/scratch/rx32940/minION/polyA_cDNA/vnp/out"
BED="/scratch/rx32940/minION/polyA_cDNA/map_full/stats"
ref="/scratch/rx32940/minION/polyA_cDNA/map_full/reference"

ml BEDTools/2.30.0-GCC-8.3.0

for dir in $PYCHOPPER/*/;
do
sample=$(basename $dir)
mkdir -p $WORK/$sample

# from pychopper's output, get all reads with VNP primer hit
cat $dir/aln_hits_DCS109.bed $dir/aln_hits_PCB109.bed | grep VNP | grep -v $'VNP\t0' > $WORK/$sample/vnp_aln_hit.bed

# get the read id's for reads with VNP primer hit
awk '{print $1}' $WORK/$sample/vnp_aln_hit.bed > $WORK/$sample/vnp_read_id.txt

# from alignment, get the position of where reads with VNP hits were mapped to (*.anot.stats file from anot_read_transcript.py)
grep -F -f $WORK/$sample/vnp_read_id.txt $BED/$sample.anot.stats > $WORK/$sample/vnp_read_mapped.bed

# create BED file for positions with 30 bps before/after 3' end
python $WORK/../get_position_after_3_prime.py $WORK/$sample/vnp_read_mapped.bed $WORK/$sample/position_after_3_prime.bed $WORK/$sample/position_before_3_prime.bed

cat $WORK/$sample/position_after_3_prime.bed | grep $'\t+\t' | grep -e "partial(3')" > $WORK/$sample/after_3_prime.forward.bed
cat $WORK/$sample/position_after_3_prime.bed | grep $'\t-\t' | grep -e "partial(3')"  > $WORK/$sample/after_3_prime.reverse.bed

cat $WORK/$sample/position_before_3_prime.bed |  grep $'\t+\t' | grep -e "partial(3')" > $WORK/$sample/before_3_prime.forward.bed
cat $WORK/$sample/position_before_3_prime.bed |  grep $'\t-\t' | grep -e "partial(3')"  > $WORK/$sample/before_3_prime.reverse.bed


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

bedtools getfasta -fi $ref/$REF.fna -bed $WORK/$sample/after_3_prime.forward.bed -fo $WORK/$sample/$sample.after.forward.fasta
bedtools getfasta -fi $ref/$REF.fna -bed $WORK/$sample/after_3_prime.reverse.bed -fo $WORK/$sample/$sample.after.reverse.fasta

bedtools getfasta -fi $ref/$REF.fna -bed $WORK/$sample/before_3_prime.forward.bed -fo $WORK/$sample/$sample.before.forward.fasta
bedtools getfasta -fi $ref/$REF.fna -bed $WORK/$sample/before_3_prime.reverse.bed -fo $WORK/$sample/$sample.before.reverse.fasta

done

