#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=long_reads_bed-FL
#SBATCH --ntasks=1                      
#SBATCH --cpus-per-task=1           
#SBATCH --time=100:00:00
#SBATCH --mem=50G
#SBATCH --output=%x.%j.out       
#SBATCH --error=%x.%j.out        
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL 

ml SAMtools/1.14-GCC-8.3.0
ml BEDTools/2.30.0-GCC-8.3.0
ml Pysam/0.16.0.1-iccifort-2020.4.304

MAIN="/scratch/rx32940/minION/polyA_cDNA/longGap"
BAM="/scratch/rx32940/minION/polyA_cDNA/map/genome/bam"
FLBAM="/scratch/rx32940/minION/polyA_cDNA/map_full/genome/bam"
# # get read id with "N" first (BED file obtained from R script)
# for bed in $MAIN/long_read_list/*.FL.bed;
# do
# sample=$(basename $bed ".FL.bed")
# cat $bed | awk '{print $4}' | tail -n +2 > $MAIN/long_read_list/$sample.FL.txt

# done

# # raw cDNA longGap reads
# for bam in $BAM/*.sorted.bam;
# do
# sample=$(basename $bam ".sorted.bam")
# echo $sample

# time python $MAIN/extract_from_bam.py -b $bam -n $MAIN/long_read_list/$sample.txt -o $MAIN/bam/$sample.LG.bam
# time bedtools bamtobed -split -i $MAIN/bam/$sample.LG.bam > $MAIN/bed/$sample.LG.bed
# done

# FL cDNA longGap reads
for bam in $FLBAM/*.withClip.sorted.bam;
do
sample=$(basename $bam ".withClip.sorted.bam")
echo $sample

time python $MAIN/extract_from_bam.py -b $bam -n $MAIN/long_read_list/$sample.FL.txt -o $MAIN/bam/$sample.LG.FL.bam
time bedtools bamtobed -split -i $MAIN/bam/$sample.LG.FL.bam > $MAIN/bed/$sample.LG.FL.bed
done