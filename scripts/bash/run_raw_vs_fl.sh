#!/bin/sh
#SBATCH --partition=batch
#SBATCH --job-name=find_raw_fl_diff
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=100:00:00
#SBATCH --mem=30G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL

MAIN="/scratch/rx32940/minION/polyA_cDNA/raw_vs_fl"
RAW="/scratch/rx32940/minION/polyA_cDNA/fastq"
FL="/scratch/rx32940/minION/polyA_cDNA/pychopper/output"

############### 1) get read id from both raw and fl cDNA reads###############
# for file in $RAW/*fastq;
# do
# sample=$(basename $file ".fastq")
# cat $file | grep "^@" | awk '{print substr($1,2)}'> $MAIN/raw_reads/$sample.txt
# done

# for file in $FL/*/*.final.fastq;
# do
# sample=$(basename $file ".final.fastq")
# cat $file | grep -E "^@[0-9]{1,10}\:[0-9]{1,10}\|" | awk -F'[| ]' '{print $2}' > $MAIN/fl_reads/$sample.txt
# done

############### 2) grep raw reads filtered during FL identification###############
for file in $MAIN/raw_reads/*.txt;
do
sample=$(basename $file ".txt")
echo $sample
diff <(sort $MAIN/raw_reads/$sample.txt) <(sort $MAIN/fl_reads/$sample.txt) --suppress-common-line -y | awk '{print $1}' | grep -v ">" > $MAIN/raw_reads_filtered/$sample.txt
done

############### 3) raw/FL cDNA reads mapped to their corresponding reference genome###############
# raw cDNA reads mapped to genome 
RBED="/scratch/rx32940/minION/polyA_cDNA/map/genome/bed"
FBED="/scratch/rx32940/minION/polyA_cDNA/map_full/genome/bed"

for bed in $RBED/*bed;
do
sample=$(basename $bed ".bed")
echo $sample
awk '{print $4}' $bed | sort | uniq > $MAIN/map_raw/$sample.txt

done

## FL cDNA reads mapped to genome 
for bed in $FBED/*.withClip.bed;
do
sample=$(basename $bed ".withClip.bed")
echo $sample
awk '{print $4}' $bed | awk -F"|" '{print $2}'| sort | uniq > $MAIN/map_fl/$sample.txt

done

############### 4) reads mapped by raw cDNA reads but not FL reads (and the raw reads was not filtered out due to not identified as FL)###############

for file in $MAIN/map_raw/*.txt;
do
sample=$(basename $file ".txt")
echo $sample

# 1) reads mapped by raw but not fl reads
diff <(sort $MAIN/map_raw/$sample.txt) <(sort $MAIN/map_fl/$sample.txt) --suppress-common-line -y | awk '{print $1}' | grep -v ">" > $MAIN/map_raw_vs_fl/$sample.tmp.txt
# 2) raw reads mapped only (also could have FL reads identified)
diff <(sort $MAIN/map_raw_vs_fl/$sample.tmp.txt) <(sort $MAIN/raw_reads_filtered/$sample.txt) --suppress-common-line -y | awk '{print $1}' | grep -v ">" > $MAIN/map_raw_vs_fl/$sample.txt
rm $MAIN/map_raw_vs_fl/$sample.tmp.txt
done
