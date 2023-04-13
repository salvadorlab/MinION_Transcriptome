#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=get_trim_pos
#SBATCH --ntasks=1                      
#SBATCH --cpus-per-task=1           
#SBATCH --time=100:00:00
#SBATCH --mem=50G
#SBATCH --output=%x.%j.out       
#SBATCH --error=%x.%j.out        
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL 

OUT="/scratch/rx32940/minION/polyA_cDNA/pychopper/output"
FASTQ="/scratch/rx32940/minION/polyA_cDNA/fastq"
ml seqtk/1.3-GCC-8.3.0

# readid kept after pychopper (this step can replace by "get position where each read was trimmed at" below)
for file in $OUT/*; 
do
    sample=$(basename $file)
    echo $sample    
  cat $file/$sample.final.fastq | grep '^@[0-9]' | awk '{print $1}' | awk -F'|' '{print $2}' | grep . > $file/kept_readid.txt; 
done

# grep kept read id in aln_hits_combined.bed (combined from aln_hits_PCB109.bed & aln_hits_PCB109.bed )
for file in $OUT/*; 
do
    echo $(basename $file)
    sort $file/kept_readid.txt > $file/tmp.txt
    sort $file/aln_hits_combined.bed | uniq > $file/tmp2.txt
    grep -F -f $file/tmp.txt $file/tmp2.txt > $file/aln_hits_kept.bed; 
    rm $file/tmp.txt $file/tmp2.txt
done

# # get position where each read was trimmed at
# for file in $OUT/*; 
# do

# sample=$(basename $file)
# echo $file
# echo $sample
#  cat $file/$sample.final.fastq | grep -E "^@[0-9]{1,10}\:[0-9]{1,10}\|" | awk '{print $1}' | awk -F'[@:|]' '{print $4"\t"$2"\t"$3}' > $file/trim_pos.bed; 
# done


# # add the total length of each read to the trim position bed file (* this trask takes really long time)
for file in $OUT/*;
do

sample=$(basename $file)
echo $sample

cat $file/trim_pos.bed | awk '{print $1}' > $file/trim_list.txt

seqtk subseq $FASTQ/$sample.fastq $file/trim_list.txt > $file/trimmed.filtered.fastq

cat $file/trimmed.filtered.fastq | awk 'NR%4==2{ print length($0) }' > $file/tmp.txt

cat $file/trimmed.filtered.fastq | awk 'NR%4==1{ print substr($1,2) }' > $file/tmp2.txt

paste $file/tmp2.txt $file/tmp.txt | sort > $file/tmp3.txt

paste $file/tmp3.txt <(sort $file/trim_pos.bed) | awk '{print $1"\t"$4"\t"$5"\t"$2}' > $file/trim_pos_sorted.bed

rm $file/trim_list.txt $file/trimmed.filtered.fastq $file/tmp.txt $file/tmp2.txt $file/tmp3.txt 

done