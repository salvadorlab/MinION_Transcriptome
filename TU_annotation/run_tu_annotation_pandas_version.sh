#!/bin/sh
#SBATCH --partition=batch
#SBATCH --job-name=cDNA_TTS
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=100:00:00
#SBATCH --mem=30G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL

source activate tu_annotation
ml SAMtools/1.10-GCC-8.3.0


output="/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/output"
input="/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/input"
genefile="GCF_000007685.1_gene_table.txt"
direct_bam="/scratch/rx32940/minION/polyA_directRNA/map/genome/bam"
direct2_bam="/scratch/rx32940/minION/dRNA_Comparison/map/genome/bam"
cdna_bam="/scratch/rx32940/minION/polyA_cDNA/map/genome/bam"
genome_file="GCA_000007685.1_genome_file.txt"
combined="/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/combined"

# # echo -e "sample\tseq_method\tcluster_window\tnum_tts" > $output/TTS/dRNA_samples_window_search_TTS.tab
# # for i in 1 3 5 8 {10..100..10};
# # do
# for file in $direct_bam/*_rna_filtered.linear.bam;
# do
# echo "$file"
# sample=$(basename $file '.linear.bam')
# # # filter out chimeric reads, reads with N in CIGAR means Skipped region from the reference, we grep CIGAR string with N to get only non-chimeric reads
# # samtools view -h -F 4 $file | awk '$6 !~ /N/ || $1 ~ /@/' | samtools view -bS - > $direct_bam/$sample.linear.bam
# # samtools index $direct_bam/$sample.linear.bam

# file1="$direct_bam/LIC_NOPOLYA_rna_filtered.linear.bam" # 90
# file2="$direct_bam/LIC_POLYA_DRNA_CDNA_rna_filtered.linear.bam" # 0
# file3="$direct_bam/LIC_POLYA_rna_filtered.linear.bam" # 0
# file4="$direct_bam/Q29_Copenhageni_Basecalled-June_11_2020_Repeat_Direct-RNA_rna_filtered.linear.bam" # 80
# file5="$direct_bam/Q29_Copenhageni_Basecalled_May_22_2020_Direct-RNA_rna_filtered.linear.bam" # 60
# file6="$direct_bam/Q36_Copenhageni_Basecalled_June_9_2020-Repeat_Direct-RNA_rna_filtered.linear.bam" # 70
# file7="$direct_bam/Q36_Copenhageni_Basecalled_May_31_2020_Direct-RNA_rna_filtered.linear.bam" # 70

# drna_files=($file1 $file2 $file3 $file4 $file5 $file6 $file7)
# window_size=(90 0 0 80 60 70 70)

# for i in {0..6..1};
# do

# echo -e "Current file ${drna_files[$i]};\nCurrent window size ${window_size[$i]}" 

# python tu_annotation.py TSS -o $output/TSS \
# -g $input/$genefile \
# -b ${drna_files[$i]} \
# -c $input/$genome_file \
# -t ${window_size[$i]}
# done

# # echo -e "$sample\tdRNA\t$i\t$(cat $output/TSS/$sample.combined_tss.tab | wc -l)" >> $output/TSS/dRNA_samples_window_search_TSS.tab

# file1="$direct_bam/LIC_NOPOLYA_rna_filtered.linear.bam" # 70
# file2="$direct_bam/LIC_POLYA_DRNA_CDNA_rna_filtered.linear.bam" # 0
# file3="$direct_bam/LIC_POLYA_rna_filtered.linear.bam" # 0
# file4="$direct_bam/Q29_Copenhageni_Basecalled-June_11_2020_Repeat_Direct-RNA_rna_filtered.linear.bam" # 80
# file5="$direct_bam/Q29_Copenhageni_Basecalled_May_22_2020_Direct-RNA_rna_filtered.linear.bam" # 50
# file6="$direct_bam/Q36_Copenhageni_Basecalled_June_9_2020-Repeat_Direct-RNA_rna_filtered.linear.bam" # 60
# file7="$direct_bam/Q36_Copenhageni_Basecalled_May_31_2020_Direct-RNA_rna_filtered.linear.bam" # 60

# drna_files=($file1 $file2 $file3 $file4 $file5 $file6 $file7)
# window_size=(70 0 0 80 50 60 60)

# for i in {0..6..1};
# do
# echo -e "Current file ${drna_files[$i]};\nCurrent window size ${window_size[$i]}" 

# python tu_annotation.py TTS -o $output/TTS \
# -g $input/$genefile \
# -b ${drna_files[$i]} \
# -c $input/$genome_file \
# -t ${window_size[$i]}
# done

# # echo -e "$sample\tdRNA\t$i\t$(cat $output/TTS/$sample.combined_tts.tab | wc -l)" >> $output/TTS/dRNA_samples_window_search_TTS.tab
# done
# # done

# # # echo -e "sample\tseq_method\tcluster_window\tnum_tts" > $output/TTS/cDNA_samples_window_search_TTS.tab
# # # for i in 1 3 5 8 {10..150..10};
# # # do
# for file in $cdna_bam/Copenhageni*_rna_filtered.linear.bam;
# do
# echo "$file"
# sample=$(basename $file '.linear.bam')

# # # # filter out chimeric reads, reads with N in CIGAR means Skipped region from the reference, we grep CIGAR string with N to get only non-chimeric reads
# # # samtools view -h -F 4 $file | awk '$6 !~ /N/ || $1 ~ /@/' | samtools view -bS - > $cdna_bam/$sample.linear.bam
# # # samtools index $cdna_bam/$sample.linear.bam

# file8="$cdna_bam/Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_NoPolyATail_Qiagen_rna_filtered.linear.bam" # 150
# file9="$cdna_bam/Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_NoPolyATail_rna_filtered.linear.bam" # 140
# file10="$cdna_bam/Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_PolyATail_rna_filtered.linear.bam" # 120

# cdna_files=($file8 $file9 $file10)
# window_size=(150 140 120)

# for i in {0..2..1};
# do

# echo -e "Current file ${cdna_files[$i]};\nCurrent window size ${window_size[$i]}" 

# python tu_annotation.py TSS -o $output/TSS \
# -g $input/$genefile \
# -b ${cdna_files[$i]} \
# -c $input/$genome_file \
# -t ${window_size[$i]}
# done

# # echo -e "$sample\tcDNA\t$i\t$(cat $output/TSS/$sample.combined_tss.tab | wc -l)" >> $output/TSS/cDNA_samples_window_search_TSS.tab

file8="$cdna_bam/Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_NoPolyATail_Qiagen_rna_filtered.linear.bam" # 120
file9="$cdna_bam/Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_NoPolyATail_rna_filtered.linear.bam" # 120
file10="$cdna_bam/Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_PolyATail_rna_filtered.linear.bam" # 120

cdna_files=($file8 $file9 $file10)
window_size=(120 120 120)

for i in {0..2..1};
do

echo -e "Current file ${cdna_files[$i]};\nCurrent window size ${window_size[$i]}" 

python tu_annotation.py TTS -o $output/TTS \
-g $input/$genefile \
-b ${cdna_files[$i]} \
-c $input/$genome_file \
-t ${window_size[$i]}

done

# # echo -e "$sample\tcDNA\t$i\t$(cat $output/TTS/$sample.combined_tts.tab | wc -l)" >> $output/TTS/cDNA_samples_window_search_TTS.tab
# done
# # done

# ### combine samples from each sequencing platforms 
# # first argument is input path to the tab files for each inidividual sample,  
# # second argument is feature to combine, can be either tss or tts
# ###################
# ml R/4.1.0-foss-2019b


# Rscript combine_samples.r $output/TSS/ $combined/TSS/raw tss
# Rscript combine_samples.r $output/TTS/ $combined/TTS/raw tts
######################################################
#### cluster TSS identified from all samples (cDNA and dRNA)
######################################################
# for i in {0..300..10};
# do
# for chrom in $combined/TSS/raw/chromI*_raw_TSS.csv;
# do

# sample=$(basename $chrom '.csv')

# python cluster_tss.py -i $chrom -o $combined/TSS/clustered/${sample}_clustered_$i.csv -t $i

# done
# done 

# for i in {0..300..10};
# do
# for chrom in $combined/TTS/raw/chromI*_raw_TTS.csv;
# do

# sample=$(basename $chrom '.csv')

# python cluster_tss.py -i $chrom -o $combined/TTS/clustered/${sample}_clustered_$i.csv -t $i

# done
# done

# # cluster dRNA and cDNA combined TSS- final TSS

for i in {0..150..10};
do

output_file=$output/TSS/cDNA_dRNA_lead_TSS_clustered_$i.csv
python cluster_tss.py -i $output/TSS/cDNA_dRNA_lead_TSS.csv -o $output_file -t $i -c 4
done
# conda deactivate