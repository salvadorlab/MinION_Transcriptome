#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=find_tss
#SBATCH --ntasks=1                    	
#SBATCH --cpus-per-task=1            
#SBATCH --time=100:00:00
#SBATCH --mem=10G
#SBATCH --output=%x.%j.out       
#SBATCH --error=%x.%j.out        
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL 

source activate tu_annotation

code_path="/home/rx32940/github/MinION_Transcriptome/TU_annotation"
bam_path="/scratch/rx32940/minION/polyA_directRNA/map/genome/bam"
output="/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/direct_output"
REF="/scratch/rx32940/minION/polyA_cDNA/map/genome/reference"

cd $output/log


for bam in $bam_path/*.linear.bam;
do
(

sample=$(basename $bam '.linear.bam')
header="#!/bin/bash\n#SBATCH --partition=batch\n\
#SBATCH --job-name=${sample}\n\
#SBATCH --ntasks=1\n\
#SBATCH --cpus-per-task=1\n\
#SBATCH --time=100:00:00\n\
#SBATCH --mem=10G\n\
#SBATCH --output=%x.%j.out\n\
#SBATCH --error=%x.%j.out\n\
#SBATCH --mail-user=rx32940@uga.edu\n\
#SBATCH --mail-type=ALL\n"

echo -e $header > $code_path/sub.sh
echo -e "source activate tu_annotation\n" >> $code_path/sub.sh

echo -e "echo \"Transcritption analysis for $sample\"" >> $code_path/sub.sh

# echo -e "python $code_path/find_feature.py TSS -b $bam -a $REF/GCF_000007685.1_ASM768v1_genomic.gff -o $output/tss -uf 0.9 \n" >> $code_path/sub.sh

# echo -e "python $code_path/find_feature.py TTS -b $bam -a $REF/GCF_000007685.1_ASM768v1_genomic.gff -o $output/tts\n" >> $code_path/sub.sh

echo -e "python $code_path/find_tu.py -b $bam -a $REF/GCF_000007685.1_ASM768v1_genomic.gff -o $output/tu\n" >> $code_path/sub.sh

sbatch $code_path/sub.sh
) & 
wait
done

#####################################################################
### combine features from all samples together
#####################################################################
# ml R/4.1.3-foss-2020b

# code_path="/home/rx32940/github/MinION_Transcriptome/TU_annotation"
# bam_path="/scratch/rx32940/minION/polyA_directRNA/map/genome/bam"
# output="/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/direct_output"
# REF="/scratch/rx32940/minION/polyA_cDNA/map/genome/reference"

# mkdir -p $output/tss/Q29 $output/tss/Q36 $output/tss/combined

# mv $output/tss/Q36_Copen* $output/tss/Q36 

# mv $output/tss/*tab $output/tss/Q29

# Rscript $code_path/combine_samples.r $output/tss/Q29 $output/tss/Q29 tss Q29_dRNA
# Rscript $code_path/combine_samples.r $output/tss/Q36 $output/tss/Q36 tss Q36_dRNA

# #####################################################################
# ### cluster TSS from all samples
# #####################################################################


# # mkdir -p $output/tss/combined/clustered/
# # for i in {0..100..10};
# # do
# # for com in $output/tss/Q*/*.combined.TSS.tab;
# # do
# # sample=$(basename $com '.combined.TSS.tab')
# # # python $code_path/cluster_tss.py -i $com -o $output/tss/$sample.clustered.TSS.tab -c 1 # 10 bp cluster
# # # python $code_path/cluster_tss.py -i $com -o $output/tss/$sample.clustered20.TSS.tab -c 1 -t 20 # 20 bp cluster
# # python $code_path/cluster_tss.py -i $com -o $output/tss/combined/clustered/$sample.clustered$i.TSS.tab -c 1 -t $i # 20 bp cluster
# # echo $sample
# # done
# # done

# #####################################################################
# ### Filter TSS from each cluster
# ### one TSS taken from each cluster (min on lead, max on lag)
# ### get the furthest end of each cluster (max on lead, min on lag)
# #####################################################################

# # based on overlapping size between all samples, TSS within 20bps were used as threshold to be seen as the same tss
# for file in $output/tss/combined/clustered/Q*.clustered20.TSS.tab;
# do
# sample=$(basename $file '.clustered20.TSS.tab')

# python $code_path/filter_tss.py -i $file -o $output/tss/combined/filtered/$sample.filtered20.TSS.tab
# done

# #####################################################################
# ### annotate TSS based on gTSS (pTSS & sTSS) iTSS asTSS oTSS
# #####################################################################

source activate tu_annotation

code_path="/home/rx32940/github/MinION_Transcriptome/TU_annotation"
bam_path="/scratch/rx32940/minION/polyA_directRNA/map/genome/bam"
output="/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/direct_output"
REF="/scratch/rx32940/minION/polyA_cDNA/map/genome/reference"



# combined annotation
for file in $output/tss/combined/filtered/Q*_dRNA.filtered20.TSS.tab;
do
sample=$(basename $file '.filtered20.TSS.tab')
python $code_path/annotate_tss.py \
-i $file \
-g $REF/GCF_000007685.1_ASM768v1_genomic.gff \
-o $output/tss/combined/filtered/$sample.filtered20.TSS.anot.tab \
-s $output/tss/combined/clustered/$sample.clustered20.TSS.tab -c True
done

# individual sample annotation
for file in $output/tss/Q*/*_rna_filtered.linear.TSS.merged.tab;
do
sample=$(basename $file '_rna_filtered.linear.TSS.merged.tab')
dirname=$(dirname $file)
python $code_path/annotate_tss.py \
-i $file \
-g $REF/GCF_000007685.1_ASM768v1_genomic.gff \
-o $dirname/$sample.TSS.anot.tab 
done

# #####################################################################
# ### check overlapping TSS for LIC under 29C and 36C
# #####################################################################


# cat $output/tss/combined/filtered/Q29_dRNA.filtered20.TSS.anot.tab | awk -F'\t' '{for (i=1; i<=14; i++) printf $i"\t"; printf "Q29\n"}' > test1.tab
# cat $output/tss/combined/filtered/Q36_dRNA.filtered20.TSS.anot.tab | tail -n +2 | awk -F'\t' '{for (i=1; i<=14; i++) printf $i"\t"; printf "Q36\n"}' > test.tab

# cat test1.tab test.tab > $output/tss/combined/filtered/combined_dRNA.filtered20.TSS.anot.tab 

# # cluster q29 and q36 tss
# python $code_path/cluster_tss.py -i $output/tss/combined/filtered/combined_dRNA.filtered20.TSS.anot.tab \
# -o $output/tss/combined/filtered/combined_dRNA.filtered20.TSS.anot.clustered.tab \
# -c 1 -t 20 -s 5 # strand on fifth



