#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=find_tss_cdna
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
bam_path="/scratch/rx32940/minION/polyA_cDNA/map/genome/bam"
output="/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/cdna_output"
REF="/scratch/rx32940/minION/polyA_cDNA/map/genome/reference"

cd $output/log


for bam in $bam_path/Copen*.linear.bam;
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

echo -e "python $code_path/find_feature.py TSS -b $bam -a $REF/GCF_000007685.1_ASM768v1_genomic.gff -o $output/tss/raw -uf 0.9 \n" >> $code_path/sub.sh

echo -e "python $code_path/find_feature.py TTS -b $bam -a $REF/GCF_000007685.1_ASM768v1_genomic.gff -o $output/tts\n" >> $code_path/sub.sh

echo -e "python $code_path/find_tu.py -b $bam -a $REF/GCF_000007685.1_ASM768v1_genomic.gff -o $output/tu\n" >> $code_path/sub.sh

sbatch $code_path/sub.sh
) & 
wait
done

#####################################################################
### combine features from all samples together
#####################################################################
ml R/4.1.3-foss-2020b

code_path="/home/rx32940/github/MinION_Transcriptome/TU_annotation"
bam_path="/scratch/rx32940/minION/polyA_cDNA/map/genome/bam"
output="/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/cdna_output"
REF="/scratch/rx32940/minION/polyA_cDNA/map/genome/reference"

mkdir -p $output/tss/combined

# mv $output/tss/Q36_Copen* $output/tss/Q36 

Rscript $code_path/combine_samples.r $output/tss/raw $output/tss/combined tss LIC_cDNA

# #####################################################################
# ### cluster TSS from all samples
# #####################################################################


python $code_path/cluster_tss.py -i $output/tss/combined/LIC_cDNA.combined.TSS.tab -o $output/tss/combined/LIC_cDNA.clustered20.TSS.tab -c 1 -t 20

# #####################################################################
# ### Filter TSS from each cluster
# ### one TSS taken from each cluster (min on lead, max on lag)
# ### get the furthest end of each cluster (max on lead, min on lag)
# #####################################################################

# ### based on overlapping size between all samples, TSS within 20bps were used as threshold to be seen as the same tss
python $code_path/filter_tss.py -i $output/tss/combined/LIC_cDNA.clustered20.TSS.tab -o $output/tss/combined/LIC_cDNA.filtered20.TSS.tab


# #####################################################################
# ### annotate TSS based on gTSS (pTSS & sTSS) iTSS asTSS oTSS
# #####################################################################

source activate tu_annotation

code_path="/home/rx32940/github/MinION_Transcriptome/TU_annotation"
bam_path="/scratch/rx32940/minION/polyA_cDNA/map/genome/bam"
output="/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/cdna_output"
REF="/scratch/rx32940/minION/polyA_cDNA/map/genome/reference"

###########combined annotation
python $code_path/annotate_tss.py \
-i $output/tss/combined/LIC_cDNA.filtered20.TSS.tab \
-g $REF/GCF_000007685.1_ASM768v1_genomic.gff \
-o $output/tss/combined/LIC_cDNA.anot.TSS.tab \
-s $output/tss/combined/LIC_cDNA.clustered20.TSS.tab \
-c True


# ########individual sample annotation
# python $code_path/annotate_tss.py \
# -i $output/tss/raw/Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_NoPolyATail_Qiagen_rna_filtered.linear.TSS.merged.tab \
# -g $REF/GCF_000007685.1_ASM768v1_genomic.gff \
# -o $output/tss/anot/Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_NoPolyATail_Qiagen_rna_filtered.linear.TSS.anot.tab

# #####################################################################
# ### check overlapping TSS for LIC under 29C and 36C
# #####################################################################


# cat $output/tss/combined/filtered/Q29_dRNA.filtered30.TSS.anot.tab | awk -F'\t' '{for (i=1; i<=14; i++) printf $i"\t"; printf "Q29\n"}' > test1.tab
# cat $output/tss/combined/filtered/Q36_dRNA.filtered30.TSS.anot.tab | tail -n +2 | awk -F'\t' '{for (i=1; i<=14; i++) printf $i"\t"; printf "Q36\n"}' > test.tab

# cat test1.tab test.tab > $output/tss/combined/filtered/combined_dRNA.filtered30.TSS.anot.tab 

# # cluster q29 and q36 tss
# python $code_path/cluster_tss.py -i $output/tss/combined/filtered/combined_dRNA.filtered30.TSS.anot.tab \
# -o $output/tss/combined/filtered/combined_dRNA.filtered30.TSS.anot.clustered.tab \
# -c 1 -t 30 -s 5 # strand on fifth


