#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=batch
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
#SBATCH --job-name=$sample\n\
#SBATCH --ntasks=1\n\
#SBATCH --cpus-per-task=1\n\
#SBATCH --time=100:00:00\n\
#SBATCH --mem=10G\n\
#SBATCH --output=%x.%j.out\n\
#SBATCH --error=%x.%j.out\n\
#SBATCH --mail-user=rx32940@uga.edu\n\
#SBATCH --mail-type=ALL\n"

echo -e $header > $code_path/sub.sh

echo -e "echo \"Transcritption analysis for $sample\"" >> $code_path/sub.sh

echo -e "python $code_path/find_feature.py TSS -b $bam -a $REF/GCF_000007685.1_ASM768v1_genomic.gff -o $output/tss\n" >> $code_path/sub.sh

echo -e "python $code_path/find_feature.py TTS -b $bam -a $REF/GCF_000007685.1_ASM768v1_genomic.gff -o $output/tts\n" >> $code_path/sub.sh

echo -e "python $code_path/find_tu.py -b $bam -a $REF/GCF_000007685.1_ASM768v1_genomic.gff -o $output/tu\n" >> $code_path/sub.sh

sbatch $code_path/sub.sh
) & 
wait
done

