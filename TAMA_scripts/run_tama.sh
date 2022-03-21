#!/bin/sh
#SBATCH --partition=batch
#SBATCH --job-name=tama_merge
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=100:00:00
#SBATCH --mem=30G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL


source activate tama
ml SAMtools/1.10-GCC-8.3.0

MAIN="/scratch/rx32940/minION/polyA_directRNA/TAMA"
TAMA="$MAIN/tama"
cdna_bam="/scratch/rx32940/minION/polyA_cDNA/map/genome/bam"
REF="/scratch/rx32940/minION/polyA_cDNA/map/genome/reference"
direct_bam="/scratch/rx32940/minION/polyA_directRNA/map/genome/bam"


for file in $cdna_bam/Copen*_rna_filtered.bam;
do
sample=$(basename $file '.bam')
echo $sample
# samtools view -h -o $cdna_bam/$sample.sam $file

python $TAMA/tama_collapse.py -s $cdna_bam/$sample.sam \
-f $REF/GCF_000007685.1_ASM768v1_genomic.fna \
-p $MAIN/output/$sample \
-x no_cap \
-icm ident_map \
-c 70 # cDNA reads has low fraction of reads mapping reference, change coverage fraction from 99 -> 40
done

for file in $direct_bam/*_rna_filtered.bam;
do
sample=$(basename $file '.bam')
echo $sample
# samtools view -h -o $direct_bam/$sample.sam $file

python $TAMA/tama_collapse.py -s $direct_bam/$sample.sam \
-f $REF/GCF_000007685.1_ASM768v1_genomic.fna \
-p $MAIN/output/$sample \
-x no_cap \
-icm ident_map -c 90
done

################################################################
#### merge all samples from same sequencing platform together
################################################################

### dRNA
python $TAMA/tama_merge.py -f $MAIN/input/drna_merge_input.tab -p $MAIN/merged/dRNA_all_samples 

### cDNA
python $TAMA/tama_merge.py -f $MAIN/input/cdna_merge_input.tab -p $MAIN/merged/cDNA_all_samples 


