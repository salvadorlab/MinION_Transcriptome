#!/bin/sh
#SBATCH --partition=batch
#SBATCH --job-name=nanopolishTrans
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=100:00:00
#SBATCH --mem=100G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL

ml Miniconda3/4.9.2
source activate nanopolish

fastq="/scratch/rx32940/minION/merged_data"
bam="/scratch/rx32940/minION/Guppy/map/transcriptome/bam"
ref="/scratch/rx32940/minION/map/reference/GCF_000007685.1_ASM768v1_genomic.fna"
ref_trans="/scratch/rx32940/minION/map/reference/GCF_000007685.1_ASM768v1_cds_from_genomic.fna"
fast5="/scratch/rx32940/minION/Guppy/Direct_RNA_FAST5_Pass"


# nanopolish index --directory=$fast5/Q29_Copenhageni_Direct-RNA_FAST5_Pass/workspace \
# --sequencing-summary=$fast5/Q29_Copenhageni_Direct-RNA_FAST5_Pass/sequencing_summary.txt \
# $fast5/Q29_Copenhageni_Direct-RNA_FAST5_Pass/Q29_Copenhageni_Direct-RNA_FAST5_Pass.fastq.gz

# nanopolish index --directory=$fast5/Q29_Copenhageni_Direct-RNA_Repeat_FAST5_Pass/workspace \
# --sequencing-summary=$fast5/Q29_Copenhageni_Direct-RNA_Repeat_FAST5_Pass/sequencing_summary.txt \
# $fast5/Q29_Copenhageni_Direct-RNA_Repeat_FAST5_Pass/Q29_Copenhageni_Direct-RNA_Repeat_FAST5_Pass.fastq.gz

# nanopolish index --directory=$fast5/Q36_Copenhageni_Direct-RNA_FAST5_Pass/workspace \
# --sequencing-summary=$fast5/Q36_Copenhageni_Direct-RNA_FAST5_Pass/sequencing_summary.txt \
# $fast5/Q36_Copenhageni_Direct-RNA_FAST5_Pass/Q36_Copenhageni_Direct-RNA_FAST5_Pass.fastq.gz

# nanopolish index --directory=$fast5/Q36_Copenhageni_Direct-RNA_Repeat_FAST5_Pass/workspace \
# --sequencing-summary=$fast5/Q36_Copenhageni_Direct-RNA_Repeat_FAST5_Pass/sequencing_summary.txt \
# $fast5/Q36_Copenhageni_Direct-RNA_Repeat_FAST5_Pass/Q36_Copenhageni_Direct-RNA_Repeat_FAST5_Pass.fastq.gz

# for sample in $fastq/*PolyATail.fastq;
# do
# nanopolish index --directory=$fast5/PolyA_cDNA_FAST5_Pass \
# $sample
# done

for sample in $fast5/*;
do

sample_name=$(basename $sample)

nanopolish polya --reads=$sample/$sample_name.fastq.gz --bam=$bam/$sample_name.trans.sorted.bam --genome=$ref_trans > nanopolish_output/$sample_name.nanopolish.trans.polyA.tsv

done


conda deactivate