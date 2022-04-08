#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_PolyATail_rna_filtered
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=100:00:00
#SBATCH --mem=10G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL

source activate tu_annotation

echo "Transcritption analysis for Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_PolyATail_rna_filtered"
python /home/rx32940/github/MinION_Transcriptome/TU_annotation/find_feature.py TSS -b /scratch/rx32940/minION/polyA_cDNA/map/genome/bam/Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_PolyATail_rna_filtered.linear.bam -a /scratch/rx32940/minION/polyA_cDNA/map/genome/reference/GCF_000007685.1_ASM768v1_genomic.gff -o /scratch/rx32940/minION/polyA_directRNA/TU_Annotation/cdna_output/tss/raw -uf 0.9 

python /home/rx32940/github/MinION_Transcriptome/TU_annotation/find_feature.py TTS -b /scratch/rx32940/minION/polyA_cDNA/map/genome/bam/Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_PolyATail_rna_filtered.linear.bam -a /scratch/rx32940/minION/polyA_cDNA/map/genome/reference/GCF_000007685.1_ASM768v1_genomic.gff -o /scratch/rx32940/minION/polyA_directRNA/TU_Annotation/cdna_output/tts

python /home/rx32940/github/MinION_Transcriptome/TU_annotation/find_tu.py -b /scratch/rx32940/minION/polyA_cDNA/map/genome/bam/Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_PolyATail_rna_filtered.linear.bam -a /scratch/rx32940/minION/polyA_cDNA/map/genome/reference/GCF_000007685.1_ASM768v1_genomic.gff -o /scratch/rx32940/minION/polyA_directRNA/TU_Annotation/cdna_output/tu

