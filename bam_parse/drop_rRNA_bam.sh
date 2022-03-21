#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=droprRNAcDNA
#SBATCH --ntasks=1                      
#SBATCH --cpus-per-task=3            
#SBATCH --time=100:00:00
#SBATCH --mem=10G
#SBATCH --output=%x.%j.out       
#SBATCH --error=%x.%j.out        
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL

ml BEDOPS/2.4.39-foss-2019b
ml SAMtools/1.10-GCC-8.3.0
ml BEDTools/2.30.0-GCC-8.3.0
ml Qualimap/2.2.1-foss-2019b-R-3.6.2

REF="/scratch/rx32940/minION/polyA_cDNA/map/genome/reference"
BAM="/scratch/rx32940/minION/polyA_cDNA/map/genome/bam"
QUAL="/scratch/rx32940/minION/polyA_cDNA/map/genome/qualimap"
STATS="/scratch/rx32940/minION/polyA_cDNA/map/genome/sam_stats"
FLAG="/scratch/rx32940/minION/polyA_cDNA/map/genome/flagstat"

# # convert gff to bed
# sortBed -i $REF/GCF_000007685.1_ASM768v1_genomic.gff | gff2bed > $REF/GCF_000007685.1_ASM768v1_genomic.bed
# sortBed -i $REF/GCF_014858895.1_ASM1485889v1_genomic.gff | gff2bed > $REF/GCF_014858895.1_ASM1485889v1_genomic.bed
# sortBed -i $REF/GCF_000017685.1_ASM1768v1_genomic.gff | gff2bed > $REF/GCF_000017685.1_ASM1768v1_genomic.bed

# # # grep 16S and 23S rRNA and all tRNA coding region from bed file (other rRNA ex. 30S and 50S rRNA are included in the CDS file)
# cat $REF/GCF_000007685.1_ASM768v1_genomic.bed | grep 'ID=rna-LIC_'  > $REF/GCF_000007685.1_ASM768v1_rRNA_region.bed
# cat $REF/GCF_014858895.1_ASM1485889v1_genomic.bed | grep 'ID=rna-LeptoLang_'  > $REF/GCF_014858895.1_ASM1485889v1_rRNA_region.bed
# cat $REF/GCF_000017685.1_ASM1768v1_genomic.bed | grep 'ID=rna-LEPBI_'  > $REF/GCF_000017685.1_ASM1768v1_rRNA_region.bed

# # filter out reads mapped to either 16S rRNA region or 23S rRNA region
# for file in $BAM/Copen*.sorted.bam;
# do
# sample=$(basename $file ".sorted.bam")
# echo $sample

# samtools view \
# -L $REF/GCF_000007685.1_ASM768v1_rRNA_region.bed \
# -U $BAM/${sample}_rna_filtered.bam \
# -o /dev/null \
# $file


# mkdir -p $QUAL/${sample}_rna_filtered
# qualimap bamqc -bam $BAM/${sample}_rna_filtered.bam -nw 5000 -nt 3 -c -outdir $QUAL/${sample}_rna_filtered -gff $REF/GCF_000007685.1_ASM768v1_genomic.gtf
# samtools stats -d -@ 3 $BAM/${sample}_rna_filtered.bam > $STATS/${sample}_rna_filtered.stats
# samtools flagstats $BAM/${sample}_rna_filtered.bam  > $FLAG/${sample}_rna_filtered.flagstat
# done

# Icterohaemorrhagiae
for file in $BAM/Ictero*.sorted.bam;
do
sample=$(basename $file ".sorted.bam")
echo $sample

samtools view \
-L $REF/GCF_014858895.1_ASM1485889v1_rRNA_region.bed \
-U $BAM/${sample}_rna_filtered.bam \
-o /dev/null \
$file

mkdir -p $QUAL/${sample}_rna_filtered
qualimap bamqc -bam $BAM/${sample}_rna_filtered.bam -nw 5000 -nt 3 -c -outdir $QUAL/${sample}_rna_filtered -gff $REF/GCF_014858895.1_ASM1485889v1_genomic.gtf
samtools stats -d -@ 3 $BAM/${sample}_rna_filtered.bam > $STATS/${sample}_rna_filtered.stats
samtools flagstats $BAM/${sample}_rna_filtered.bam  > $FLAG/${sample}_rna_filtered.flagstat
done

# Mankarso
for file in $BAM/Man*.sorted.bam;
do
sample=$(basename $file ".sorted.bam")
echo $sample

samtools view \
-L $REF/GCF_014858895.1_ASM1485889v1_rRNA_region.bed \
-U $BAM/${sample}_rna_filtered.bam \
-o /dev/null \
$file

mkdir -p $QUAL/${sample}_rna_filtered
qualimap bamqc -bam $BAM/${sample}_rna_filtered.bam -nw 5000 -nt 3 -c -outdir $QUAL/${sample}_rna_filtered -gff $REF/GCF_014858895.1_ASM1485889v1_genomic.gtf
samtools stats -d -@ 3 $BAM/${sample}_rna_filtered.bam > $STATS/${sample}_rna_filtered.stats
samtools flagstats $BAM/${sample}_rna_filtered.bam  > $FLAG/${sample}_rna_filtered.flagstat
done

# patoc
for file in $BAM/Patoc*.sorted.bam;
do
sample=$(basename $file ".sorted.bam")
echo $sample

samtools view \
-L $REF/GCF_000017685.1_ASM1768v1_rRNA_region.bed \
-U $BAM/${sample}_rna_filtered.bam \
-o /dev/null \
$file

mkdir -p $QUAL/${sample}_rna_filtered
qualimap bamqc -bam $BAM/${sample}_rna_filtered.bam -nw 5000 -nt 3 -c -outdir $QUAL/${sample}_rna_filtered -gff $REF/GCF_000017685.1_ASM1768v1_genomic.gtf
samtools stats -d -@ 3 $BAM/${sample}_rna_filtered.bam > $STATS/${sample}_rna_filtered.stats
samtools flagstats $BAM/${sample}_rna_filtered.bam  > $FLAG/${sample}_rna_filtered.flagstat
done