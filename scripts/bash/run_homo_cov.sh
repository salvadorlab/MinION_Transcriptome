#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=get_nuc5_cov_cDNA-raw
#SBATCH --ntasks=1                      
#SBATCH --cpus-per-task=1           
#SBATCH --time=100:00:00
#SBATCH --mem=50G
#SBATCH --output=%x.%j.out       
#SBATCH --error=%x.%j.out        
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL 

MAIN="/scratch/rx32940/minION/homopolymer_cov_genome"
REF="/scratch/rx32940/minION/polyA_cDNA/map/genome/reference/"

##########################################################
#
# get regions in reference genome with homopolymers
# this script will return the region +- 100 from the homopolymer region
# if the exceed region of transcript, then the begining and end of transcript position will be used
# the 4th column has the exact region of the homopolymer
# the 5th column has the length of the homopolymer
#
#########################################################

# ml Biopython/1.78-intel-2020b
# # arg1: nucleotide to find, arg2: min length of the homopolymers, arg3: input reference fasta, arg4: output bed file
# for nuc in A T C G;
# do
# # python $MAIN/get_homo_coord.py $nuc 5 $REF/GCF_000007685.1_ASM768v1_genomic.fna $MAIN/bed/GCF_000007685.1_ASM768v1_homo${nuc}_5.bed
# python $MAIN/isHomo_in_cDNA.py $MAIN/bed/GCF_000007685.1_ASM768v1_homo${nuc}_5.bed $REF/GCF_000007685.1_ASM768v1_genomic.gff $MAIN/bed/GCF_000007685.1_ASM768v1_homo${nuc}_5.anot.bed

# # python $MAIN/get_homo_coord.py $nuc 5 $REF/GCF_000017685.1_ASM1768v1_genomic.fna $MAIN/bed/GCF_000017685.1_ASM1768v1_homo${nuc}_5.bed
# python $MAIN/isHomo_in_cDNA.py $MAIN/bed/GCF_000017685.1_ASM1768v1_homo${nuc}_5.bed $REF/GCF_000017685.1_ASM1768v1_genomic.gff $MAIN/bed/GCF_000017685.1_ASM1768v1_homo${nuc}_5.anot.bed


# # python $MAIN/get_homo_coord.py $nuc 5 $REF/GCF_014858895.1_ASM1485889v1_genomic.fna $MAIN/bed/GCF_014858895.1_ASM1485889v1_homo${nuc}_5.bed
# python $MAIN/isHomo_in_cDNA.py $MAIN/bed/GCF_014858895.1_ASM1485889v1_homo${nuc}_5.bed $REF/GCF_014858895.1_ASM1485889v1_genomic.gff $MAIN/bed/GCF_014858895.1_ASM1485889v1_homo${nuc}_5.anot.bed

# done

# # get only CDS region from GFF converted BED file
# cat $MAIN/reference/GCF_000007685.1_ASM768v1_genomic.bed | grep "cds-" > $MAIN/reference/GCF_000007685.1_ASM768v1_genomic.CDS.bed
# cat $MAIN/reference/GCF_014858895.1_ASM1485889v1_genomic.bed | grep "cds-" > $MAIN/reference/GCF_014858895.1_ASM1485889v1_genomic.CDS.bed
# cat $MAIN/reference/GCF_000017685.1_ASM1768v1_genomic.bed | grep "cds-" > $MAIN/reference/GCF_000017685.1_ASM1768v1_genomic.CDS.bed


##########################################################
#
# get coverage from every base around homopolymer region
#
#########################################################
# ml Anaconda3/2022.05
# ml BEDTools/2.30.0-GCC-8.3.0

# BAM="/scratch/rx32940/minION/polyA_directRNA/map/genome/bam"
# OUT="/scratch/rx32940/minION/homopolymer_cov_genome/output"
# lower_limit=5 # at least how many nucleotide to count as homopolymer

# ml BEDTools/2.30.0-GCC-8.3.0

# mkdir -p $OUT/overall_trans_cov

# for HOMO in A T C G;
# do
#     echo "current homo$HOMO"
#     mkdir -p $OUT/nuc_$lower_limit/homo$HOMO/
#     for file in $BAM/*sorted.bam;
#     do
#     sample=$(basename $file ".sorted.bam")
#     echo $sample

#     REF="GCF_000007685.1_ASM768v1_genomic"
#     REF2="GCF_000007685.1_ASM768v1"

#     echo "      - get overall coverage within CDS region"
#     # # get overall coverage for each transcript
#     bedtools coverage -s -a $MAIN/reference/$REF.CDS.bed -b $file -d > $OUT/overall_trans_cov/$sample.cov

#     echo "      - get max coverage within each CDS region"
#     # # get max cov from each transcript
#     python trans_max_cov_pos.py $OUT/overall_trans_cov/$sample.cov $OUT/overall_trans_cov/$sample.max.cov

#     echo "      - only get homo region identified with CDS region"
#     cat $MAIN/bed/${REF2}_homo${HOMO}_$lower_limit.anot.bed | grep -v "noncoding" > $MAIN/bed/${REF2}_homo${HOMO}_$lower_limit.anot.CDS.bed

#     echo "      - get coverage 100bps within CDS homo region"
#     # # get coverage around each homopolymer region
#     bedtools coverage -s -a $MAIN/bed/${REF2}_homo${HOMO}_$lower_limit.anot.CDS.bed -b $file -d > $OUT/nuc_$lower_limit/homo$HOMO/$sample.$HOMO$lower_limit.CDS.cov

#     echo "      - get relative coverage 100bps within CDS homo region (based on maximum coverage of each CDS region)"
#     # # get relative cov within the flanking regions of homopolymer by divide each base' cov with the max cov of the transcript
#     python relative_cov_per_base.py $OUT/overall_trans_cov/$sample.max.cov $OUT/nuc_$lower_limit/homo$HOMO/$sample.$HOMO$lower_limit.CDS.cov $OUT/nuc_$lower_limit/homo$HOMO/$sample.$HOMO$lower_limit.relative.CDS.cov True

#     done
# done

# # ######## cDNA #################################
# ml Anaconda3/2022.05
# ml BEDTools/2.30.0-GCC-8.3.0

# BAM="/scratch/rx32940/minION/polyA_cDNA/map/genome/bam" #if map_full also need to change search pattern to *.withClip.sorted.bam 
# OUT="/scratch/rx32940/minION/homopolymer_cov_genome/output"
# lower_limit=5 # at least how many nucleotide to count as homopolymer

# mkdir -p $OUT/overall_trans_cov

# for HOMO in A T C G;
# do
#     echo "current homo$HOMO"
#     mkdir -p $OUT/nuc_$lower_limit/homo$HOMO/
#     for file in $BAM/*.sorted.bam;
#     do
#     sample=$(basename $file ".sorted.bam")
#     echo $sample

#     if [[ $sample == Copen* ]]
#     then
#         REF="GCF_000007685.1_ASM768v1_genomic"
#         REF2="GCF_000007685.1_ASM768v1"
#     elif [[ $sample == Man* ]] || [[ $sample == Icter* ]]
#     then
#         REF="GCF_014858895.1_ASM1485889v1_genomic"
#         REF2="GCF_014858895.1_ASM1485889v1"
#     elif [[ $sample == Pa* ]] 
#     then
#         REF="GCF_000017685.1_ASM1768v1_genomic"
#         REF2="GCF_000017685.1_ASM1768v1"
#     else
#         REF=""
#         REF2=""
#     fi

#     echo $REF

#     echo "      - get overall coverage within CDS region"
#     # # get overall coverage for each transcript
#     bedtools coverage -a $MAIN/reference/$REF.CDS.bed -b $file -d > $OUT/overall_trans_cov/$sample.cov

#     echo "      - get max coverage within each CDS region"
#     # # get max cov from each transcript
#     python trans_max_cov_pos.py $OUT/overall_trans_cov/$sample.cov $OUT/overall_trans_cov/$sample.max.cov

#     echo "      - only get homo region identified with CDS region"
#     cat $MAIN/bed/${REF2}_homo${HOMO}_$lower_limit.anot.bed | grep -v "noncoding" > $MAIN/bed/${REF2}_homo${HOMO}_$lower_limit.anot.CDS.bed

#     echo "      - get coverage 100bps within CDS homo region"
#     # get coverage around each homopolymer region
#     bedtools coverage -a $MAIN/bed/${REF2}_homo${HOMO}_$lower_limit.anot.CDS.bed -b $file -d > $OUT/nuc_$lower_limit/homo$HOMO/$sample.$HOMO$lower_limit.CDS.cov

#     echo "      - get relative coverage 100bps within CDS homo region (based on maximum coverage of each CDS region)"
#     # # get relative cov within the flanking regions of homopolymer by divide each base' cov with the max cov of the transcript
#     python relative_cov_per_base.py $OUT/overall_trans_cov/$sample.max.cov $OUT/nuc_$lower_limit/homo$HOMO/$sample.$HOMO$lower_limit.CDS.cov $OUT/nuc_$lower_limit/homo$HOMO/$sample.$HOMO$lower_limit.CDS.relative.cov True

#     done
# done

# ######## cDNA single strand cov #################################
# ml Anaconda3/2022.05
# ml BEDTools/2.30.0-GCC-8.3.0

# BAM="/scratch/rx32940/minION/polyA_cDNA/map/genome/bam" #if map_full also need to change search pattern to *.withClip.sorted.bam 
# OUT="/scratch/rx32940/minION/homopolymer_cov_genome/output"
# lower_limit=5 # at least how many nucleotide to count as homopolymer

# mkdir -p $OUT/overall_trans_cov

# for HOMO in A T C G;
# do
#     echo "current homo$HOMO"
#     mkdir -p $OUT/nuc_$lower_limit/homo$HOMO/
#     for file in $BAM/*.sorted.bam;
#     do
#     sample=$(basename $file ".sorted.bam")
#     echo $sample

#     if [[ $sample == Copen* ]]
#     then
#         REF="GCF_000007685.1_ASM768v1_genomic"
#         REF2="GCF_000007685.1_ASM768v1"
#     elif [[ $sample == Man* ]] || [[ $sample == Icter* ]]
#     then
#         REF="GCF_014858895.1_ASM1485889v1_genomic"
#         REF2="GCF_014858895.1_ASM1485889v1"
#     elif [[ $sample == Pa* ]] 
#     then
#         REF="GCF_000017685.1_ASM1768v1_genomic"
#         REF2="GCF_000017685.1_ASM1768v1"
#     else
#         REF=""
#         REF2=""
#     fi

#     echo $REF

#     echo "      - get overall coverage within CDS region"
#     # # get overall coverage for each transcript
#     bedtools coverage -s -a $MAIN/reference/$REF.CDS.bed -b $file -d > $OUT/overall_trans_cov/$sample.ss.cov

#     echo "      - get max coverage within each CDS region"
#     # # get max cov from each transcript
#     python trans_max_cov_pos.py $OUT/overall_trans_cov/$sample.ss.cov $OUT/overall_trans_cov/$sample.ss.max.cov

#     # echo "      - only get homo region identified with CDS region"
#     # cat $MAIN/bed/${REF2}_homo${HOMO}_$lower_limit.anot.bed | grep -v "noncoding" > $MAIN/bed/${REF2}_homo${HOMO}_$lower_limit.anot.CDS.bed

#     echo "      - get coverage 100bps within CDS homo region"
#     # get coverage around each homopolymer region
#     bedtools coverage -s -a $MAIN/bed/${REF2}_homo${HOMO}_$lower_limit.anot.CDS.bed -b $file -d > $OUT/nuc_$lower_limit/homo$HOMO/$sample.$HOMO$lower_limit.ss.CDS.cov

#     echo "      - get relative coverage 100bps within CDS homo region (based on maximum coverage of each CDS region)"
#     # # get relative cov within the flanking regions of homopolymer by divide each base' cov with the max cov of the transcript
#     python relative_cov_per_base.py $OUT/overall_trans_cov/$sample.ss.max.cov $OUT/nuc_$lower_limit/homo$HOMO/$sample.$HOMO$lower_limit.ss.CDS.cov $OUT/nuc_$lower_limit/homo$HOMO/$sample.$HOMO$lower_limit.ss.CDS.relative.cov True

#     done
# done


######## cDNA FL #################################
# ml Anaconda3/2022.05
# ml BEDTools/2.30.0-GCC-8.3.0

# BAM="/scratch/rx32940/minION/polyA_cDNA/map_full/genome/bam" #if map_full also need to change search pattern to *.withClip.sorted.bam 
# OUT="/scratch/rx32940/minION/homopolymer_cov_genome/output"
# lower_limit=5 # at least how many nucleotide to count as homopolymer


# mkdir -p $OUT/overall_trans_cov

# for HOMO in A T C G;
# do
#     echo "current homo$HOMO"
#     mkdir -p $OUT/nuc_$lower_limit/homo$HOMO/
#     for file in $BAM/*.withClip.sorted.bam;
#     do
#     sample=$(basename $file ".withClip.sorted.bam")
#     echo $sample

#     if [[ $sample == Copen* ]]
#     then
#         REF="GCF_000007685.1_ASM768v1_genomic"
#         REF2="GCF_000007685.1_ASM768v1"
#     elif [[ $sample == Man* ]] || [[ $sample == Icter* ]]
#     then
#         REF="GCF_014858895.1_ASM1485889v1_genomic"
#         REF2="GCF_014858895.1_ASM1485889v1"
#     elif [[ $sample == Pa* ]] 
#     then
#         REF="GCF_000017685.1_ASM1768v1_genomic"
#         REF2="GCF_000017685.1_ASM1768v1"
#     else
#         REF=""
#         REF2=""
#     fi

#     echo $REF

#     echo "      - get overall coverage within CDS region"
#     # # get overall coverage for each transcript
#     bedtools coverage -s -a $MAIN/reference/$REF.CDS.bed -b $file -d > $OUT/overall_trans_cov/$sample.FL.cov

#     echo "      - get max coverage within each CDS region"
#     # # get max cov from each transcript
#     python trans_max_cov_pos.py $OUT/overall_trans_cov/$sample.FL.cov $OUT/overall_trans_cov/$sample.FL.max.cov

#     echo "      - only get homo region identified with CDS region"
#     cat $MAIN/bed/${REF2}_homo${HOMO}_$lower_limit.anot.bed | grep -v "noncoding" > $MAIN/bed/${REF2}_homo${HOMO}_$lower_limit.anot.CDS.bed

#     echo "      - get coverage 100bps within CDS homo region"
#     # get coverage around each homopolymer region
#     bedtools coverage -s -a $MAIN/bed/${REF2}_homo${HOMO}_$lower_limit.anot.CDS.bed -b $file -d > $OUT/nuc_$lower_limit/homo$HOMO/$sample.$HOMO$lower_limit.FL.CDS.cov

#     echo "      - get relative coverage 100bps within CDS homo region (based on maximum coverage of each CDS region)"
#     # # get relative cov within the flanking regions of homopolymer by divide each base' cov with the max cov of the transcript
#     python relative_cov_per_base.py $OUT/overall_trans_cov/$sample.FL.max.cov $OUT/nuc_$lower_limit/homo$HOMO/$sample.$HOMO$lower_limit.FL.CDS.cov $OUT/nuc_$lower_limit/homo$HOMO/$sample.$HOMO$lower_limit.FL.CDS.relative.cov True

#     done
# done

##########################################################
#
# plot relative expression
#
#########################################################

ml R/4.2.1-foss-2020b
ml LibTIFF/4.1.0-GCCcore-10.2.0
main_path="/scratch/rx32940/minION/homopolymer_cov_genome/output/nuc_5"
WORK="/scratch/rx32940/minION/homopolymer_cov_genome"

for file in $main_path/homo*/;
do
nuc=$(basename $file)

echo $nuc
Rscript $WORK/homo_plot.r $main_path $nuc

done