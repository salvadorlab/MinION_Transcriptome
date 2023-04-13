#!/bin/bash

ml SAMtools/1.14-GCC-8.3.0
BAM="/scratch/rx32940/minION/polyA_cDNA/map/genome/bam"
WORK="/scratch/rx32940/minION/polyA_cDNA/ladder_operon"

samtools view -N $WORK/flagellar_operon_reads_LIC-cDNA-polyA.txt -o $WORK/flagellar_operon_reads_LIC-cDNA-polyA.bam $BAM/Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_PolyATail.sorted.bam
samtools index $WORK/flagellar_operon_reads_LIC-cDNA-polyA.bam