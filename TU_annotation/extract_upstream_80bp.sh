#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=motif_finding
#SBATCH --ntasks=1                    	
#SBATCH --cpus-per-task=1            
#SBATCH --time=100:00:00
#SBATCH --mem=10G
#SBATCH --output=%x.%j.out       
#SBATCH --error=%x.%j.out        
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL 

BED="/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/direct_output/tss/upstream80"
REF="/scratch/rx32940/minION/polyA_cDNA/map/genome/reference"
# ml BEDTools/2.30.0-GCC-8.3.0
# ml GenomeTools/1.6.1-GCC-8.3.0-Python-2.7.16
# ml seqkit/0.16.1

# #### extract sequences 80 base pairs upstream of pTSS identified
# bedtools getfasta -bed $BED/Q29_TSS_upstream_80bp.bed -fi $REF/GCF_000007685.1_ASM768v1_genomic.fna -fo $BED/Q29_TSS_upstream_80bp.fasta
# bedtools getfasta -bed $BED/Q36_TSS_upstream_80bp.bed -fi $REF/GCF_000007685.1_ASM768v1_genomic.fna -fo $BED/Q36_TSS_upstream_80bp.fasta

# ### only keep the unique sequences in the upstream fasta
# gt sequniq $BED/Q29_TSS_upstream_80bp.fasta > $BED/Q29_TSS_upstream_80bp.uniq.fasta
# gt sequniq $BED/Q36_TSS_upstream_80bp.fasta > $BED/Q36_TSS_upstream_80bp.uniq.fasta

# ### remove fasta sequences shorter than 8bp (won't take by MEME)
# seqkit seq -m 8 $BED/Q29_TSS_upstream_80bp.uniq.fasta > $BED/Q29_TSS_upstream_80bp.uniq.fil8.fasta
# seqkit seq -m 8 $BED/Q36_TSS_upstream_80bp.uniq.fasta > $BED/Q36_TSS_upstream_80bp.uniq.fil8.fasta

### MOTIF finding
ml MEME/5.4.1-foss-2019b-Python-3.7.4

# xstreme will use MEME and streme algorithm to identify motif and do motif enrichment analysis.
meme $BED/Q29_TSS_upstream_80bp.uniq.fil8.fasta -dna -oc . -nostatus -time 14400 -mod zoops -nmotifs 5 -minw 5 -maxw 50 -objfun classic -revcomp -markov_order 0

treme --verbosity 1 --oc . --dna --totallength 4000000 --time 14400 --minw 5 --maxw 15 --thresh 0.05 --align right --p Q29_TSS_upstream_80bp.uniq.fil8.fasta
# xstreme -o $BED/Q29_streme --align right --desc --dna --p $BED/Q29_TSS_upstream_80bp.uniq.fil8.fasta --fimo-skip 
# xstreme -o $BED/Q36_streme --align right --desc --dna --p $BED/Q36_TSS_upstream_80bp.uniq.fil8.fasta --fimo-skip 


