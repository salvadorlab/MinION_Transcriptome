#!/bin/sh
#SBATCH --partition=batch
#SBATCH --job-name=tama_go
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=100:00:00
#SBATCH --mem=50G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL

source activate tama
ml SAMtools/1.10-GCC-8.3.0
ml BEDTools/2.30.0-GCC-8.3.0

MAIN="/scratch/rx32940/minION/polyA_directRNA/TAMA"
TAMA="$MAIN/tama"
cdna_bam="/scratch/rx32940/minION/polyA_cDNA/map/genome/bam"
REF="/scratch/rx32940/minION/polyA_cDNA/map/genome/reference"
direct_bam="/scratch/rx32940/minION/polyA_directRNA/map/genome/bam"

#########################################################################################
######## extract transcript regions merged from different seq_method's samples from reference fasta
###########################################################################################
bedtools getfasta -name -split -s -fi $REF/GCF_000007685.1_ASM768v1_genomic.fna -bed $MAIN/merged/cDNA_all_samples.bed -fo $MAIN/merged/cDNA_all_samples.fasta

bedtools getfasta -name -split -s -fi $REF/GCF_000007685.1_ASM768v1_genomic.fna -bed $MAIN/merged/dRNA_all_samples.bed -fo $MAIN/merged/dRNA_all_samples.fasta

#########################################################################################
######## search for ORF from the transcripts fasta extracted
###########################################################################################

python $TAMA/tama_go/orf_nmd_predictions/tama_orf_seeker.py -f $MAIN/merged/cDNA_all_samples.fasta -o $MAIN/merged/cDNA_all_samples.faa

python $TAMA/tama_go/orf_nmd_predictions/tama_orf_seeker.py -f $MAIN/merged/dRNA_all_samples.fasta -o $MAIN/merged/dRNA_all_samples.faa


#########################################################################################
######## blastp; blast open reading frame against db
###########################################################################################
# wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz -P $MAIN/uniprot
ml BLAST+/2.11.0-gompi-2020b

# makeblastdb -dbtype prot -in $MAIN/uniprot/uniprot_sprot.fasta -input_type fasta -out $MAIN/uniprot/uniprot_sprot

blastp -evalue 1e-10 -num_threads 16 -db $MAIN/uniprot/uniprot_sprot -query $MAIN/merged/cDNA_all_samples.faa -out $MAIN/merged/cDNA_all_samples.blastp.txt

blastp -evalue 1e-10 -num_threads 16 -db $MAIN/uniprot/uniprot_sprot -query $MAIN/merged/dRNA_all_samples.faa -out $MAIN/merged/dRNA_all_samples.blastp.txt


#########################################################################################
######## parse blastp output from previous step
###########################################################################################


python $TAMA/tama_go/orf_nmd_predictions/tama_orf_blastp_parser.py -b $MAIN/merged/cDNA_all_samples.blastp.txt -o $MAIN/merged/cDNA_all_samples.blastp.parsed.txt

python $TAMA/tama_go/orf_nmd_predictions/tama_orf_blastp_parser.py -b $MAIN/merged/dRNA_all_samples.blastp.txt -o $MAIN/merged/dRNA_all_samples.blastp.parsed.txt

#########################################################################################
######## parse blastp to bed
###########################################################################################

python $TAMA/tama_go/orf_nmd_predictions/tama_cds_regions_bed_add.py -p $MAIN/merged/cDNA_all_samples.blastp.parsed.txt -a $MAIN/merged/cDNA_all_samples.bed -f $MAIN/merged/cDNA_all_samples.fasta -o $MAIN/merged/cDNA_all_samples_orf.bed

python $TAMA/tama_go/orf_nmd_predictions/tama_cds_regions_bed_add.py -p $MAIN/merged/dRNA_all_samples.blastp.parsed.txt -a $MAIN/merged/dRNA_all_samples.bed -f $MAIN/merged/dRNA_all_samples.fasta -o $MAIN/merged/dRNA_all_samples_orf.bed

cat $MAIN/merged/cDNA_all_samples_orf.bed | grep -v ";no_hit;" > $MAIN/merged/cDNA_all_samples_orf_hit.bed 
cat $MAIN/merged/dRNA_all_samples_orf.bed | grep -v ";no_hit;" > $MAIN/merged/dRNA_all_samples_orf_hit.bed