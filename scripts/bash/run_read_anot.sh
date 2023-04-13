#!/bin/bash
#!/bin/sh
#SBATCH --partition=batch
#SBATCH --job-name=anot_reads
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --time=100:00:00
#SBATCH --mem=10G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL

ml tqdm/4.62.2-GCCcore-8.3.0-Python-3.8.2
MAP="/scratch/rx32940/minION/polyA_cDNA/map/genome"
ref="$MAP/reference"
BED="$MAP/bed"
for file in $BED/*PolyATail.bed;
do
sample=$(basename $file '.bed')
python $MAP/calculate_aln_identity.py $file $MAP/stats/$sample.stats 
echo $sample
if [[ $sample == Copen* ]]
then
    REF="GCF_000007685.1_ASM768v1_genomic";
elif [[ $sample == Man* ]] || [[ $sample == Icter* ]]
then
    REF="GCF_014858895.1_ASM1485889v1_genomic"
elif [[ $sample == Pa* ]] 
then
    REF="GCF_000017685.1_ASM1768v1_genomic"
else
    REF=""
fi
echo $REF
python $MAP/anot_read_transcript_v2.py $MAP/stats/$sample.stats $MAP/stats/$sample.anot.raw.stats $ref/$REF.gff
rm $MAP/stats/$sample.stats
done