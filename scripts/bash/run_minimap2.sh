#!/bin/bash
#!/bin/sh
#SBATCH --partition=batch
#SBATCH --job-name=map_raw
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --time=100:00:00
#SBATCH --mem=10G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL

# OUT="/scratch/rx32940/minION/polyA_cDNA/map/genome"
# FASTQ="/scratch/rx32940/minION/polyA_cDNA/fastq"
# ref="$OUT/reference"
# for file in $FASTQ/*.fastq;
# do
# (
#     sample=$(basename $file ".fastq")
#     header="
#     #!/bin/sh\n#SBATCH --partition=batch\n#SBATCH --job-name=minimap2_$sample\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=3\n#SBATCH --time=100:00:00\n#SBATCH --mem=10G\n#SBATCH --output=%x.%j.out\n#SBATCH --error=%x.%j.out\n#SBATCH --mail-user=rx32940@uga.edu\n#SBATCH --mail-type=ALL\n 

#     ml minimap2/2.17-GCC-8.3.0\n
#     ml SAMtools/1.10-GCC-8.3.0\n
#     ml BEDTools/2.30.0-GCC-8.3.0\n
#     mkdir -p $OUT/bam $OUT/sam $OUT/bed\n
#     source activate pychopper\n
#     "

#     if [[ $sample == Copen* ]]
#     then
#         REF="GCF_000007685.1_ASM768v1_genomic";
#     elif [[ $sample == Man* ]] || [[ $sample == Icter* ]]
#     then
#         REF="GCF_014858895.1_ASM1485889v1_genomic"
#     elif [[ $sample == Pa* ]] 
#     then
#         REF="GCF_000017685.1_ASM1768v1_genomic"
#     else
#         REF=""
#     fi

#     echo -e $header > sub.sh
#     echo -e "
#     # minimap2 -ax splice -k14 -p 0.99 --MD $ref/$REF.fna $file > $OUT/sam/$sample.sam\n
#     samtools view -bS $OUT/sam/$sample.sam > $OUT/bam/$sample.bam \n
#     samtools view -@ 3 -bS $OUT/sam/$sample.sam | samtools sort - -@ 3 -o $OUT/bam/$sample.sorted.bam\n 
#     samtools index $OUT/bam/$sample.sorted.bam \n
#     bedtools bamtobed -i $OUT/bam/$sample.sorted.bam -cigar -tag NM > $OUT/bed/$sample.bed\n
#     " >> sub.sh

#     sbatch sub.sh
# ) &

# wait
# done

cd $SLURM_SUBMIT_DIR

ml SAMtools/1.10-GCC-8.3.0
ml Qualimap/2.2.1-foss-2019b-R-3.6.2
mkdir -p sam_stats qualimap sam_flagstat

for file in /scratch/rx32940/minION/polyA_cDNA/map/genome/bam/*sorted.bam;
do 
sample=$(basename $file ".sorted.bam")
echo $sample
samtools stats -d -@ 3 sam/$sample.sam > sam_stats/$sample.sam.stats
samtools flagstats sam/$sample.sam  > sam_flagstat/$sample.flagstat
mkdir -p qualimap/$sample
qualimap bamqc -bam bam/$sample.sorted.bam -nw 5000 -nt 3 -c -outdir qualimap/$sample -gff /scratch/rx32940/minION/map/reference/GCF_000007685.1_ASM768v1_genomic.gtf
done

