#!/bin/bash

THREAD=12

input="/scratch/rx32940/minION/Guppy/fast5"
output="/scratch/rx32940/minION/polyA_cDNA/tailfindr/tailfindr_2021_basecall"



for dir in $input/*PolyATail;
do

(

sample_name=$(basename "$dir")
echo $sample_name

header="
#!/bin/sh\n\
#SBATCH --partition=batch\n\
#SBATCH --job-name=${sample_name}_tailfindr\n\
#SBATCH --ntasks=1\n\
#SBATCH --cpus-per-task=$THREAD\n\
#SBATCH --time=168:00:00\n\
#SBATCH --mem-per-cpu=10G\n\
#SBATCH --output=%x.%j.out\n\
#SBATCH --error=%x.%j.out\n\
#SBATCH --mail-user=rx32940@uga.edu\n\
#SBATCH --mail-type=ALL\n 
"

echo -e $header > sub.sh

echo -e "source activate tailfindr" >> sub.sh

echo -e "Rscript run_tailfindr.r $dir $output $sample_name.csv $THREAD" >> sub.sh

echo -e "conda deactivate" >> sub.sh

sbatch sub.sh

) &

wait
done
