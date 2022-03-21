#!/bin/bash

THREAD=24

input="/scratch/rx32940/minION/Guppy/Direct_RNA_FAST5_Pass"
output="/scratch/rx32940/minION/polyA_directRNA/tailfindr/tailfindr_out/"

header="
#!/bin/sh\n\
#SBATCH --partition=batch\n\
#SBATCH --job-name=tailfindr\n\
#SBATCH --ntasks=1\n\
#SBATCH --cpus-per-task=$THREAD\n\
#SBATCH --time=168:00:00\n\
#SBATCH --mem-per-cpu=10G\n\
#SBATCH --output=%x.%j.out\n\
#SBATCH --error=%x.%j.out\n\
#SBATCH --mail-user=rx32940@uga.edu\n\
#SBATCH --mail-type=ALL\n 
"

for dir in $input/*;
do

(

sample_name=$(basename "$dir")
echo $sample_name

echo -e $header > sub.sh

echo -e "source activate tailfindr" >> sub.sh

echo -e "Rscript run_tailfindr.r $dir $output $sample_name.csv $THREAD" >> sub.sh

echo -e "conda deactivate" >> sub.sh

sbatch sub.sh

) &

wait
done
