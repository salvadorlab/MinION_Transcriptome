#!/bin/sh
#SBATCH --partition=gpu_p
#SBATCH --job-name=demultiplex
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=100:00:00
#SBATCH --mem=10G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out
#SBATCH --mail-user=rx32940@uga.edu
#SBATCH --mail-type=ALL

cd $SLURM_SUBMIT_DIR

ml ont-guppy/4.4.2-GPU

# #############################################################################
#
# rebasecalling raw fast5 files for direct RNA samples to perform tailfindr
#    - direct RNA no need to demultiplex, already in separete folder
#
# ##################################################################################

# input="/scratch/rx32940/minION/MinION.Fast5.Data/Direct_RNA_FAST5_Pass"
# output="/scratch/rx32940/minION/Guppy/Direct_RNA_FAST5_Pass"

# guppy_basecaller -x "cuda:0" --compress_fastq \
# --input_path $input/Q36_Copenhageni_Direct-RNA_Repeat_FAST5_Pass \
# --save_path $output/Q36_Copenhageni_Direct-RNA_Repeat_FAST5_Pass \
# --config rna_r9.4.1_70bps_hac.cfg \
# --reverse_sequence true --u_substitution true \
# --fast5_out --num_callers 1 --gpu_runners_per_device 8

# #############################################################################
#
# demultiplex and rebasecalling raw fast5 files for cDNA samples to perform tailfindr
#    
# ##################################################################################

#### STEP 1: basecalling ##########

input="/scratch/rx32940/minION/Guppy/demultiplex"
output="/scratch/rx32940/minION/Guppy/demultiplex/rebasecall"

for dir in $input/output/*;
do

sample=$(basename $dir)

mkdir -p $output/$sample

guppy_basecaller -x "cuda:0"  --compress_fastq --input_path $dir --save_path $output/$sample --config dna_r9.4.1_450bps_hac.cfg --fast5_out –cpu_threads_per_caller 5 --num_callers 1

done
#### STEP 2: demultiplex ##########

# input="/scratch/rx32940/minION/Guppy/PolyA_cDNA_FAST5_Pass"
# output="/scratch/rx32940/minION/Guppy/PolyA_cDNA_demultiplex_FAST5_Pass"

# guppy_barcoder --input_path $input --save_path $output --barcode_kits EXP-NBD104 --num_barcoding_buffers 8

### STEP 3: merge each sample's fastq files together

# # base on the barcode to sample information provided by Dr. Rajeev, combined fastq is also renamed
# input="/scratch/rx32940/minION/Guppy/PolyA_cDNA_demultiplex_FAST5_Pass"

# for file in $input/*/; # only iterate dirs
# do
# barcode=$(basename $file)
# echo $barcode
# cat $file/* > $file/$barcode.fastq
# done

