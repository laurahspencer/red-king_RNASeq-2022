#!/bin/bash

#SBATCH --job-name=redking_align-v2
#SBATCH --output=/home/lspencer/2022-redking-OA/sbatch_logs/redking_align-v2.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 10-0:0:0

# This script is for aligning trimmed RNA-Seq data to the blue king crab genome
# This version 2 uses the Dryad version of the blue king crab genome

#source /home/lspencer/venv/bin/activate
module load aligners/bowtie2/2.4.2

REF=/home/lspencer/references/bluekingcrab
IN=/home/lspencer/2022-redking-OA/trimmed
OUT=/scratch/lspencer/2022-redking-OA/aligned-bowtie
VER=2

# Move to the directory containing blue king crab genome.
# Bowtie2 looks for index in current directory - easier just to be in it.
cd ${REF}

### Create bowtie2 index for blue king crab genome
bowtie2-build \
--large-index \
--threads 20 \
-f Paralithodes_platypus_genome.fasta.gz \
Paralithodes_platypus_dryad

# Run Bowtie2 over each RNASeq paired sample
for file in ${IN}/*.trimmed.R1.v2.fastq.gz
do
sample="$(basename -a ${file} | cut -d "." -f 1)"
file_R1="${sample}.trimmed.R1.v2.fastq.gz"
file_R2="${sample}.trimmed.R2.v2.fastq.gz"
map_file="${sample}.bowtie.v${VER}.sam"

# run Bowtie2 on each file
bowtie2 \
-x Paralithodes_platypus_dryad \
--sensitive \
--threads 20 \
--no-unal \
-1 ${IN}/${file_R1} \
-2 ${IN}/${file_R2} \
-S ${OUT}/${map_file}; \
done >> ${OUT}/bowtieout.v${VER}.txt 2>&1
