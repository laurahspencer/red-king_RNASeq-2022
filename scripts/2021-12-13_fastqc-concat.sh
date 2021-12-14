#!/bin/bash

#SBATCH --job-name=redking_raw-fastqc
#SBATCH --output=/home/lspencer/2022-redking-OA/sbatch_logs/redking_raw-fastqc.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 8
#SBATCH -D /share/nwfsc/ggoetz/red_king_crab/illumina/

# This script is for running FastQC and MultiQC on raw (but concatenated)
# red king crab RNA-Seq data

# Load modules and virtual environments
module load bio/fastqc
source /home/lspencer/venv/bin/activate

# run fastqc on each raw read file
fastqc \
--threads 8 \
*.fastq.gz \
--outdir /home/lspencer/2022-redking-OA/fastqc/concat/

# Run multiqc to summarize fastqc reports
multiqc \
/home/lspencer/2022-redking-OA/fastqc/concat/ \
--outdir /home/lspencer/2022-redking-OA/fastqc/concat/
