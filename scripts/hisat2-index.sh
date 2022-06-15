#!/bin/bash

#SBATCH --job-name=hisat-index
#SBATCH --output=/home/lspencer/2022-redking-OA/sbatch_logs/hisat_index.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 21-0:0:0

# Load modules
source /home/lspencer/venv/bin/activate
module load bio/hisat2/2.2.1

ref_dir="/home/lspencer/references/redkingcrab"

# Generate HiSat2 index
hisat2-build -p 20 ${ref_dir}/Paralithodes.camtschaticus.genome.fasta \
${ref_dir}/hisat2/Paralithodes.camtschaticus
