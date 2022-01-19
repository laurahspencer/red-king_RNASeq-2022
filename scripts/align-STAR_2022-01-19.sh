#!/bin/bash

#SBATCH --job-name=redking_STAR
#SBATCH --output=/home/lspencer/2022-redking-OA/sbatch_logs/redking_STAR.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 24
#SBATCH -t 10-0:0:0

BASE=/scratch/lspencer/2022-redking-OA
REF=/home/lspencer/references/bluekingcrab
IN=${BASE}/trimmed
OUT=${BASE}/aligned/STAR

SAMPLES=$(ls ${IN}/*.trimmed.R1.v2.fastq.gz | \
	awk -F "/" '{print $NF}' | \
	awk -F "." '{print $1}')
