#!/bin/bash

#SBATCH --job-name=wgcna
#SBATCH --output=/home/lspencer/2022-redking-OA/sbatch_logs/wgcna.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH --mem=180GB
#SBATCH -t 7-0:0:0

module load R
R CMD BATCH /home/lspencer/2022-redking-OA/scripts/wgcna.R /home/lspencer/2022-redking-OA/wgcna/wgcna-output.txt
