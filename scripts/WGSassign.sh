#!/bin/bash

#SBATCH --job-name=WGSassign
#SBATCH --output=/home/lspencer/pcod-lcwgs-2023/analysis-20231002/sbatch-WGSassign.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 5-0:0:0

conda activate /home/lspencer/programs/WGSassign/WGSassign

IN=/home/lspencer/pcod-lcwgs-2023/analysis-20231002
OUT=/home/lspencer/pcod-lcwgs-2023/analysis-20231002/WGSassign

reference_beagle=${IN}/gls/FILE.beagle.gz
reference_IDs=${OUT}/reference-IDs.txt
outname=${OUT}/wgsassign-reference-AF.out

# Estimate reference population allele frequencies using 20 threads
# Output = ${outname}.pop_af.npy (numpy binary matrix of size L (# loci) rows x K (ref pops) columns)
WGSassign --beagle ${reference_beagle} --pop_af_IDs ${reference_IDs} --get_reference_af --out ${outname} --threads 20
