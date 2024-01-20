#!/bin/bash

#SBATCH --job-name=blast
#SBATCH --output=/home/lspencer/2022-redking-OA/sbatch_logs/blast_BKC.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 10
#SBATCH -t 1000

# This script is for blasting the blue king crab coding sequences against
# the Uniqprot/Swissprot database

module load bio/blast/2.11.0+

BASE=/home/lspencer/references

blastx \
-query ${BASE}/bluekingcrab/EVM.out_new.cds \
-db ${BASE}/blast/uniprot_sprot_20220111_protein \
-out ${BASE}/blast/EVM.out_new_blastx.tab \
-evalue 1E-20 \
-num_threads 10 \
-outfmt 6
