#!/bin/bash

#SBATCH --job-name=blast
#SBATCH --output=/home/lspencer/2022-redking-OA/sbatch_logs/blast_BKC_genefasta.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 21-0:0:0

# This script is for blasting the blue king crab coding sequences against
# the Uniqprot/Swissprot database

module load bio/blast/2.11.0+

BASE=/home/lspencer/references

blastx \
-query ${BASE}/bluekingcrab/Paralithodes_platypus_genome.fasta \
-db ${BASE}/blast/uniprot_sprot_20220111_protein \
-out ${BASE}/bluekingcrab/P.platypus.genome_blastx.tab \
-evalue 1E-5 \
-num_threads 20 \
-outfmt 6
