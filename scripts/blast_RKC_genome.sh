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
module load bio/bedtools/2.29.2

BASE=/home/lspencer/references

# Extract sequences for just genes in fasta format
bedtools getfasta \
-fi ${BASE}/redkingcrab/Paralithodes.camtschaticus.genome.fasta \
-bed ${BASE}/redkingcrab/Paralithodes.camtschaticus_genes.bed \
-fo ${BASE}/redkingcrab/Paralithodes.camtschaticus_genes.fasta

# Blast genes against uniprot/swissprot
blastx \
-query ${BASE}/redkingcrab/Paralithodes.camtschaticus_genes.fasta \
-db ${BASE}/blast/uniprot_sprot_20220111_protein \
-out ${BASE}/redkingcrab/Paralithodes.camtschaticus_genes_blastx.tab \
-evalue 1E-5 \
-num_threads 20 \
-outfmt 6
