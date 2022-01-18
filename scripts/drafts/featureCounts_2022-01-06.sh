#!/bin/bash

#SBATCH --job-name=redking_feature-counts
#SBATCH --output=/home/lspencer/2022-redking-OA/sbatch_logs/redking_feature-counts_gene-2022-01-07.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 10
#SBATCH -t 5-0:0:0

# This script is for sorting then deduplicating aligned RNASeq data

module load bio/subread/2.0.3

IN=/scratch/lspencer/2022-redking-OA/aligned-bowtie/sorted-deduped
OUT=/home/lspencer/2022-redking-OA/counts
GFF=/home/lspencer/references/bluekingcrab

# Summarize paired-end reads and count fragments:
featureCounts \
-p --countReadPairs \
-T 10 \
-t gene \
-g ID \
-a ${GFF}/EVM.out_new.gff3 \
-o ${OUT}/featurecounts_redking_gene-v1 \
${IN}/*.sorted.bam


# use featureCounts to count the number of fragments (paired-end) which have both ends successfully aligned
featureCounts \
-p --countReadPairs \
-B \
-T 10 \
-t gene \
-g ID \
-a ${GFF}/EVM.out_new.gff3 \
-o ${OUT}/featurecounts_redking_gene-v2 \
${IN}/*.sorted.bam
