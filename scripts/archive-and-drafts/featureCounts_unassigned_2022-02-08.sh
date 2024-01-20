#!/bin/bash

#SBATCH --job-name=redking_feature-counts
#SBATCH --output=/home/lspencer/2022-redking-OA/sbatch_logs/redking_feature-counts_unassigned.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 5-0:0:0

# This script is for summarizing the number of fragments (paired-end reads) align to each blue king crab gene using various settings

module load bio/subread/2.0.3

IN=/home/lspencer/2022-redking-OA/aligned/bowtie2-sorted
OUT=/home/lspencer/2022-redking-OA/counts
GFF=/home/lspencer/references/bluekingcrab

# Summarize paired-end reads and count fragments (instead of reads)
# Include singletons
# Don't include chimeric reads (those mapping to diff. chromosomes)
# Output detailed read assignment results for each fragment
featureCounts \
-p --countReadPairs \
-T 20 \
-t gene \
-g ID \
-C \
-R BAM \
-a ${GFF}/EVM.out_new.gff3 \
-o ${OUT}/featurecounts_redking_gene-v7 \
${IN}/*.sorted.bam
