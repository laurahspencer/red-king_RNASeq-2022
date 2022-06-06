#!/bin/bash

#SBATCH --job-name=feature-counts
#SBATCH --output=/home/lspencer/2022-redking-OA/sbatch_logs/feature-counts_star.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 5-0:0:0

# This script is for summarizing the number of fragments (paired-end reads) align to each blue king crab gene using various settings

module load bio/subread/2.0.3

IN=/scratch/lspencer/2022-redking-OA/aligned/star
OUT=/home/lspencer/2022-redking-OA/counts
GFF=/home/lspencer/references/bluekingcrab

### Count fragments, don't count chimeras, count all multimapped
featureCounts \
-p --countReadPairs \
-T 20 \
-C \
-M \
-t gene \
-g ID \
-a ${GFF}/EVM.out_new.gff3 \
-o ${OUT}/featurecounts_redking_gene_star_multimapped-all \
${IN}/*.bam

### Count fragments, don't count chimeras, count multimapped fragments fractionally
featureCounts \
-p --countReadPairs \
-T 20 \
-C \
-M --fraction \
-t gene \
-g ID \
-a ${GFF}/EVM.out_new.gff3 \
-o ${OUT}/featurecounts_redking_gene_star_multimapped-fraction \
${IN}/*.bam

### Count fragments, don't count chimeras, don't count multimapped fragments
featureCounts \
-p --countReadPairs \
-T 20 \
-C \
-t gene \
-g ID \
-a ${GFF}/EVM.out_new.gff3 \
-o ${OUT}/featurecounts_redking_gene_star_multimapped-none \
${IN}/*.bam
