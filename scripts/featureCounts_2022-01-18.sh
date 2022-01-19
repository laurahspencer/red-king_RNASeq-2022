#!/bin/bash

#SBATCH --job-name=redking_feature-counts
#SBATCH --output=/home/lspencer/2022-redking-OA/sbatch_logs/redking_feature-counts_gene-2022-01-18.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 5-0:0:0

# This script is for summarizing the number of fragments (paired-end reads) align to each blue king crab gene using various settings

module load bio/subread/2.0.3

IN=/scratch/lspencer/2022-redking-OA/aligned-bowtie/sorted-deduped
OUT=/home/lspencer/2022-redking-OA/counts
GFF=/home/lspencer/references/bluekingcrab

# Summarize paired-end reads and count fragments (instead of reads)
featureCounts \
-p --countReadPairs \
-T 20 \
-t gene \
-g ID \
-a ${GFF}/EVM.out_new.gff3 \
-o ${OUT}/featurecounts_redking_gene-v1 \
${IN}/*.sorted.bam


# Summarize paired-end reads and count fragments (instead of reads) which have both ends successfully aligned
featureCounts \
-p --countReadPairs \
-B \
-T 20 \
-t gene \
-g ID \
-a ${GFF}/EVM.out_new.gff3 \
-o ${OUT}/featurecounts_redking_gene-v2 \
${IN}/*.sorted.bam


### Version 3 settings: same as version 1, but do NOT count chimeric fragments (fragments mapping to different chromosomes)
featureCounts \
-p --countReadPairs \
-T 20 \
-C \
-t gene \
-g ID \
-a ${GFF}/EVM.out_new.gff3 \
-o ${OUT}/featurecounts_redking_gene-v3 \
${IN}/*.sorted.bam


### Version 4 settings: same as version 1, but count multi-mapping reads fractionally
featureCounts \
-p --countReadPairs \
-T 20 \
-M --fraction \
-t gene \
-g ID \
-a ${GFF}/EVM.out_new.gff3 \
-o ${OUT}/featurecounts_redking_gene-v4 \
${IN}/*.sorted.bam


### Version 5 settings: same as version 2 (which requires both ends to be aligned), but increase minimum bp overlap to 2.
featureCounts \
-p --countReadPairs \
-B \
--minOverlap 2 \
-T 20 \
-t gene \
-g ID \
-a ${GFF}/EVM.out_new.gff3 \
-o ${OUT}/featurecounts_redking_gene-v5 \
${IN}/*.sorted.bam
