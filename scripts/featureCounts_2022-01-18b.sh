#!/bin/bash

#SBATCH --job-name=redking_feature-counts
#SBATCH --output=/home/lspencer/2022-redking-OA/sbatch_logs/redking_feature-counts_gene-2022-01-18b.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 5-0:0:0

# This script is a continuation of another script dated same day, since that job was interrupted.

module load bio/subread/2.0.3

IN=/scratch/lspencer/2022-redking-OA/aligned/bowtie2-sorted
OUT=/home/lspencer/2022-redking-OA/counts
GFF=/home/lspencer/references/bluekingcrab

### Version 6 settings, a combination of settings that seem to be a good idea:
# - Summarize paired-end reads and count fragments (instead of reads) which have both ends successfully aligned
# - Increase minimum bp overlap to 2 (summed across both ends for paired-end data)
# - Do NOT count chimeric fragments (those mapping to multiple chromosomes)
# - Count multi-mapping reads fractionally

featureCounts \
-p --countReadPairs \
-B \
--minOverlap 2 \
-C \
-M --fraction \
-T 20 \
-t gene \
-g ID \
-a ${GFF}/EVM.out_new.gff3 \
-o ${OUT}/featurecounts_redking_gene-v6 \
${IN}/*.sorted.bam
