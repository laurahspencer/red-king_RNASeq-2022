#!/bin/bash

#SBATCH --job-name=redking_sort-dedupe
#SBATCH --output=/home/lspencer/2022-redking-OA/sbatch_logs/redking_sort-dedupe-v2.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 10
#SBATCH -t 5-0:0:0

# This script is for sorting then deduplicating aligned RNASeq data

#source /home/lspencer/venv/bin/activate
module load bio/samtools/1.11
module load bio/gatk/4.2.0.0
conda activate gatk

IN=/scratch/lspencer/2022-redking-OA/aligned-bowtie
OUT=/scratch/lspencer/2022-redking-OA/aligned-bowtie/sorted-deduped

#Convert .sam files to .bam files & sort by coordinate
for file in ${IN}/*.bowtie.v2.sam
do
sample="$(basename -a $file | cut -d "." -f 1)"
sorted_file="$sample.sorted.bam"
samtools view -@10 -b $file | samtools sort -@10 -o ${OUT}/$sorted_file
done


# Deduplicate
for file in ${OUT}/*.sorted.bam
do
sample="$(basename -a $file | cut -d "." -f 1)"
dedup_file="$sample.dedup.bam"

# deduplicate using picard, output will have duplicates removed
gatk MarkDuplicates \
I=$file \
O="${OUT}/$dedup_file" \
M="${OUT}/$sample.dup_metrics.txt" \
REMOVE_DUPLICATES=true
done >> "${OUT}/dedup_stout.txt" 2>&1
