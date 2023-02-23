#!/bin/bash

#SBATCH --job-name=feature-counts
#SBATCH --output=/home/lspencer/2022-redking-OA/sbatch_logs/redking_feature-counts_RKC-genome.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 5-0:0:0

module load bio/subread/2.0.3
module load bio/samtools/1.11

IN=/scratch/lspencer/2022-redking-OA/aligned/bowtie-rkc
OUT=/home/lspencer/2022-redking-OA/aligned/bowtie-rkc
GFF=/home/lspencer/references/redkingcrab

#Convert .sam files to .bam files & sort by coordinate
for file in ${IN}/*.bowtie.v2.sam
do
sample="$(basename -a $file | cut -d "." -f 1)"
sorted_file="$sample.sorted.bam"
samtools view -@10 -b $file | samtools sort -@10 -o ${OUT}/$sorted_file
done

#v3
### # Summarize paired-end reads and count fragments (instead of reads), do NOT count chimeric fragments (fragments mapping to different chromosomes)
featureCounts \
-p --countReadPairs \
-T 20 \
-C \
-t gene \
-g ID \
-a ${GFF}/Paralithodes.camtschaticus.genome.gff \
-o ${OUT}/featurecounts_rkc-genome_gene-v3 \
${OUT}/*.sorted.bam

#v2
# Summarize paired-end reads and count fragments (instead of reads) which have both ends successfully aligned, i.e. don't count singletons
featureCounts \
-p --countReadPairs \
-B \
-T 20 \
-C \
-t gene \
-g ID \
-a ${GFF}/Paralithodes.camtschaticus.genome.gff \
-o ${OUT}/featurecounts_rkc-genome_gene-v2 \
${OUT}/*.sorted.bam

# v7
# Summarize paired-end reads and count fragments (instead of reads) which have both ends successfully aligned, i.e. don't count singletons, and don't count chimeras
featureCounts \
-p --countReadPairs \
-B \
-T 20 \
-C \
-t gene \
-g ID \
-a ${GFF}/Paralithodes.camtschaticus.genome.gff \
-o ${OUT}/featurecounts_rkc-genome_gene-v7 \
${OUT}/*.sorted.bam
