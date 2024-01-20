#!/bin/bash

#SBATCH --job-name=redking_align-STAR
#SBATCH --output=/home/lspencer/2022-redking-OA/sbatch_logs/redking_STAR.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 10-0:0:0

module load aligners/star/2.7.10a
module load bio/rsem/1.3.3

BASE=/scratch/lspencer/2022-redking-OA
REF=/home/lspencer/references/bluekingcrab
IN=${BASE}/trimmed
OUT=${BASE}/aligned/star-rsem

# ========= Build STAR genome index
# Use 20 threads, -runThreadN 20
# specify that I want to generate genome, --runMode genomeGenerate
# specify path to save STAR genome directory. Must already exist, --genomeDir
# specify path to genome, --genomeFastaFiles
# specify path to annotation file, --sjdbGTFfile
# specify length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database.
#    Ideally, this length should be equal to the ReadLength-1, where ReadLength is the length of the reads.
#    My reads are 100bp, so I'll use 99 (Default is 100). --sjdbOverhang 99

#STAR \
#--runThreadN 20 \
#--runMode genomeGenerate \
#--genomeDir ${REF}/STAR/ \
#--genomeFastaFiles ${REF}/Paralithodes_platypus_genome.fasta \
#--sjdbGTFfile ${REF}/EVM.out_new.gtf \
#--sjdbOverhang 99

## Build an RSEM index
#rsem-prepare-reference --gtf ${REF}/EVM.out_new.gtf \
#${REF}/Paralithodes_platypus_genome.fasta \
#${REF}/RSEM/Paralithodes_platypus_genome

# ========= Run STAR alignment
# Use most of STAR's default settings, which include (but aren't limited to) ...
# Max number of multiple alignments is 10, if exceeded read is considered unmapped --outFilterMultimapNmax 10
# Min numer of bp overlap to assign read to a gene is 1nt
# Don't use 2-pass mode (not looking for novel splice variants)
# Also ... instruct STAR to output genomic alignments in transcriptomic coordinates (i.e. Aligned.toTranscriptome.out.bam
# this then allows us to use RSEM to quantify fragments per gene in.

# FASTQ=$(ls ${IN}/*.trimmed.R1.v2.fastq.gz | \
# 	awk -F "/" '{print $NF}' | \
# 	awk -F "." '{print $1}')
#
# for sample in ${FASTQ}
# do
# echo "Started mapping ${sample}"
# STAR \
# --runThreadN 20 \
# --genomeDir ${REF}/STAR/ \
# --readFilesIn ${IN}/${sample}.trimmed.R1.v2.fastq.gz ${IN}/${sample}.trimmed.R2.v2.fastq.gz \
# --readFilesCommand gunzip -c \
# --outFilterMultimapNmax 50 \
# --outFileNamePrefix ${OUT}/${sample}. \
# --outSAMtype BAM SortedByCoordinate \
# --quantMode TranscriptomeSAM GeneCounts
# done

# ======= Use RSEM to quantify the number of fragments per gene
# options
# --bam Input file is in BAM format.
# --no-bam-output Do not output any BAM file.
# -p Number of threads to use.
# --paired-end Input reads are paired-end reads.
# --forward-prob Probability of generating a read from the forward strand of a transcript. 1: strand-specific protocol where all (upstream) reads are derived from the forward strand; 0: strand-specific protocol where all (upstream) read are derived from the reverse strand; 0.5: non-strand-specific protocol

BAMS=$(ls ${OUT}/*.Aligned.toTranscriptome.out.bam | \
	awk -F "/" '{print $NF}' | \
	awk -F "." '{print $1}')

for sample in ${BAMS}
do

rsem-calculate-expression \
--no-bam-output \
-p 20 \
--strandedness reverse \
--alignments \
--paired-end \
${OUT}/${sample}.Aligned.toTranscriptome.out.bam \
${REF}/RSEM/Paralithodes_platypus_genome \
${sample}
done >> "${OUT}/RSEM-calc-expr_stout.txt" 2>&1
