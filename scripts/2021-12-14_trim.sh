#!/bin/bash

#SBATCH --job-name=redking_raw-trim
#SBATCH --output=/home/lspencer/2022-redking-OA/sbatch_logs/redking_trim.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 8
#SBATCH -t 1000

# This script is for trimming raw (but concatenated) RNA-Seq data and
# filtering for quality and length.
# It has been slightly adapted from code written by Giles Goetz

module load bio/fastqc
source /home/lspencer/venv/bin/activate

IN=/share/nwfsc/ggoetz/red_king_crab/illumina
OUT=/home/lspencer/2022-redking-OA/trimmed
FASTQC=/home/lspencer/2022-redking-OA/fastqc/trimmed/v2
VER=2

SAMPLES=$(ls ${IN}/*_R1.fastq.gz | \
awk -F "/" '{print $NF}' | \
awk -F "." '{print $1}' | \
sed -e 's/_R1//')

for sample in ${SAMPLES}
do
    # Trimming the Illumina adapters
    # Quality-trim  5’ end with cutoff=20 & 3’ end with cutoff=15
    # Trimming out leftover N's
    # Filtering out sequences shorter then 50bp
    cutadapt \
        -o ${OUT}/${sample}.trimmed.R1.v${VER}.fastq.gz \
        -p ${OUT}/${sample}.trimmed.R2.v${VER}.fastq.gz \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -q 20,15 \
        -m 50 \
        --trim-n \
        --cores=8 \
        ${IN}/${sample}_R1.fastq.gz \
        ${IN}/${sample}_R2.fastq.gz \
        &> ${OUT}/cutadapt.${sample}.v${VER}.log

# Run fastqc on trimmed data files
    fastqc \
        --threads 2 \
        -o ${FASTQC} \
        ${OUT}/${sample}.trimmed.R1.v${VER}.fastq.gz \
        ${OUT}/${sample}.trimmed.R2.v${VER}.fastq.gz \
        &> ${FASTQC}/fastqc.${sample}.v${VER}.log
done

# Run multiqc to summarize fastqc reports
multiqc \
${FASTQC} \
--outdir ${FASTQC}
