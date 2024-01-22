#!/bin/bash

#SBATCH --job-name=redking_feature-counts
#SBATCH --output=/home/lspencer/2022-redking-OA/sbatch_logs/extract-chimeras-bam.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 5-0:0:0

module load bio/samtools/1.11

BASE=/home/lspencer/2022-redking-OA/counts

for file in ${BASE}/*.sorted.bam.featureCounts.bam
do
sample="$(basename -a $file | cut -d "." -f 1)"

# Extract all rows that contain the tag "Unassigned_Chimera"
samtools view --threads 20 $file | \
grep "Unassigned_Chimera" > ${BASE}/$sample.chimeras.txt

samtools view -H --threads 20 $file > ${BASE}/$sample.header.txt

cat ${BASE}/$sample.header.txt ${BASE}/$sample.chimeras.txt > ${BASE}/$sample.chimeras.sam

done
rclone --progress --tpslimit 10 --fast-list copy /home/lspencer/2022-redking-OA/counts/*.chimeras.sam remote:Testing/
