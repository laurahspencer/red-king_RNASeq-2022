#!/bin/bash

#SBATCH --job-name=STAR
#SBATCH --output=/home/lspencer/2022-redking-OA/sbatch_logs/align-STAR-rkc-genome.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -p himem
#SBATCH --mem=610GB
#SBATCH -t 10-0:0:0

module load aligners/star/2.7.10a

REF=/home/lspencer/references/redkingcrab
IN=/home/lspencer/2022-redking-OA/trimmed/v2
OUT=/scratch/lspencer/2022-redking-OA/aligned/star-rkc

# ========= Build STAR genome index
# Use 20 threads, -runThreadN 20
# specify that I want to generate genome, --runMode genomeGenerate
# specify path to save STAR genome directory. Must already exist, --genomeDir
# specify path to genome, --genomeFastaFiles
# specify path to annotation file, --sjdbGTFfile
# specify length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database.
#    Ideally, this length should be equal to the ReadLength-1, where ReadLength is the length of the reads.
#    My reads are 100bp, so I'll use 99 (Default is 100). --sjdbOverhang 99

STAR \
--runThreadN 20 \
--runMode genomeGenerate \
--limitGenomeGenerateRAM 601057045429 \
--genomeDir ${REF}/STAR \
--genomeFastaFiles ${REF}/Paralithodes.camtschaticus.genome.fasta \
--sjdbGTFfile ${REF}/Paralithodes.camtschaticus.genome.gtf \
--sjdbOverhang 99 \
done

# ========= Run STAR alignment
# Use most of STAR's default settings, which include (but aren't limited to) ...
# Max number of multiple alignments is 10, if exceeded read is considered unmapped --outFilterMultimapNmax 10
# Min numer of bp overlap to assign read to a gene is 1nt
# Also ... count number of reads per gene while mapping, using --quantMode GeneCounts
# I'd like to use DEXSeq for differential exon usage (i.e. splice differences) analysis. So, instead of outputting sorted bam files, I output SAM files

SAMPLES=$(ls ${IN}/*.trimmed.R1.v2.fastq.gz | \
	awk -F "/" '{print $NF}' | \
	awk -F "." '{print $1}')

for sample in ${SAMPLES}
do
echo "Started mapping ${sample}"
STAR \
--runThreadN 20 \
--genomeDir ${REF}/STAR \
--readFilesIn ${IN}/${sample}.trimmed.R1.v2.fastq.gz ${IN}/${sample}.trimmed.R2.v2.fastq.gz \
--readFilesCommand gunzip -c \
--outFilterMultimapNmax 50 \
--outFileNamePrefix ${OUT}/${sample}. \
#--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
&> ${OUT}/star-rkc.${sample}.log
done
