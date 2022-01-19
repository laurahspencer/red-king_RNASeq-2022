#!/bin/bash

#SBATCH --job-name=redking_Bowtie2-RSEM
#SBATCH --output=/home/lspencer/2022-redking-OA/sbatch_logs/redking_Bowtie2-RSEM.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 24
#SBATCH -t 10-0:0:0

module load bio/rsem/
module load aligners/bowtie2/2.4.2

BASE=/scratch/lspencer/2022-redking-OA
REF=/home/lspencer/references/bluekingcrab
IN=${BASE}/trimmed
OUT=${BASE}/aligned/bowtie2-rsem

SAMPLES=$(ls ${IN}/*.trimmed.R1.v2.fastq.gz | \
	awk -F "/" '{print $NF}' | \
	awk -F "." '{print $1}')

rsem-prepare-reference \
--gff3 ${REF}/EVM.out_new.gff3 \
	--bowtie2 \
	-p 24 \
	${REF}/Paralithodes_platypus_genome.fasta \
	${REF}/P.plat_rsem_ref \
	&> ${REF}/prep_rsem_ref.log

if [ ! -d ${OUT} ]; then
	mkdir -p ${OUT}
fi

for SAMPLE in ${SAMPLES}
do
	echo ${SAMPLE}
	rsem-calculate-expression \
    	--bowtie2 \
    	-p 24 \
    	--paired-end \
    	--calc-pme \
    	--calc-ci \
    	${IN}/${SAMPLE}.trimmed.R1.v2.fastq.gz \
    	${IN}/${SAMPLE}.trimmed.R2.v2.fastq.gz \
    	${REF}/P.plat_rsem_ref \
    	${OUT}/${SAMPLE} \
    	&> ${OUT}/rsem.${SAMPLE}.log
done
