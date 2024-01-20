#!/bin/bash

#SBATCH --job-name=redking_htseq
#SBATCH --output=/home/lspencer/2022-redking-OA/sbatch_logs/htseq-from-STAR.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 10-0:0:0

IN=/scratch/lspencer/2022-redking-OA/trimmed
OUT=/scratch/lspencer/2022-redking-OA/aligned/star
GTF=/home/lspencer/references/bluekingcrab

# For multimapped fragments, count all instances

for file in ${INPUT}/*.Aligned.sortedByCoord.out.bam
do
sample="$(basename -a $file | cut -d "." -f 1)"

htseq-count \
--format=bam \
--order=pos \
--stranded=yes \
--nonunique=all \
--counts_output=${OUT}/$file-htseq-multimapped-all \
--nprocesses=20 \
${IN}/*.bam \
${GTF}/EVM.out_new.gtf




# For multimapped fragments, count them at random
htseq-count \
--format=bam \
--order=pos \
--stranded=yes \
--nonunique=random \
--counts_output=${OUT}/star-htseq-multimapped-random \
--nprocesses=20 \
${IN}/*.bam \
${GTF}/EVM.out_new.gtf

# For multimapped fragments, count them fractionally
htseq-count \
--format=bam \
--order=pos \
--stranded=yes \
--nonunique=fraction \
--counts_output=${OUT}/star-htseq-multimapped-fraction \
--nprocesses=20 \
${IN}/*.bam \
${GTF}/EVM.out_new.gtf

done
