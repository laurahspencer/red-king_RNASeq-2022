#!/bin/bash

#SBATCH --job-name=salmon-pseudoalign
#SBATCH --output=/home/lspencer/2022-redking-OA/sbatch_logs/salmon-pseudoalign.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 7-0:0:0

module load aligners/salmon/1.4.0

REF=/home/lspencer/references/bluekingcrab
IN=/home/lspencer/2022-redking-OA/trimmed/v2
OUT=/scratch/lspencer/2022-redking-OA/aligned/salmon

# Build a decoy-aware salmon index from transcriptome + genome

# Salmon indexing requires the names of the genome targets, which is extractable by using the grep command:

#grep "^>" < ${REF}/Paralithodes_platypus_genome.fasta | cut -d " " -f 1 > ${REF}/salmon/decoys.txt
#sed -i.bak -e 's/>//g' ${REF}/salmon/decoys.txt

# Along with the list of decoys salmon also needs the concatenated transcriptome and genome reference file for index. NOTE: the genome targets (decoys) should come after the transcriptome targets in the reference
#cat ${REF}/EVM.out_new.cds ${REF}/Paralithodes_platypus_genome.fasta > ${REF}/salmon/Paralithodes_platypus_gentrome.fa

### Salmon indexing
salmon index -t ${REF}/EVM.out_new.cds \
-i ${REF}/salmon/salmon_index \
-k 31 \
-p 20
#salmon index -t ${REF}/salmon/Paralithodes_platypus_gentrome.fa \
#-d ${REF}/salmon/decoys.txt \

# Run salmon
FASTQ=$(ls ${IN}/*.trimmed.R1.v2.fastq.gz | \
	awk -F "/" '{print $NF}' | \
	awk -F "." '{print $1}')
for sample in ${FASTQ}
do

salmon quant \
-i ${REF}/salmon/salmon_index \
--libType A \
-1 ${IN}/${sample}.trimmed.R1.v2.fastq.gz \
-2 ${IN}/${sample}.trimmed.R2.v2.fastq.gz \
--validateMappings \
--gcBias \
--threads 20 \
--output ${OUT}/${sample}.salmon.counts
done
