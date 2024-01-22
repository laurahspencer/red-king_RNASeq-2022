#!/bin/bash

#SBATCH --job-name=RNAseq-genotyping
#SBATCH --output=/home/lspencer/2022-redking-OA/sbatch_logs/RNAseq-genotyping-testing.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 21-0:0:0

module load bio/gatk/4.2.0.0
module load bio/samtools/1.11
source /home/lspencer/venv/bin/activate

#REF=/home/lspencer/references/bluekingcrab
REF=/home/lspencer/references/redkingcrab
#INPUT=/home/lspencer/2022-redking-OA/aligned/star #to use star aligned
#INPUT=/home/lspencer/2022-redking-OA/aligned/bowtie2-sorted/ #to use bowtie2 aligned to BKC genome
INPUT=/home/lspencer/2022-redking-OA/aligned/bowtie-rkc #to use bowtie2 aligned to RKG genome
#OUTPUT=/scratch/lspencer/2022-redking-OA/genotypes-star
OUTPUT=/scratch/lspencer/2022-redking-OA/genotypes-bowtie-rkc

# Starting this pipeline with pre-processed bam files from upstream of the pipeline

# Move to output directory
cd ${OUTPUT}

# Call variants
echo "Calling variants using HaplotypeCaller"
for file in *dedup-split-RG.bam
do
sample="$(basename -a $file | cut -d "." -f 1)"

gatk HaplotypeCaller \
-R ${REF}/Paralithodes_platypus_genome.fasta \
-I $sample.dedup-split-RG.bam \
-O $sample.variants.g.vcf \
-L intervals.bed \
--native-pair-hmm-threads 20 \
-ERC GVCF
done >> "${OUTPUT}/06-HaplotypeCaller_stout.txt" 2>&1

# create sample map of all gvcfs
rm sample_map.txt  #remove if already exists
echo "Creating sample map of all gvcfs"
for file in *variants.g.vcf
do
sample="$(basename -a $file | cut -d "." -f 1)"
echo -e "$sample\t$file" >> sample_map.txt
done

# Aggregate single-sample GVCFs into GenomicsDB
# Note: the intervals file requires a specific name - e.g. for .bed format, it MUST be "intervals.bed"
echo "Aggregating single-sample GVCFs into GenomicsDB"
rm -r GenomicsDB/ #can't already have a GenomicsDB directory, else will fail
gatk GenomicsDBImport \
--genomicsdb-workspace-path GenomicsDB/ \
-L intervals.bed \
--sample-name-map sample_map.txt \
--reader-threads 20 >> "${OUTPUT}/07-GenomicsDBImport_stout.txt" 2>&1

# Joint genotype
echo "Joint genotyping"
gatk GenotypeGVCFs \
-R ${REF}/Paralithodes_platypus_genome.fasta \
-V gendb://GenomicsDB \
-O rkc_rnaseq_genotypes.vcf.gz \
>> "${OUTPUT}/08-GenotypeGVCFs_stout.txt" 2>&1

# Hard filter variants
echo "Hard filtering variants"
gatk VariantFiltration \
-R ${REF}/Paralithodes_platypus_genome.fasta \
-V rkc_rnaseq_genotypes.vcf.gz \
-O rkc_rnaseq_genotypes-filtered.vcf.gz \
--filter-name "FS" \
--filter "FS > 60.0" \
--filter-name "QD" \
--filter "QD < 2.0" \
--filter-name "QUAL30" \
--filter "QUAL < 30.0" \
--filter-name "SOR3" \
--filter "SOR > 3.0" \
--filter-name "DP15" \
--filter "DP < 15" \
--filter-name "DP150" \
--filter "DP > 150" \
--filter-name "AF30" \
--filter "AF < 0.30" >> "${OUTPUT}/09-GenotypeGVCFs_stout.txt" 2>&1

# Select only SNPs that pass filtering
echo "Selecting SNPs that pass fitering"
gatk SelectVariants \
-R ${REF}/Paralithodes_platypus_genome.fasta \
-V rkc_rnaseq_genotypes-filtered.vcf.gz \
--exclude-filtered TRUE \
--select-type-to-include SNP \
-O rkc_rnaseq_genotypes-filtered-true.vcf.gz \
 >> "${OUTPUT}/10-SelectVariants_stout.txt" 2>&1

echo "complete!"
