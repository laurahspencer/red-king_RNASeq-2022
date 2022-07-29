#!/bin/bash

#SBATCH --job-name=genotype-rkc
#SBATCH --output=/home/lspencer/2022-redking-OA/sbatch_logs/genotype-rnaseq-concat-RKC-genome-20220627.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 24
#SBATCH -p himem
#SBATCH -t 30-0:0:0
#SBATCH --mem=1200GB

module load aligners/bowtie2/2.4.2
module load bio/gatk/4.2.0.0
module load bio/samtools/1.15.1
source /home/lspencer/venv/bin/activate
module load bio/vcftools/0.1.16

REF=/home/lspencer/references/redkingcrab #red king crab genome directory
INPUT=/home/lspencer/2022-redking-OA/trimmed/v2 #trimmed reads
OUTPUT=/scratch/lspencer/2022-redking-OA/genotype/genotypes-bowtie-rkc-concat # destination for all output files

# # Move to the directory containing red king crab genome.
# # Bowtie2 looks for index in current directory - easier just to be in it from the get-go
# cd ${REF}

# # Create bowtie2 index for concatenated Oly genome
# echo "Building bowtie2 genome index"
# bowtie2-build \
# --threads 1 \
# Paralithodes.camtschaticus.genome_concat.fa \
# Paralithodes.camtschaticus.genome_concat.fa >> "${OUTPUT}/01-bowtie2-build.txt" 2>&1
#
# # Align trimmed reads to concatenated genome
# echo "Aligning reads"
#
# # Run Bowtie2 over each RNASeq paired sample
# for file in ${INPUT}/*.trimmed.R1.v2.fastq.gz
# do
# sample="$(basename -a ${file} | cut -d "." -f 1)"
# file_R1="${sample}.trimmed.R1.v2.fastq.gz"
# file_R2="${sample}.trimmed.R2.v2.fastq.gz"
# map_file="${sample}.bowtie.concat.sam"
#
# # run Bowtie2 on each file
# bowtie2 \
# -x Paralithodes.camtschaticus.genome_concat.fa \
# --sensitive \
# --threads 24 \
# --no-unal \
# -1 ${INPUT}/${file_R1} \
# -2 ${INPUT}/${file_R2} \
# -S ${OUTPUT}/${map_file}; \
# done >> ${OUTPUT}/02-bowtieout.txt 2>&1
#
# Move to output directory where all files are to be written
cd ${OUTPUT}

# # Convert sam to bam and sort
# echo "Convert aligned .sam to .bam"
#
# for file in *.bowtie.concat.sam
# do
# sample="$(basename -a $file | cut -d "." -f 1)"
# samtools view --threads 24 -b $file | samtools sort -o $sample.sorted.bam
# done >> "03-sam2sortedbam.txt" 2>&1
#
# # Deduplicate using picard (within gatk), output will have duplicates removed
# echo "Deduplicating bams"
#
# # Using Bowtie2 aligned to concatenated genome
# for file in *.sorted.bam
# do
# sample="$(basename -a $file | cut -d "." -f 1)"
#
# gatk MarkDuplicates \
# I=$file \
# O="$sample.dedup.bam" \
# M="$sample.dup_metrics.txt" \
# REMOVE_DUPLICATES=true
# done >> "01_dedup_stout.txt" 2>&1
#
# # Create a FASTA sequence dictionary file for genome (needed by gatk)
# echo "Creating sequence dictionary (.dict)"
# gatk CreateSequenceDictionary \
# -R ${REF}/Paralithodes.camtschaticus.genome_concat.fa \
# -O ${REF}/Paralithodes.camtschaticus.genome_concat.dict \
#  >> "02-CreateSequenceDictionary.txt" 2>&1
#
# # Split reads spanning splicing events
# echo "Splitting reads spanning splice junctions (SplitNCigarReads)"
# for file in *dedup.bam
# do
# sample="$(basename -a $file | cut -d "." -f 1)"
#
# # split CigarN reads
# # NOTE: with concatenated genome the contigs are too long for indexing within gatk, so need to specify `--create-output-bam false` (we don't need them anyway for this intermediate bam file)
# gatk SplitNCigarReads \
# --create-output-bam-index false \
# -R ${REF}/Paralithodes.camtschaticus.genome_concat.fa \
# -I $file \
# -O $sample.dedup-split.bam
# done >> "03-CigarNSplit_stout.txt" 2>&1
#
# # Remove interim .bam files to conserve space
# rm *.*dedup.bam
#
# # Add read group ID to bams (needed by gatk)
# echo "Adding read group to bams"
# for file in *dedup-split.bam
# do
# sample="$(basename -a $file | cut -d "." -f 1)"
#
# # add read group info to headers, specifying sample names
# gatk AddOrReplaceReadGroups \
# I=$sample.dedup-split.bam \
# O=$sample.dedup-split-RG.bam \
# RGID=1 \
# RGLB=$sample \
# RGPL=ILLUMINA \
# RGPU=unit1 \
# RGSM=$sample
# done >> "04-AddReadGroup_stout.txt" 2>&1
#
# # Remove interim .bam files to conserve space
# rm *.dedup-split.bam

# Index the final .bam files (that have been deduplicated, split, read-group added)
echo "Indexing variant-call ready .bam files"
for file in *dedup-split-RG.bam
do
samtools index $file
done >> "05-index-bams.txt" 2>&1

# create interval list (just a list of all contigs in genome)
# This will be used in HaplotypeCaller and GenomicsDBImport to increase speed
# Note: the intervals file requires a specific name - e.g. for .bed format, it MUST be "intervals.bed"
echo "Creating intervals list"
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' ${REF}/Paralithodes.camtschaticus.genome_concat.fa.fai > intervals.bed

# Call variants
echo "Calling variants using HaplotypeCaller"
for file in *dedup-split-RG.bam
do
sample="$(basename -a $file | cut -d "." -f 1)"

gatk HaplotypeCaller \
-R ${REF}/Paralithodes.camtschaticus.genome_concat.fa \
-I $sample.dedup-split-RG.bam \
-O $sample.variants.g.vcf \
-L intervals.bed \
--native-pair-hmm-threads 20 \
-ERC GVCF
done >> "06-HaplotypeCaller_stout.txt" 2>&1

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
--java-options '-Xmx1200g' \
--genomicsdb-workspace-path GenomicsDB/ \
-L intervals.bed \
--sample-name-map sample_map.txt \
--reader-threads 20 >> "07-GenomicsDBImport_stout.txt" 2>&1

# Joint genotype
echo "Joint genotyping"
gatk GenotypeGVCFs \
-R ${REF}/Paralithodes.camtschaticus.genome_concat.fa \
-V gendb://GenomicsDB \
-O rkc_rnaseq_genotypes.vcf.gz \
>> "08-GenotypeGVCFs_stout.txt" 2>&1

# Hard filter variants
echo "Hard filtering variants"
gatk VariantFiltration \
-R ${REF}/Paralithodes.camtschaticus.genome_concat.fa \
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
--filter "AF < 0.30" >> "09-Filter-variants_stout.txt" 2>&1

# Select only SNPs that pass filtering
echo "Selecting SNPs that pass fitering"
gatk SelectVariants \
-R ${REF}/Paralithodes.camtschaticus.genome_concat.fa \
-V rkc_rnaseq_genotypes-filtered.vcf.gz \
--exclude-filtered TRUE \
--select-type-to-include SNP \
-O rkc_rnaseq_genotypes-filtered-true.vcf.gz \
 >> "10-SelectVariants_stout.txt" 2>&1

 # Create another vcf of SNPs filtered for loci with max 10%, 15%, and 20% missing rate, and remove loci with <5% minor allele frequency

 vcftools --gzvcf \
 "rkc_rnaseq_genotypes-filtered-true.vcf.gz" \
 --max-missing 0.75 --maf 0.05 --recode --recode-INFO-all --out \
 "rkc_rnaseq_genotypes-final_miss25"

 vcftools --gzvcf \
 "rkc_rnaseq_genotypes-filtered-true.vcf.gz" \
 --max-missing 0.80 --maf 0.05 --recode --recode-INFO-all --out \
 "rkc_rnaseq_genotypes-final_miss20"

 vcftools --gzvcf \
 "rkc_rnaseq_genotypes-filtered-true.vcf.gz" \
 --max-missing 0.85 --maf 0.05 --recode --recode-INFO-all --out \
 "rkc_rnaseq_genotypes-final_miss15"

echo "complete!"
