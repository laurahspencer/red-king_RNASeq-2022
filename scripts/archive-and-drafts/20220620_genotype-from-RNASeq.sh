#!/bin/bash

#SBATCH --job-name=genotype-rkc
#SBATCH --output=/home/lspencer/2022-redking-OA/sbatch_logs/RNAseq-genotyping-Bowtie2-RKC-genome_from-step-07.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 24
#SBATCH -p himem
#SBATCH -t 30-0:0:0
#SBATCH --mem=1200GB

module load bio/gatk/4.2.0.0
module load bio/samtools/1.11
source /home/lspencer/venv/bin/activate
module load bio/vcftools/0.1.16

#REF=/home/lspencer/references/bluekingcrab
REF=/home/lspencer/references/redkingcrab
#INPUT=/home/lspencer/2022-redking-OA/aligned/star #to use star aligned
#INPUT=/home/lspencer/2022-redking-OA/aligned/bowtie2-sorted/ #to use bowtie2 aligned to BKC genome
INPUT=/home/lspencer/2022-redking-OA/aligned/bowtie-rkc #to use bowtie2 aligned to RKG genome
#OUTPUT=/scratch/lspencer/2022-redking-OA/genotypes-star
OUTPUT=/scratch/lspencer/2022-redking-OA/genotypes-bowtie-rkc

# Starting this pipeline with aligned and coordinate-sorted .bam files

# Move to output directory
cd ${OUTPUT}

# # Deduplicate using picard (within gatk), output will have duplicates removed
# echo "Deduplicating bams"
#
# ## Using STAR aligned
# #for file in ${INPUT}/*.Aligned.sortedByCoord.out.bam
#
# # Using Bowtie2 aligned
# for file in ${INPUT}/*.sorted.bam
# do
# sample="$(basename -a $file | cut -d "." -f 1)"
#
# gatk MarkDuplicates \
# I=$file \
# O="${OUTPUT}/$sample.dedup.bam" \
# M="${OUTPUT}/$sample.dup_metrics.txt" \
# REMOVE_DUPLICATES=true
# done >> "${OUTPUT}/01_dedup_stout.txt" 2>&1
#
# # Create a FASTA sequence dictionary file for genome (needed by gatk)
# echo "Creating sequence dictionary (.dict)"
# gatk CreateSequenceDictionary \
# -R ${REF}/Paralithodes.camtschaticus.genome.fasta \
# -O ${REF}/Paralithodes.camtschaticus.genome.dict \
#  >> "${OUTPUT}/02-CreateSequenceDictionary.txt" 2>&1
#
# # Split reads spanning splicing events
# echo "Splitting reads spanning splice junctions (SplitNCigarReads)"
# for file in *dedup.bam
# do
# sample="$(basename -a $file | cut -d "." -f 1)"
#
# # split CigarN reads
# gatk SplitNCigarReads \
# -R ${REF}/Paralithodes.camtschaticus.genome.fasta \
# -I $file \
# -O $sample.dedup-split.bam
# done >> "${OUTPUT}/03-CigarNSplit_stout.txt" 2>&1
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
# done >> "${OUTPUT}/04-AddReadGroup_stout.txt" 2>&1
#
# # Remove interim .bam files to conserve space
# rm *.dedup-split.bam
# rm *.dedup-split.bam.bai
#
# # Index the final .bam files (that have been deduplicated, split, read-group added)
# echo "Indexing variant-call ready .bam files"
# for file in *dedup-split-RG.bam
# do
# samtools index $file
# done >> "${OUTPUT}/05-index-bams.txt" 2>&1
#
# # create interval list (just a list of all contigs in genome)
# # This will be used in HaplotypeCaller and GenomicsDBImport to increase speed
# # Note: the intervals file requires a specific name - e.g. for .bed format, it MUST be "intervals.bed"
# echo "Creating intervals list"
# awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' ${REF}/Paralithodes.camtschaticus.genome.fasta.fai > intervals.bed
#
# # Call variants
# echo "Calling variants using HaplotypeCaller"
# for file in *dedup-split-RG.bam
# do
# sample="$(basename -a $file | cut -d "." -f 1)"
#
# gatk HaplotypeCaller \
# -R ${REF}/Paralithodes.camtschaticus.genome.fasta \
# -I $sample.dedup-split-RG.bam \
# -O $sample.variants.g.vcf \
# -L intervals.bed \
# --native-pair-hmm-threads 20 \
# -ERC GVCF
# done >> "${OUTPUT}/06-HaplotypeCaller_stout.txt" 2>&1
#
# # create sample map of all gvcfs
# rm sample_map.txt  #remove if already exists
# echo "Creating sample map of all gvcfs"
# for file in *variants.g.vcf
# do
# sample="$(basename -a $file | cut -d "." -f 1)"
# echo -e "$sample\t$file" >> sample_map.txt
# done

# Aggregate single-sample GVCFs into GenomicsDB
# Note: the intervals file requires a specific name - e.g. for .bed format, it MUST be "intervals.bed"
echo "Aggregating single-sample GVCFs into GenomicsDB"
rm -r GenomicsDB/ #can't already have a GenomicsDB directory, else will fail
gatk GenomicsDBImport \
--java-options '-Xmx1200g' \
--genomicsdb-workspace-path GenomicsDB/ \
-L intervals.bed \
--sample-name-map sample_map.txt \
--reader-threads 20 >> "${OUTPUT}/07-GenomicsDBImport_stout.txt" 2>&1

# Joint genotype
echo "Joint genotyping"
gatk GenotypeGVCFs \
-R ${REF}/Paralithodes.camtschaticus.genome.fasta \
-V gendb://GenomicsDB \
-O rkc_rnaseq_genotypes.vcf.gz \
>> "${OUTPUT}/08-GenotypeGVCFs_stout.txt" 2>&1

# Hard filter variants
echo "Hard filtering variants"
gatk VariantFiltration \
-R ${REF}/Paralithodes.camtschaticus.genome.fasta \
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
--filter "AF < 0.30" >> "${OUTPUT}/09-Filter-variants_stout.txt" 2>&1

# Select only SNPs that pass filtering
echo "Selecting SNPs that pass fitering"
gatk SelectVariants \
-R ${REF}/Paralithodes.camtschaticus.genome.fasta \
-V rkc_rnaseq_genotypes-filtered.vcf.gz \
--exclude-filtered TRUE \
--select-type-to-include SNP \
-O rkc_rnaseq_genotypes-filtered-true.vcf.gz \
 >> "${OUTPUT}/10-SelectVariants_stout.txt" 2>&1

 # Create another vcf of SNPs filtered for loci with max 10%, 15%, and 20% missing rate, and remove loci with <5% minor allele frequency

 vcftools --gzvcf \
 "${OUTPUT}/rkc_rnaseq_genotypes-filtered-true.vcf.gz" \
 --max-missing 0.75 --maf 0.05 --recode --recode-INFO-all --out \
 "${OUTPUT}/rkc_rnaseq_genotypes-final_miss25"

 vcftools --gzvcf \
 "${OUTPUT}/rkc_rnaseq_genotypes-filtered-true.vcf.gz" \
 --max-missing 0.80 --maf 0.05 --recode --recode-INFO-all --out \
 "${OUTPUT}/rkc_rnaseq_genotypes-final_miss20"

 vcftools --gzvcf \
 "${OUTPUT}/rkc_rnaseq_genotypes-filtered-true.vcf.gz" \
 --max-missing 0.85 --maf 0.05 --recode --recode-INFO-all --out \
 "${OUTPUT}/rkc_rnaseq_genotypes-final_miss15"

echo "complete!"
