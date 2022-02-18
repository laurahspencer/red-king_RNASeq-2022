#!/bin/bash
## Job Name
#SBATCH --job-name=QuantSeq-genotyping
## Allocation Definition 
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes 
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=7-00:00:00
## Memory per node
#SBATCH --mem=100G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lhs3@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/srlab/lhs3/jobs/20210314_Olurida_QuantSeq-Genotype-all/

# Set program paths 
bowtie2="/gscratch/srlab/programs/bowtie2-2.1.0/"
gatk="/gscratch/srlab/programs/gatk-4.1.9.0/gatk"
samtools="/gscratch/srlab/programs/samtools-1.10/samtools"

# Here are various important paths
# trimmed QuantSeq fastq files: /gscratch/srlab/lhs3/data/QuantSeq2020
# O.lurida genome: /gscratch/srlab/lhs3/data/Olurida_QuantSeq2020-trimmed/
# Aligned QuantSeq reads: /gscratch/scrubbed/lhs3/20210314_Olurida_QuantSeq-Genotype-all/mapped/
# gatk files:  /gscratch/scrubbed/lhs3/20210314_Olurida_QuantSeq-Genotype-all/gatk/

homedir="/gscratch/srlab/lhs3/"
jobdir="/gscratch/srlab/lhs3/jobs/20210314_Olurida_QuantSeq-Genotype-all/"
outputdir="/gscratch/scrubbed/lhs3/20210314_Olurida_QuantSeq-Genotype-all/"

# Make directories that don't yet exist
mkdir ${homedir}data/Olurida_v081_concat/
mkdir ${outputdir}mapped/
mkdir ${outputdir}gatk/

# Concatenate Oly genome contigs into ~30 super-contigs
echo "Starting genome concatenation"
awk -v i=5314 -f ${homedir}code/merge_contigs.awk.txt ${homedir}data/Olurida_v081.fa > ${homedir}data/Olurida_v081_concat/Olurida_v081_concat.fa

# Edit fasta headers in place from concatenated IDs to the range (e.g. "Contig0 Contig1 Contig2" to "Contig0:2")

# First remove all horizontal white spaces between contigIDs, save to new fasta file 
cat ${homedir}data/Olurida_v081_concat/Olurida_v081_concat.fa | tr -d "[:blank:]" > ${homedir}data/Olurida_v081_concat/Olurida_v081_concat_rehead.fa #remove all horizontal white spaces between contigIDs

n_contigs_old=$(grep ">" ${homedir}data/Olurida_v081.fa | wc -l) #number of contigs in original fasta 
d="${n_contigs_old//[^[:digit:]]/}" #Return ID of last/highest number contig
digits=$(echo ${#d}) #number of digits in the last contig
n_contigs_new=$(grep ">" ${homedir}data/Olurida_v081_concat/Olurida_v081_concat.fa | wc -l) #number of contigs in concatenated genome
n=$(seq 1 1 $n_contigs_new) #vector of 1:# of contigs in concat. genome

# loop over each new contig header, paste first and last ID # together and replace concat header with new ID 
for i in $n
do
old=$(grep ">" ${homedir}data/Olurida_v081_concat/Olurida_v081_concat_rehead.fa | sed "${i}q;d")
start=$(echo $old | grep -o -E '[0-9]{1,6}' | head -n 1) #If needed, change the 6 to match the value in $digits 
startID=$(echo ">Contig"$start)
end=$(echo $old | grep -o -E '[0-9]{1,6}' | tail -n 1) #If needed, change the 6 to match the value in $digits 
new=$(echo ">Contig"$start"_"$end)
sed -i.bak "s|$startID.*|$new|" ${homedir}data/Olurida_v081_concat/Olurida_v081_concat_rehead.fa # edit fasta in place. 
done

# Index concatenated genome 
echo "Indexing concatenated and re-headed fasta" 
${samtools} faidx ${homedir}data/Olurida_v081_concat/Olurida_v081_concat_rehead.fa

# Create bowtie2 index for concatenated Oly genome 
echo "Building bowtie2 genome index"
${bowtie2}bowtie2-build \
${homedir}data/Olurida_v081_concat/Olurida_v081_concat_rehead.fa \
${homedir}data/Olurida_v081_concat/Olurida_v081_concat_rehead.fa >> "${jobdir}01-bowtie2-build.txt" 2>&1

# Align trimmed reads to concatenated genome 
cd ${homedir}data/QuantSeq2020/
echo "Aligning reads" 

for file in *.trim.fastq
do
sample="$(basename -a $file | cut -d "." -f 1)"

# run Bowtie2 on each file
${bowtie2}bowtie2 \
--local \
-x ${homedir}data/Olurida_v081_concat/Olurida_v081_concat_rehead.fa \
--sensitive-local \
--threads 8 \
--no-unal \
-k 5 \
-U $file \
-S ${outputdir}mapped/$sample.bowtie.sam; \
done >> ${jobdir}02-bowtieout.txt 2>&1

cd ${outputdir}mapped/

# Convert sam to bam and sort  
echo "Convert aligned .sam to .bam" 

for file in *.bowtie.sam
do
sample="$(basename -a $file | cut -d "." -f 1)"
${samtools} view -b $file | ${samtools} sort -o $sample.sorted.bam
done >> "${jobdir}03-sam2sortedbam.txt" 2>&1

# Deduplicate using picard (within gatk), output will have duplicates removed 

echo "Deduplicating bams"
for file in *sorted.bam
do
sample="$(basename -a $file | cut -d "." -f 1)"

${gatk} MarkDuplicates \
I=$file \
O="${outputdir}gatk/$sample.dedup.bam" \
M="${outputdir}gatk/$sample.dup_metrics.txt" \
REMOVE_DUPLICATES=true
done >> "${jobdir}04-dedup_stout.txt" 2>&1

# Move to gatk directory and execute gatk tools 
cd ${outputdir}gatk/

# Create a FASTA sequence dictionary file for O.lurida genome (needed by gatk)
echo "Creating sequence dictionary (.dict)" 
${gatk} CreateSequenceDictionary \
-R ${homedir}data/Olurida_v081_concat/Olurida_v081_concat_rehead.fa \
-O ${homedir}data/Olurida_v081_concat/Olurida_v081_concat_rehead.dict \
 >> "${jobdir}05-CreateSequenceDictionary.txt" 2>&1

# Split reads spanning splicing events 
echo "Splitting reads spanning splice junctions (SplitNCigarReads)"
for file in *dedup.bam
do
sample="$(basename -a $file | cut -d "." -f 1)"

# split CigarN reads
${gatk} SplitNCigarReads \
-R ${homedir}data/Olurida_v081_concat/Olurida_v081_concat_rehead.fa \
-I $file \
-O $sample.dedup-split.bam
done >> "${jobdir}06-CigarNSplit_stout.txt" 2>&1

# Add read group ID to bams (needed by gatk)
echo "Adding read group to bams" 
for file in *dedup-split.bam
do
sample="$(basename -a $file | cut -d "." -f 1)"

# add read group info to headers, specifying sample names 
${gatk} AddOrReplaceReadGroups \
I=$sample.dedup-split.bam \
O=$sample.dedup-split-RG.bam \
RGID=1 \
RGLB=$sample \
RGPL=ILLUMINA \
RGPU=unit1 \
RGSM=$sample
done >> "${jobdir}07-AddReadGroup_stout.txt" 2>&1

# Index the final .bam files (that have been deduplicated, split, read-group added)
echo "Indexing variant-call ready .bam files"
for file in *dedup-split-RG.bam
do
${samtools} index $file 
done >> "${jobdir}08-index-bams.txt" 2>&1

# Call variants 
echo "Calling variants using HaplotypeCaller"
for file in *dedup-split-RG.bam
do
sample="$(basename -a $file | cut -d "." -f 1)"

${gatk} HaplotypeCaller \
-R ${homedir}data/Olurida_v081_concat/Olurida_v081_concat_rehead.fa \
-I $sample.dedup-split-RG.bam \
-O $sample.variants.g.vcf \
-ERC GVCF
done >> "${jobdir}09-HaplotypeCaller_stout.txt" 2>&1

# create sample map of all gvcfs
rm sample_map.txt  #remove if already exists 
echo "Creating sample map of all gvcfs"
for file in *variants.g.vcf
do
sample="$(basename -a $file | cut -d "." -f 1)"
echo -e "$sample\t$file" >> sample_map.txt
done

# create interval list (just a list of all contigs in genome)
echo "Creating intervals list"
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' ${homedir}data/Olurida_v081_concat/Olurida_v081_concat_rehead.fa.fai > intervals.bed 

# Aggregate single-sample GVCFs into GenomicsDB
# Note: the intervals file requires a specific name - e.g. for .bed format, it MUST be "intervals.bed"
echo "Aggregating single-sample GVCFs into GenomicsDB"
rm -r GenomicsDB/ #can't already have a GenomicsDB directory, else will fail 
${gatk} GenomicsDBImport \
--genomicsdb-workspace-path GenomicsDB/ \
-L intervals.bed \
--sample-name-map sample_map.txt \
--reader-threads 40 >> "${jobdir}10-GenomicsDBImport_stout.txt" 2>&1

# Joint genotype 
echo "Joint genotyping"
${gatk} GenotypeGVCFs \
-R ${homedir}data/Olurida_v081_concat/Olurida_v081_concat_rehead.fa \
-V gendb://GenomicsDB \
-O Olurida_QuantSeq2020_genotypes.vcf.gz \
>> "${jobdir}11-GenotypeGVCFs_stout.txt" 2>&1

# Hard filter variants 
echo "Hard filtering variants"
${gatk} VariantFiltration \
-R ${homedir}data/Olurida_v081_concat/Olurida_v081_concat_rehead.fa \
-V Olurida_QuantSeq2020_genotypes.vcf.gz \
-O Olurida_QuantSeq2020_genotypes-filtered.vcf.gz \
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
--filter "AF < 0.30" >> "${jobdir}12-GenotypeGVCFs_stout.txt" 2>&1

# Select only SNPs that pass filtering
echo "Selecting SNPs that pass fitering"
${gatk} SelectVariants \
-R ${homedir}data/Olurida_v081_concat/Olurida_v081_concat_rehead.fa \
-V Olurida_QuantSeq2020_genotypes-filtered.vcf.gz \
--exclude-filtered TRUE \
--select-type-to-include SNP \
-O Olurida_QuantSeq2020_genotypes-filtered-true.vcf.gz \
 >> "${jobdir}13-SelectVariants_stout.txt" 2>&1

echo "complete!"