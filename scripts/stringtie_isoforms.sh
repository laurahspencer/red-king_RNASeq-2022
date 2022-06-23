#!/bin/bash

#SBATCH --job-name=stringtie
#SBATCH --output=/home/lspencer/2022-redking-OA/sbatch_logs/stringtie.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 24
#SBATCH -t 21-0:0:0
#SBATCH -p himem
#SBATCH --mem=1400GB

# Script produced by Sam White, https://github.com/RobertsLab/sams-notebook/blob/master/sbatch_scripts/20220225_cvir_stringtie_GCF_002022765.2_isoforms.sh
## Script using Stringtie with Red king crab draft genome
## Expects FastQ input filenames to match <sample name>.trimmed.R[12].v2.fastq.gz

# Load modules
source /home/lspencer/venv/bin/activate
module load bio/hisat2/2.2.1
module load bio/samtools/1.11
module load bio/stringtie/2.2.0

###################################################################################
# These variables need to be set by user

## Specify and move to the working directory for this job
cd /scratch/lspencer/2022-redking-OA/stringtie

# Set number of CPUs to use
threads=24

ref_dir="/home/lspencer/references/redkingcrab/hisat2"

## Assign Variables
# Index name for Hisat2 use
# Needs to match index naem used in previous Hisat2 indexing step
genome_index_name="Paralithodes.camtschaticus"

# Location of Hisat2 index files
# Must keep variable name formatting, as it's used by HiSat2
HISAT2_INDEXES=${ref_dir}
export HISAT2_INDEXES

# Input/output files
genome_index_dir="/home/lspencer/references/redkingcrab/hisat2"
genome_gff="/home/lspencer/references/redkingcrab/Paralithodes.camtschaticus.genome.gff"
fastq_dir="/home/lspencer/2022-redking-OA/trimmed/v2/"
gtf_list="gtf_list.txt"
merged_bam="stringtie_GCF_sorted-bams-merged.bam"

# Declare associative array of sample names and metadata
declare -A samples_associative_array=()

# Set total number of samples (NOT number of FastQ files)
total_samples=43

# Programs associative array
declare -A programs_array
programs_array=(
[hisat2]="hisat2" \
[samtools_index]="samtools index" \
[samtools_merge]="samtools merge" \
[samtools_sort]="samtools sort" \
[samtools_view]="samtools view" \
[stringtie]="stringtie"
)

echo "${programs_array[*]}"
###################################################################################################

# Exit script if any command fails
set -e

# # Load Python Mox module for Python module availability
#
# module load intel-python3_2017

## Load associative array
## Only need to use one set of reads to capture sample name

# Set sample counter for array verification
sample_counter=0

# Load array
for fastq in "${fastq_dir}"*.trimmed.R1.v2.fastq.gz
do
  # Increment counter
  ((sample_counter+=1))

  # # Remove path
  # sample_name="${fastq##*/}"

  # Get sample name from first _-delimited field
  sample_name="$(basename -a ${fastq} | cut -d "." -f 1)"

  # Set treatment condition for each sample
  if [[ "${sample_name}" == "Tank_1_Crab_1" ]] \
  || [[ "${sample_name}" == "Tank_1_Crab_2" ]] \
  || [[ "${sample_name}" == "Tank_1_Crab_3" ]] \
  || [[ "${sample_name}" == "Tank_2_Crab_1" ]] \
  || [[ "${sample_name}" == "Tank_2_Crab_2" ]] \
  || [[ "${sample_name}" == "Tank_2_Crab_3" ]] \
  || [[ "${sample_name}" == "Tank_7_Crab_1" ]] \
  || [[ "${sample_name}" == "Tank_7_Crab_3" ]] \
  || [[ "${sample_name}" == "Tank_7_Crab_4" ]] \
  || [[ "${sample_name}" == "Tank_10_Crab_1" ]] \
  || [[ "${sample_name}" == "Tank_10_Crab_2" ]] \
  || [[ "${sample_name}" == "Tank_10_Crab_3" ]] \
  || [[ "${sample_name}" == "Tank_11_Crab_1" ]] \
  || [[ "${sample_name}" == "Tank_11_Crab_2" ]] \
  || [[ "${sample_name}" == "Tank_11_Crab_3" ]]
  then
    treatment="ambient"
  elif [[ "${sample_name}" == "Tank_3_Crab_1" ]] \
  || [[ "${sample_name}" == "Tank_3_Crab_2" ]] \
  || [[ "${sample_name}" == "Tank_3_Crab_3" ]] \
  || [[ "${sample_name}" == "Tank_4_Crab_1" ]] \
  || [[ "${sample_name}" == "Tank_4_Crab_2" ]] \
  || [[ "${sample_name}" == "Tank_4_Crab_3" ]] \
  || [[ "${sample_name}" == "Tank_9_Crab_1" ]] \
  || [[ "${sample_name}" == "Tank_9_Crab_2" ]] \
  || [[ "${sample_name}" == "Tank_9_Crab_3" ]] \
  || [[ "${sample_name}" == "Tank_9_Crab_4" ]] \
  || [[ "${sample_name}" == "Tank_20_Crab_1" ]] \
  || [[ "${sample_name}" == "Tank_20_Crab_2" ]] \
  || [[ "${sample_name}" == "Tank_20_Crab_3" ]]
  then
    treatment="moderate"

  else
    treatment="low"
  fi

  echo ${sample_name}
  echo ${treatment}

  # Append to associative array
  samples_associative_array+=(["${sample_name}"]="${treatment}")

done

echo "This is the final list of array indices:" "${!samples_associative_array[@]}"
echo "This is the final list of array elements:" "${samples_associative_array[@]}"

# Check array size to confirm it has all expected samples
# Exit if mismatch
if [[ "${#samples_associative_array[@]}" != "${sample_counter}" ]] \
|| [[ "${#samples_associative_array[@]}" != "${total_samples}" ]]
  then
    echo "samples_associative_array does not have all samples."
    echo ""
    echo "samples_associative_array contents:"
    echo ""
    for item in "${!samples_associative_array[@]}"
    do
      printf "%s\t%s\n" "${item}" "${samples_associative_array[${item}]}"
    done

    exit
fi

# Copy Hisat2 genome index files
rsync -av "${genome_index_dir}"/${genome_index_name}*.ht2 .

for sample in "${!samples_associative_array[@]}"
do

  ## Inititalize arrays
  fastq_array_R1=()
  fastq_array_R2=()

  # Create array of fastq R1 files
  # and generated MD5 checksums file.
  for fastq in "${fastq_dir}""${sample}"*.trimmed.R1.v2.fastq.gz
  do
    fastq_array_R1+=("${fastq}")
    echo "Generating checksum for ${fastq}..."
    md5sum "${fastq}" >> input_fastqs_checksums.md5
    echo "Checksum for ${fastq} completed."
    echo ""
  done

  # Create array of fastq R2 files
  for fastq in "${fastq_dir}""${sample}"*.trimmed.R2.v2.fastq.gz
  do
    fastq_array_R2+=("${fastq}")
    echo "Generating checksum for ${fastq}..."
    md5sum "${fastq}" >> input_fastqs_checksums.md5
    echo "Checksum for ${fastq} completed."
    echo ""
  done

  # Create comma-separated lists of FastQs for Hisat2
  printf -v joined_R1 '%s,' "${fastq_array_R1[@]}"
  fastq_list_R1=$(echo "${joined_R1%,}")

  printf -v joined_R2 '%s,' "${fastq_array_R2[@]}"
  fastq_list_R2=$(echo "${joined_R2%,}")

  # Create and switch to dedicated sample directory
  mkdir "${sample}" && cd "$_"

  # Hisat2 alignments
  # Sets read group info (RG) using samples array
  "${programs_array[hisat2]}" \
  -x "${genome_index_name}" \
  -1 "${fastq_list_R1}" \
  -2 "${fastq_list_R2}" \
  -S "${sample}".sam \
  --rg-id "${sample}" \
  --rg "SM:""${samples_associative_array[$sample]}" \
  2> "${sample}"_hisat2.err

  # Sort SAM files, convert to BAM, and index
  ${programs_array[samtools_view]} \
  -@ "${threads}" \
  -Su "${sample}".sam \
  | ${programs_array[samtools_sort]} - \
  -@ "${threads}" \
  -o "${sample}".sorted.bam
  # Index BAM
  ${programs_array[samtools_index]} "${sample}".sorted.bam


  # Run stringtie on alignments
  # Uses "-B" option to output tables intended for use in Ballgown
  # Uses "-e" option; recommended when using "-B" option.
  # Limits analysis to only reads alignments matching reference.
  "${programs_array[stringtie]}" "${sample}".sorted.bam \
  -p "${threads}" \
  -o "${sample}".gtf \
  -G "${genome_gff}" \
  -C "${sample}.cov_refs.gtf" \
  -B \
  -e

# Add GTFs to list file, only if non-empty
# Identifies GTF files that only have header
  gtf_lines=$(wc -l < "${sample}".gtf )
  if [ "${gtf_lines}" -gt 2 ]; then
    echo "$(pwd)/${sample}.gtf" >> ../"${gtf_list}"
  fi

  # Delete unneeded SAM files
  rm ./*.sam

  # Generate checksums
  for file in *
  do
    md5sum "${file}" >> ${sample}_checksums.md5
  done

  # Move up to orig. working directory
  cd ../

done

# Merge all BAMs to singular BAM for use in transcriptome assembly later
## Create list of sorted BAMs for merging
find . -name "*sorted.bam" > sorted_bams.list

## Merge sorted BAMs
${programs_array[samtools_merge]} \
-b sorted_bams.list \
${merged_bam} \
--threads ${threads}

## Index merged BAM
${programs_array[samtools_index]} ${merged_bam}



# Create singular transcript file, using GTF list file
"${programs_array[stringtie]}" --merge \
"${gtf_list}" \
-p "${threads}" \
-G "${genome_gff}" \
-o "${genome_index_name}".stringtie.gtf

# # Delete unneccessary index files
rm "${genome_index_name}"*.ht2


# Generate checksums
for file in *
do
  md5sum "${file}" >> checksums.md5
done

#######################################################################################################

# Capture program options
if [[ "${#programs_array[@]}" -gt 0 ]]; then
  echo "Logging program options..."
  for program in "${!programs_array[@]}"
  do
    {
    echo "Program options for ${program}: "
    echo ""
    # Handle samtools help menus
    if [[ "${program}" == "samtools_index" ]] \
    || [[ "${program}" == "samtools_sort" ]] \
    || [[ "${program}" == "samtools_view" ]]
    then
      ${programs_array[$program]}

    # Handle DIAMOND BLAST menu
    elif [[ "${program}" == "diamond" ]]; then
      ${programs_array[$program]} help

    # Handle NCBI BLASTx menu
    elif [[ "${program}" == "blastx" ]]; then
      ${programs_array[$program]} -help
    fi
    ${programs_array[$program]} -h
    echo ""
    echo ""
    echo "----------------------------------------------"
    echo ""
    echo ""
  } &>> program_options.log || true

    # If MultiQC is in programs_array, copy the config file to this directory.
    if [[ "${program}" == "multiqc" ]]; then
      cp --preserve ~/.multiqc_config.yaml multiqc_config.yaml
    fi
  done
  echo "Finished logging programs options."
  echo ""
fi


# Document programs in PATH (primarily for program version ID)
echo "Logging system $PATH..."
{
date
echo ""
echo "System PATH for $SLURM_JOB_ID"
echo ""
printf "%0.s-" {1..10}
echo "${PATH}" | tr : \\n
} >> system_path.log
echo "Finished logging system $PATH."
