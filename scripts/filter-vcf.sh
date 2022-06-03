#!/bin/bash

#SBATCH --job-name=filter-vcf
#SBATCH --output=/home/lspencer/2022-redking-OA/sbatch_logs/filter-vcf.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 3-0:0:0

module load bio/vcftools/0.1.16

vcftools --gzvcf \
"/scratch/lspencer/2022-redking-OA/genotypes/rkc_rnaseq_genotypes-filtered-true.vcf.gz" \
--max-missing 0.75 --maf 0.05 --recode --recode-INFO-all --out \
"/scratch/lspencer/2022-redking-OA/genotypes/rkc_rnaseq_genotypes-final_miss25"

vcftools --gzvcf \
"/scratch/lspencer/2022-redking-OA/genotypes/rkc_rnaseq_genotypes-filtered-true.vcf.gz" \
--max-missing 0.80 --maf 0.05 --recode --recode-INFO-all --out \
"/scratch/lspencer/2022-redking-OA/genotypes/rkc_rnaseq_genotypes-final_miss20"

vcftools --gzvcf \
"/scratch/lspencer/2022-redking-OA/genotypes/rkc_rnaseq_genotypes-filtered-true.vcf.gz" \
--max-missing 0.85 --maf 0.05 --recode --recode-INFO-all --out \
"/scratch/lspencer/2022-redking-OA/genotypes/rkc_rnaseq_genotypes-final_miss15"
