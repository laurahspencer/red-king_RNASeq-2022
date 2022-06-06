#!/bin/bash

#SBATCH --job-name=bayescan
#SBATCH --output=/home/lspencer/2022-redking-OA/sbatch_logs/bayescan_miss25.txt
#SBATCH --mail-user=laura.spencer@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH -c 20
#SBATCH -t 3-0:0:0

# Load modules
module load bio/pgdspider/2.1.1.5
module load bio/bayescan/2.1

# Note- to access PGDSpider:
  # The env variable PGDSPIDER is set to the path of the jar file
  # java -jar $PGDSPIDER <pgdspider commands>

DIR=/scratch/lspencer/2022-redking-OA/genotypes

cd ${DIR}
pwd

# Run PGDSpider to convert VCF to Bayescan format

# Using pH treatment as population id
java -jar $PGDSPIDER \
-inputfile ${DIR}/rkc_rnaseq_genotypes-final_miss25.recode.vcf -inputformat VCF \
-outputfile ${DIR}/rkc_miss25_for_bayescan_pH -outputformat GESTE_BAYE_SCAN \
-spid ${DIR}/pops-pH.spid

# Using genetic clustering as population id
java -jar $PGDSPIDER \
-inputfile ${DIR}/rkc_rnaseq_genotypes-final_miss25.recode.vcf -inputformat VCF \
-outputfile ${DIR}/rkc_miss25_for_bayescan_cluster -outputformat GESTE_BAYE_SCAN \
-spid ${DIR}/pops-cluster.spid

# Run bayescan to identify outlier loci using pH treatment as grouping variable
bayescan rkc_miss25_for_bayescan_pH \
-snp -threads 20 -n 100000 -pr_odds 10 -o rkc_miss25_bayescan_outliers_pH

# Run bayescan to identify outlier loci using pH treatment as grouping variable
bayescan rkc_miss25_for_bayescan_cluster \
-snp -threads 20 -n 100000 -pr_odds 10 -o rkc_miss25_bayescan_outliers_cluster
