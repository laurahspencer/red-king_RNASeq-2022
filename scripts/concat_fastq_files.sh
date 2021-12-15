#!/bin/bash

BASE=/home/ggoetz/rairay/nichols/201910-redking_crab-rnaseq

IN=${BASE}/raw
OUT=${BASE}/concat

FOLDERS=$(ls -1d ${IN}/3* | awk -F "/" '{ print $NF }')

for folder in ${FOLDERS}
do
    echo ${folder}
    SUB_FOLDER=${IN}/${folder}
    SAMPLES=$(ls ${SUB_FOLDER}/3*_R1_*.fastq.gz | \
        awk -F "/" '{ print $NF }' | \
        awk -F "_" '{ print "Tank_"$3"_Crab_"$5 }')

    for sample in ${SAMPLES}
    do
        echo ${sample}
        zcat ${SUB_FOLDER}/${folder}_${sample}_*_R1_*.fastq.gz \
            >> ${OUT}/${sample}_R1.fastq
        zcat ${SUB_FOLDER}/${folder}_${sample}_*_R2_*.fastq.gz \
            >> ${OUT}/${sample}_R2.fastq
    done
done
