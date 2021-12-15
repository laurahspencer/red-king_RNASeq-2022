#!/bin/bash

BASE=/home/ggoetz/rairay/nichols/201910-redking_crab-rnaseq

IN=${BASE}/concat

FILES=$(ls ${IN}/*.fastq | awk -F "/" '{print $NF}')

for file in ${FILES}
do
    echo ${file}
    gzip -c ${IN}/${file} > ${IN}/${file}.gz
done
