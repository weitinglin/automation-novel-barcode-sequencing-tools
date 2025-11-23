#!/bin/bash
#=======================================================
# author: weitinglin66
# log date  : 202205
# purpose : merge fastq files by sample ID into individual merged fastq file
# input : fastq files by sample ID
# required data : path information
# output : merged fastq files by sample
#=======================================================

fastqpath=/your/input/path/fastq_by_sample
fastq_suffix=runid_766e7ef0b15421ebbf79d5a0aef6f1e60aaeb327

outputpath=/your/output/path/merged_fastq_by_sample
outputfile_suffix=merged_20220612

cd $fastqpath
for item in $(seq 1 10)
do
echo ====== concate into sample $item ============= 
# barcode 1
cat ${fastq_suffix}_*_sample${item}.fastq > ${outputpath}/${outputfile_suffix}_sample${item}.fastq
done

