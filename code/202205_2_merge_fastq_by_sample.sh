#!/bin/bash
#=======================================================
# author: weitinglin66
# date  : 202205
# copy files from Nanopore
#=======================================================
cd /home/weitinglin66/Documents/analysis/202205_2_nanopore/202205_2_filterfasta

fastqpath=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_seperateFastq
outputpath=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_mergeFastqbySample

cd $fastqpath
for item in $(seq 1 10)
do
echo ====== concate into sample $item ============= 
# barcode 1
cat fastq_runid_766e7ef0b15421ebbf79d5a0aef6f1e60aaeb327_*_sample${item}.fastq > ${outputpath}/merged_20220612_sample${item}.fastq
done

