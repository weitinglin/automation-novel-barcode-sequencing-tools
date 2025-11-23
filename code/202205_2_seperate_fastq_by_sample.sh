#!/bin/bash
#=======================================================
# author: weitinglin66
# log date  : 202205
# purpose: use demultiplexed sample ID to seperate reads in fastq files by sample
# input : fastq files
# required data : read ID by sample after demultiplexing (which generated from R script)
# output : sample fastq files grouping by read ID
#=======================================================


fastqpath=/home/system76/extradrive/storage_weitinglin66/rawdata/20220510_oxfordnanopore/20220510_secondRecaling/fastq_pass
sampleIDpath=/home/weitinglin66/Documents/analysis/202205_2_nanopore/202205_2_sampleID
outputpath=/home/weitinglin66/Documents/analysis/202205_2_nanopore/202205_2_seperateFastq

for file in $(seq 0 52)
do
echo ============ sample $file =================
    for item in $(seq 1 10)
    do
echo ====== subsettting File $item ============= 

seqkit grep -f $sampleIDpath/20220515_2_${file}_sample${item}_SeqID.txt -o $outputpath/fastq_runid_766e7ef0b15421ebbf79d5a0aef6f1e60aaeb327_${file}_sample${item}.fastq  $fastqpath/fastq_runid_766e7ef0b15421ebbf79d5a0aef6f1e60aaeb327_${file}_0.fastq.gz

    done
done
