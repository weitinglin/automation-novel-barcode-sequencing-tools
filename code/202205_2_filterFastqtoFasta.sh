#!/bin/bash
#=======================================================
# author: weitinglin66
# log date  : 20220502
# purpose: filter the fastq file according to length threshold between 3k-4kb and convert to fasta file
# input : raw fastq file
# required data : length threshold
# output : filtered fasta file
#=======================================================

fastq_suffix=runid_766e7ef0b15421ebbf79d5a0aef6f1e60aaeb327
outputpath=/your/put/path/filtered_fasta
outputfile_suffix=FAR28891_pass_766e7ef0

for item in $(seq 0 52)
do
echo File $item 
bioawk -c fastx '{print ">"$name; print $seq}' ${fastq_suffix}_${item}_0.fastq.gz|\
bioawk -c fastx 'length($seq)>3000 {print ">"$name;print $seq}'|\
bioawk -c fastx 'length($seq)<4000 {print ">"$name;print $seq}'> $outputpath/${outputfile_suffix}_${item}_filtered.fasta
done

