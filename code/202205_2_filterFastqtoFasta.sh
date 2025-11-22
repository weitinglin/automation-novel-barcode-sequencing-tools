#!/bin/bash
#=======================================================
# author: weitinglin66
# date  : 202205_2
# filter and turn fastq to fasta
#=======================================================
cd /home/system76/extradrive/storage_weitinglin66/rawdata/20220510_oxfordnanopore/20220510_secondRecaling/fastq_pass



for item in $(seq 0 52)
do
echo File $item 
bioawk -c fastx '{print ">"$name; print $seq}' fastq_runid_766e7ef0b15421ebbf79d5a0aef6f1e60aaeb327_${item}_0.fastq.gz|\
bioawk -c fastx 'length($seq)>3000 {print ">"$name;print $seq}'|\
bioawk -c fastx 'length($seq)<4000 {print ">"$name;print $seq}'> /home/weitinglin66/Documents/analysis/202205_2_nanopore/202205_2_filterfasta/FAR28891_pass_766e7ef0_${item}_filtered.fasta
done

