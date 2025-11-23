#!/bin/bash
#=======================================================
# author: weitinglin66
# log date  : 202205
# purpose : copy files from Nanopore
# input : fasta files
# required data : path information
# output : merged fasta files by sample
#=======================================================

fastqpath=/home/system76/extradrive/storage_weitinglin66/rawdata/20220510_oxfordnanopore/20220510_secondRecaling/fastq_pass
sampleIDpath=/home/weitinglin66/Documents/analysis/202205_2_nanopore/202205_2_sampleID
outputpath=/home/weitinglin66/Documents/analysis/202205_2_nanopore/202205_2_seperateFastq


cat /home/weitinglin66/Documents/analysis/202205_2_nanopore/202205_2_assemlby/20220515_20_sample_1/assembly_20_sample1.fasta |\
 bioawk -c fastx '{print ">""20_sample1";print $seq}' > /home/weitinglin66/Documents/analysis/202205_2_nanopore/202205_2_msaFasta/msa_sample_1.fasta

for item in $(seq 21 52)
do

cat /home/weitinglin66/Documents/analysis/202205_2_nanopore/202205_2_assemlby/20220515_${item}_sample_1/assembly_${item}_sample1.fasta|\
 bioawk -v var="$item" -c fastx '{print ">"var"_sample1";print $seq}' >> /home/weitinglin66/Documents/analysis/202205_2_nanopore/202205_2_msaFasta/msa_sample_1.fasta


done