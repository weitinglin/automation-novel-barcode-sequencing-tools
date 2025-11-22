#!/bin/bash
#=======================================================
# author: weitinglin66
# date  : 202206
# use mekada to polish the draft
#=======================================================

draftfasta_path=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_assemlby
output_path=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_polish
sepfastq_path=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_seperateFastq

cd /media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_polish

for item in 21 22 23 24
do
medaka consensus $draftfasta_path/20220515_${item}_sample_1/20220515_${item}_1_output.bam $output_path/sample1_${item}_output.hdf
medaka stitch $output_path/sample1_${item}_output.hdf $draftfasta_path/20220515_${item}_sample_1/assembly_${item}_sample1.fasta $output_path/sample1_${item}_consensus.fasta
bwa index $output_path/sample1_${item}_consensus.fasta
bwa mem $output_path/sample1_${item}_consensus.fasta ${sepfastq_path}/fastq_runid_766e7ef0b15421ebbf79d5a0aef6f1e60aaeb327_${item}_sample1.fastq| samtools sort -o ${output_path}/file${item}_sample1_polish_output.bam
samtools index ${output_path}/file${item}_sample1_polish_output.bam
done