#!/bin/bash
#=======================================================
# author: weitinglin66
# date  : 20220515
# Assemlby pipeline testing
#=======================================================
#
# step 1: reads filter
# step 2: assembly with flye
#

sample_fastqPath=/home/weitinglin66/Documents/analysis/202205_2_nanopore/202205_2_seperateFastq

cd /home/weitinglin66/Documents/analysis/202205_nanopore
#mkdir 20220515_assembly

for file in $(seq 0 52)
do

# $(seq 0 28)
for item in $(seq 1 10)
do

echo Build directory =======================================================
cd /home/weitinglin66/Documents/analysis/202205_2_nanopore/202205_2_assemlby
    mkdir 20220515_${file}_sample_${item}

echo  Sample $item in  assfile $file are assembly ===========================

cd /home/weitinglin66/Documents/analysis/202205_2_nanopore/202205_2_assemlby/20220515_${file}_sample_${item}
    flye --nano-raw\
	    ${sample_fastqPath}/fastq_runid_766e7ef0b15421ebbf79d5a0aef6f1e60aaeb327_${file}_sample${item}.fastq\
        --threads 8\
        --min-overlap 1000\
        --asm-coverage 40\
        --genome-size 4k\
	    --out-dir /home/weitinglin66/Documents/analysis/202205_2_nanopore/202205_2_assemlby/20220515_${file}_sample_${item}
    mv assembly.fasta assembly_${file}_sample${item}.fasta
    bwa index assembly_${file}_sample${item}.fasta
    bwa mem assembly_${file}_sample${item}.fasta ${sample_fastqPath}/fastq_runid_766e7ef0b15421ebbf79d5a0aef6f1e60aaeb327_${file}_sample${item}.fastq| samtools sort -o 20220515_${file}_${item}_output.bam
    samtools index 20220515_${file}_${item}_output.bam
done

done

    
