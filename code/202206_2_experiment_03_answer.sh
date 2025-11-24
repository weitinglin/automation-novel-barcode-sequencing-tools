#!/bin/bash
#=======================================================
# author: weitinglin66
# date  : 20220612
# Purpose : assemble the result contig from merged fastq files by sample
#=======================================================
#
#
mergefastq_path=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_mergeFastqbySample
experiment_path=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202206_2_experiment
assembly_path=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202206_2_experiment/sample4_assembly
result_path=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202206_2_result


for sample in $(seq 1 10)
do

echo =======================================================================
echo Build directory =======================================================
echo =======================================================================
cd $result_path
    mkdir sample${sample}_result

echo =======================================================================
echo                                                                       =
echo  Sample $sample assembly
echo                                                                       =
echo =======================================================================
cd $result_path/sample${sample}_result
## asemble
    flye --nano-raw\
	    $mergefastq_path/merged_20220612_sample${sample}.fastq\
        --threads 8\
        --min-overlap 1000\
        --asm-coverage 40\
        --genome-size 4k\
	    --out-dir $result_path/sample${sample}_result
    mv assembly.fasta sample${sample}_assembly.fasta
    bwa index sample${sample}_assembly.fasta
    bwa mem sample${sample}_assembly.fasta $mergefastq_path/merged_20220612_sample${sample}.fastq| samtools sort -o sample${sample}_assembly_output.bam
    samtools index sample${sample}_assembly_output.bam
# polish
    medaka consensus sample${sample}_assembly_output.bam sample${sample}_assembly_polish.hdf
    medaka stitch sample${sample}_assembly_polish.hdf sample${sample}_assembly.fasta sample${sample}_assembly_polish_consensus.fasta
    bwa index sample${sample}_assembly_polish_consensus.fasta
    bwa mem sample${sample}_assembly_polish_consensus.fasta $mergefastq_path/merged_20220612_sample${sample}.fastq| samtools sort -o sample${sample}_assembly_polish_output.bam
    samtools index sample${sample}_assembly_polish_output.bam

done

    
