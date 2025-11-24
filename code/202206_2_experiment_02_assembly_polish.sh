#!/bin/bash
#=======================================================
# author: weitinglin66
# date  : 20220612
# Experiment different reads and the possible result to calculate the possibility
#=======================================================
#
#
mergefastq_path=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_mergeFastqbySample
experiment_path=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202206_2_experiment
assembly_path=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202206_2_experiment/sample4_assembly
# select


for reads in  $1
do
mkdir $assembly_path/sample4_${reads}_assembly

for num in $(seq 1 100)
do
echo =======================================================================
echo Build directory =======================================================
echo =======================================================================

cd $assembly_path/sample4_${reads}_assembly
    mkdir sample4_${reads}_assembly_${num}

echo =======================================================================
echo                                                                       =
echo  Sample4 reads $reads in  assfile $num are assembly ===========================
echo                                                                       =
echo =======================================================================
cd $assembly_path/sample4_${reads}_assembly/sample4_${reads}_assembly_${num}

## asemble
    flye --nano-raw\
	    $experiment_path/sample4_fastq/sample4_${reads}/sample4_${reads}_${num}.fastq\
        --threads 16\
        --min-overlap 1000\
        --asm-coverage 40\
        --genome-size 4k\
	    --out-dir $assembly_path/sample4_${reads}_assembly/sample4_${reads}_assembly_${num}
    mv assembly.fasta sample4_${reads}_assembly_${num}.fasta
    bwa index sample4_${reads}_assembly_${num}.fasta
    bwa mem sample4_${reads}_assembly_${num}.fasta $experiment_path/sample4_fastq/sample4_${reads}/sample4_${reads}_${num}.fastq| samtools sort -o sample4_${reads}_assembly_${num}_output.bam
    samtools index sample4_${reads}_assembly_${num}_output.bam
# polish
    medaka consensus sample4_${reads}_assembly_${num}_output.bam sample4_${reads}_assembly_${num}_polish.hdf
    medaka stitch sample4_${reads}_assembly_${num}_polish.hdf sample4_${reads}_assembly_${num}.fasta sample4_${reads}_assembly_${num}_polish_consensus.fasta
    bwa index sample4_${reads}_assembly_${num}_polish_consensus.fasta
    bwa mem sample4_${reads}_assembly_${num}_polish_consensus.fasta $experiment_path/sample4_fastq/sample4_${reads}/sample4_${reads}_${num}.fastq| samtools sort -o sample4_${reads}_assembly_${num}_polish_output.bam
    samtools index sample4_${reads}_assembly_${num}_polish_output.bam
done

done

    
