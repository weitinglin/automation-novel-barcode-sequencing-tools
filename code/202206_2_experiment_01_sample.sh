#!/bin/bash
#=======================================================
# author: weitinglin66
# date  : 20220612
# purpose : randomly sample reads from merged fastq files by sample for experiment to see the effect of different reads number on assembly result 
#=======================================================
#
#
mergefastq_path=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_mergeFastqbySample
experiment_path=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202206_2_experiment

# select
for reads in 51 52 53 54 55 56 57 58 59 60
do
mkdir $experiment_path/sample4_fastq/sample4_${reads}
    for num in $(seq 1 100)
    do
        seqkit sample -n $reads -s $num -o $experiment_path/sample4_fastq/sample4_${reads}/sample4_${reads}_${num}.fastq $mergefastq_path/merged_20220612_sample4.fastq
    done
done
