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
output_path=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202206_2_experiment/sample4_check
# select
for reads in 11
do
echo > $output_path/sample4_11_fastqcheck.txt
    for num in $(seq 1 100)
    do
    echo ================== $reads $num =======================
        bioawk -c fastx -v var="sample4_${reads}_$num" '{print $name"\t"var"\t"length($seq)}' $experiment_path/sample4_alt_fastq/sample4_${reads}/sample4_${reads}_${num}.fastq  >> $output_path/sample4_11_fastqcheck.txt
    done
done
