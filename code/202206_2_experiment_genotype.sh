#!/bin/bash
#=======================================================
# author: weitinglin66
# date  : 20220612
# purpose : making visulization bam files for genotype result on IGV 
#=======================================================
#
#
MLST_datapath=/media/weitinglin66/new202205/202205_MLST
#sample_assembly_path=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202206_2_experiment/sample4_assembly/sample4_10_assembly/sample4_10_assembly_1
mergefastq_path=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_mergeFastqbySample
experiment_path=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202206_2_experiment
assembly_path=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202206_2_experiment/sample4_assembly
genotyperesult_path=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202206_2_experiment/sample4_genotype
result_path=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202206_2_result



echo > $result_path/sample4_25_genotype_result.txt

for reads in 25
do
echo ==============================================================
echo ========= reads $reads test ================================== 
echo ==============================================================
    for sample in 1
    do
echo ==============================================================
echo ========= reads $reads test  in  $sample =====================
echo ==============================================================
#    minimap2 -a -x sr $assembly_path/sample4_${reads}_assembly/sample4_${reads}_assembly_${sample}/sample4_${reads}_assembly_${sample}_polish_consensus.fasta $MLST_datapath/all_MLST_20220612.fas | samtools view -d NM:0 |\
#    cut -f 1,3,4,5,6,12,13,14,15,16| awk  -v var="sample4_${reads}_$sample" '{printf "%-s\t%-s\t%-s\t%-s\t%-s\t%-s\t%-s\t%-s\n",$1, var, $5, $6, $7, $8, $9, $10}' >> $result_path/sample4_25_genotype_result.txt
minimap2 -a -x sr $assembly_path/sample4_${reads}_assembly/sample4_${reads}_assembly_${sample}/sample4_${reads}_assembly_${sample}_polish_consensus.fasta $MLST_datapath/all_MLST_20220612.fas | samtools view -S -b | samtools sort >  $result_path/sample4_${reads}_${sample}_MLSA_type.bam
samtools index $result_path/sample4_${reads}_${sample}_MLSA_type.bam
    done
done
