#!/bin/bash
#=======================================================
# author: weitinglin66
# date  : 20220612
# Build the local MLST database
#=======================================================
#
# step 1: reads filter
# step 2: assembly with flye
#
MLST_datapath=/media/weitinglin66/new202205/202205_MLST
#sample_assembly_path=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202206_2_experiment/sample4_assembly/sample4_10_assembly/sample4_10_assembly_1
mergefastq_path=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_mergeFastqbySample
experiment_path=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202206_2_experiment
assembly_path=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202206_2_experiment/sample4_assembly
genotyperesult_path=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202206_2_experiment/sample4_genotype
result_path=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202206_2_result

# arcC_20220520.fas  glpF_20220520.fas  pta_20220520.fas  yqiL_20220520.fas
# aroE_20220520.fas  gmk_20220513.fas   tpi_20220520.fas

# bioawk -v var="$item" -c fastx '{print ">item"var;print $seq}' >>  $outputpath/sample${sample}/sample${sample}_${MLSTgenes[${gene}]}.fasta

## with sequence outpu
#seqkit fish -a -f $MLST_datapath/arcC_20220520.fas -g $sample_assembly_path/sample4_10_assembly_1_polish_consensus.fasta -b $MLST_datapath/sample4_10_1_arc.bam 
## without sequence output
#item=sample4_10_1
#cat $sample_assembly_path/sample4_10_assembly_1_polish_consensus.fasta | bioawk -v var="$item" -c fastx '{print ">"var"_contig";print $seq}' |\
#seqkit fish -a -f $MLST_datapath/arcC_20220520.fas -b $MLST_datapath/sample4_10_1_arc.bam 

# use minimap2, so fast!!!!!!!!!!!!!!!
#minimap2 -ax sr $sample_assembly_path/sample4_10_assembly_1_polish_consensus.fasta $MLST_datapath/arcC_20220520.fas 

#minimap2 -a -x sr $sample_assembly_path/sample4_10_assembly_1_polish_consensus.fasta $MLST_datapath/arcC_20220520.fas | samtools view -d NM:0 | cut -f 1,2,3,4,5

#minimap2 -a -x sr $sample_assembly_path/sample4_10_assembly_1_polish_consensus.fasta $MLST_datapath/all_MLST_20220612.fas |\
# samtools view -e '[NM]<1' |\
#  cut -f 1,3,4,5,12,13,14,15,16

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
