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

cd $MLST_datapath
# arcC_20220520.fas  glpF_20220520.fas  pta_20220520.fas  yqiL_20220520.fas
# aroE_20220520.fas  gmk_20220513.fas   tpi_20220520.fas

# bioawk -v var="$item" -c fastx '{print ">item"var;print $seq}' >>  $outputpath/sample${sample}/sample${sample}_${MLSTgenes[${gene}]}.fasta

for gene in arcC glpF pta yqiL aroE gmk tpi
do

number=$(bioawk -c fastx '{print $name}' ${gene}_202205*.fas | wc -l)
printf "%-10s %-8s \n" $gene $number
done