#!/bin/bash
#=======================================================
# author: weitinglin66
# date  : 20220529
# resplit the fasta, isolate the region we want
#=======================================================
cd /media/weitinglin66/new202205/analysis/202205_nanopore/20220528_sample_MLST

declare -a PrimerSeq_gap_front
PrimerSeq_gap_front[1]='CGCTGATTGAAGCAGATACCT'
PrimerSeq_gap_front[2]='GATATCGTTATTAATACAACACC' #CATGGTGCTTTTTCAATGGAGAC 
PrimerSeq_gap_front[3]='GAATTGGACATGCGATTTTACCA'
PrimerSeq_gap_front[4]='TACGATTACGTTGTAGTTAATGA'
PrimerSeq_gap_front[5]='TTAAGCTTTTCAACAAAAGGGTC'
PrimerSeq_gap_front[6]='GTACAATTGTTAGAAGGTGCAAA'
declare -a PrimerSeq_gap_back
PrimerSeq_gap_back[1]='GAATGTGAAATAGGATTTCCGAT'
PrimerSeq_gap_back[2]='GGATTAAGATTGCAGTTCCTAG' #GGAGATTTCTACGAGCCAAA
PrimerSeq_gap_back[3]='GATGGTCCCGATAAAACGAT'
PrimerSeq_gap_back[4]='CCTTCAGGTAATACGATTTTAAC'
PrimerSeq_gap_back[5]='TTCACGACGTTCAGAATGAACGA'
PrimerSeq_gap_back[6]='GCCAATAGGTGTCCTGTATGCTG'

declare -a MLSTprimer_forward
MLSTprimer_forward[1]='TTGATTCACCAGCGCGTATTGTC' #arc_F
MLSTprimer_forward[2]='ATCGGAAATCCTATTTCACATTC' #aro_F
MLSTprimer_forward[3]='CTAGGAACTGCAATCTTAATCC'  #glp_F
MLSTprimer_forward[4]='ATCGTTTTATCGGGACCATC'    #gmk_F
MLSTprimer_forward[5]='GTTAAAATCGTATTACCTGAAGG' #pta_F
MLSTprimer_forward[6]='TCGTTCATTCTGAACGTCGTGAA' #tpi_F
MLSTprimer_forward[7]='CAGCATACAGGACACCTATTGGC' #yqi_F
declare -a MLSTprimer_reverse
MLSTprimer_reverse[1]='AGGTATCTGCTTCAATCAGCG'   #arc_R
MLSTprimer_reverse[2]='GGTGTTGTATTAATAACGATATC' #aro_R
MLSTprimer_reverse[3]='TGGTAAAATCGCATGTCCAATTC' #glp_R
MLSTprimer_reverse[4]='TCATTAACTACAACGTAATCGTA' #gmk_R
MLSTprimer_reverse[5]='GACCCTTTTGTTGAAAAGCTTAA' #pta_R
MLSTprimer_reverse[6]='TTTGCACCTTCTAACAATTGTAC' #tpi_R
MLSTprimer_reverse[7]='CGTTGAGGAATCGATACTGGAAC' #yqi_R
mismatch_parameter=7

declare -a MLSTgenes
MLSTgenes[1]='arc'
MLSTgenes[2]='aro'
MLSTgenes[3]='glp'
MLSTgenes[4]='gmk'
MLSTgenes[5]='pta'
MLSTgenes[6]='tpi'
MLSTgenes[7]='yqi'


assembly_result=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_assemlby
outputpath=/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_sample_MLST

# Make the Directory
#for sample in $(seq 1 10)
#do

#mkdir $outputpath/sample${sample}
#cd $outputpath/sample${sample}
#    for gene in $(seq 1 7)
#    do
#    echo > sample${sample}_${MLSTgenes[${gene}]}.fasta
#
#    done
#
#
#done



for item in $(seq 20 52)
do
echo Processing File number: $item =======================

 for sample in $(seq 1 10)
    do
    echo Processing sample $sample ===========================
    cd $outputpath/sample${sample}
     for gene in $(seq 1 7)
     do
        echo Processing ${MLSTgenes[${gene}]} =========================
echo seqkit amplicon -F ${MLSTprimer_forward[${gene}]} -R ${MLSTprimer_reverse[${gene}]} -m $mismatch_parameter $assembly_result/20220515_${item}_sample_${sample}/assembly_${item}_sample${sample}.fasta 
seqkit amplicon -F ${MLSTprimer_forward[${gene}]} -R ${MLSTprimer_reverse[${gene}]} -m $mismatch_parameter $assembly_result/20220515_${item}_sample_${sample}/assembly_${item}_sample${sample}.fasta |\
bioawk -v var="$item" -c fastx '{print ">item"var;print $seq}' >>  $outputpath/sample${sample}/sample${sample}_${MLSTgenes[${gene}]}.fasta

     done
 done
done

#bioawk -v var="$item" -c fastx '{print ">"var"_sample1";print $seq}' >> /home/weitinglin66/Documents/analysis/202205_2_nanopore/202205_2_msaFasta/msa_sample_1.fasta