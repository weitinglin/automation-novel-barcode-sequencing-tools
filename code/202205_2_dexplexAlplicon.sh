#!/bin/bash
#=======================================================
# author: weitinglin66
# date  : 202205
# identify amplicon barcode region
#=======================================================
cd /media/weitinglin66/new202205/analysis/202209_nanopore_sl/02_filtered_fastq

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

mismatch_parameter=7

outputpath=/home/weitinglin66/Documents/analysis/202205_2_nanopore/202205_2_amplicon

for item in $(seq 0 52)
do
echo File $item 
# barcode 1
seqkit amplicon -F CGCTGATTGAAGCAGATACCT    -R GAATGTGAAATAGGATTTCCGAT -m $mismatch_parameter -r 22:-24  FAR28891_pass_766e7ef0_${item}_filtered.fasta |\
bioawk -c fastx '{print $name,$seq,"barcode1"}' >  $outputpath/FAR28891_pass_766e7ef0_${item}_barcode.txt
# barcode 2
seqkit amplicon -F GATATCGTTATTAATACAACACC  -R GGATTAAGATTGCAGTTCCTAG  -m $mismatch_parameter -r 24:-23  FAR28891_pass_766e7ef0_${item}_filtered.fasta |\
bioawk -c fastx '{print $name,$seq,"barcode2"}' >> $outputpath/FAR28891_pass_766e7ef0_${item}_barcode.txt
# barcode 3  
seqkit amplicon -F GAATTGGACATGCGATTTTACCA  -R GATGGTCCCGATAAAACGAT    -m $mismatch_parameter -r 24:-21  FAR28891_pass_766e7ef0_${item}_filtered.fasta |\
bioawk -c fastx '{print $name,$seq,"barcode3"}' >> $outputpath/FAR28891_pass_766e7ef0_${item}_barcode.txt
# barcode 4
seqkit amplicon -F TACGATTACGTTGTAGTTAATGA  -R CCTTCAGGTAATACGATTTTAAC -m $mismatch_parameter -r 24:-24  FAR28891_pass_766e7ef0_${item}_filtered.fasta |\
bioawk -c fastx '{print $name,$seq,"barcode4"}' >> $outputpath/FAR28891_pass_766e7ef0_${item}_barcode.txt
# barcode 5
seqkit amplicon -F TTAAGCTTTTCAACAAAAGGGTC  -R TTCACGACGTTCAGAATGAACGA -m $mismatch_parameter -r 24:-24  FAR28891_pass_766e7ef0_${item}_filtered.fasta |\
bioawk -c fastx '{print $name,$seq,"barcode5"}' >> $outputpath/FAR28891_pass_766e7ef0_${item}_barcode.txt
# barcode 6
seqkit amplicon -F GTACAATTGTTAGAAGGTGCAAA  -R GCCAATAGGTGTCCTGTATGCTG -m $mismatch_parameter -r 24:-24  FAR28891_pass_766e7ef0_${item}_filtered.fasta |\
bioawk -c fastx '{print $name,$seq,"barcode6"}' >> $outputpath/FAR28891_pass_766e7ef0_${item}_barcode.txt



done

