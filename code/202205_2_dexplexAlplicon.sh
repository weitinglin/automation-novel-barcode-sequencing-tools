#!/bin/bash
#=======================================================
# author : weitinglin66
# log date : 202205
# purpose : identify amplicon barcode region
# use the primers sequence to extract the amplicon region with seqkit amplicon, which as barcode information
#
# input : filtered fasta file
# required data : primer sequences, mismatch parameter
# output : amplicon barcode sequence file
#
#=======================================================

## primer sequences ##

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

## mismatch parameter ##
mismatch_parameter=7
## output path ##
outputpath=/your/put/path/amplicon_barcode
fastaqpath=/your/put/path/filtered_fasta
fastafile_suffix=FAR28891_pass_766e7ef0

## main code ##

for item in $(seq 0 52)
do
echo File $item 
# barcode 1
seqkit amplicon -F CGCTGATTGAAGCAGATACCT    -R GAATGTGAAATAGGATTTCCGAT -m $mismatch_parameter -r 22:-24  ${fastafile_suffix}_${item}_filtered.fasta |\
bioawk -c fastx '{print $name,$seq,"barcode1"}' >  $outputpath/${fastafile_suffix}_${item}_barcode.txt
# barcode 2
seqkit amplicon -F GATATCGTTATTAATACAACACC  -R GGATTAAGATTGCAGTTCCTAG  -m $mismatch_parameter -r 24:-23  ${fastafile_suffix}_${item}_filtered.fasta |\
bioawk -c fastx '{print $name,$seq,"barcode2"}' >> $outputpath/${fastafile_suffix}_${item}_barcode.txt
# barcode 3  
seqkit amplicon -F GAATTGGACATGCGATTTTACCA  -R GATGGTCCCGATAAAACGAT    -m $mismatch_parameter -r 24:-21  ${fastafile_suffix}_${item}_filtered.fasta |\
bioawk -c fastx '{print $name,$seq,"barcode3"}' >> $outputpath/${fastafile_suffix}_${item}_barcode.txt
# barcode 4
seqkit amplicon -F TACGATTACGTTGTAGTTAATGA  -R CCTTCAGGTAATACGATTTTAAC -m $mismatch_parameter -r 24:-24  ${fastafile_suffix}_${item}_filtered.fasta |\
bioawk -c fastx '{print $name,$seq,"barcode4"}' >> $outputpath/${fastafile_suffix}_${item}_barcode.txt
# barcode 5
seqkit amplicon -F TTAAGCTTTTCAACAAAAGGGTC  -R TTCACGACGTTCAGAATGAACGA -m $mismatch_parameter -r 24:-24  ${fastafile_suffix}_${item}_filtered.fasta |\
bioawk -c fastx '{print $name,$seq,"barcode5"}' >> $outputpath/${fastafile_suffix}_${item}_barcode.txt
# barcode 6
seqkit amplicon -F GTACAATTGTTAGAAGGTGCAAA  -R GCCAATAGGTGTCCTGTATGCTG -m $mismatch_parameter -r 24:-24  ${fastafile_suffix}_${item}_filtered.fasta |\
bioawk -c fastx '{print $name,$seq,"barcode6"}' >> $outputpath/${fastafile_suffix}_${item}_barcode.txt



done

