#!/bin/bash
#=======================================================
# author: weitinglin66
# date  : 202205
# copy files from Nanopore
#=======================================================
cd /home/weitinglin66/Documents/analysis/202205_2_nanopore/202205_2_filterfasta

fastqpath=/home/system76/extradrive/storage_weitinglin66/rawdata/20220510_oxfordnanopore/20220510_secondRecaling/fastq_pass
sampleIDpath=/home/weitinglin66/Documents/analysis/202205_2_nanopore/202205_2_sampleID
outputpath=/home/weitinglin66/Documents/analysis/202205_2_nanopore/202205_2_seperateFastq

for file in $(seq 0 52)
do
echo ============ sample $file =================
    for item in $(seq 1 10)
    do
echo ====== subsettting File $item ============= 
# barcode 1
seqkit grep -f $sampleIDpath/20220515_2_${file}_sample${item}_SeqID.txt -o $outputpath/fastq_runid_766e7ef0b15421ebbf79d5a0aef6f1e60aaeb327_${file}_sample${item}.fastq  $fastqpath/fastq_runid_766e7ef0b15421ebbf79d5a0aef6f1e60aaeb327_${file}_0.fastq.gz
#seqkit grep -f $sampleIDpath/20220515_1_sample${item}_SeqID.txt -o $outputpath/FAR28891_pass_1_sample${item}.fastq  $fastqpath/FAR28891_pass_766e7ef0_1.fastq.gz
#seqkit grep -f $sampleIDpath/20220515_1_sample${item}_SeqID.txt -o $outputpath/FAR28891_pass_1_sample$item.fastq  FAR28891_pass_766e7ef0_1.fasta 
#seqkit grep -p ff323748-2d66-4108-a0db-c6025a03c2a9 -o $outputpath/FAR28891_pass_1_sample$item.fastq  FAR28891_pass_766e7ef0_1.fasta 
    done
done
