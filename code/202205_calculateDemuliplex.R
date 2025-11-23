
#=======================================================
# author : weitinglin66
# log date : 202205
# purpose : use the deplexAmplicon results to calculate the demultiplexing result and individual reads barcode combination, asign to each sample
#
# input : amplicon barcode sequence file from deplexAmplicon
# required data : primer sequences, mismatch parameter
# output : each read's sample ID assignment file
#
#=======================================================


# Loading library ---------------------------------------------------------
library(tidyverse)
library(Biostrings)
library(ggmsa)
library(msa)

# Primer sequences -------------------------------------------------------
ymarc_R <- 'CGCTGATTGAAGCAGATACCT'
aro_R <- 'GATATCGTTATTAATACAACACC' #CATGGTGCTTTTTCAATGGAGAC 
glp_R <- 'GAATTGGACATGCGATTTTACCA'
gmk_R <- 'TACGATTACGTTGTAGTTAATGA'
pta_R <- 'TTAAGCTTTTCAACAAAAGGGTC'
tpi_R <- 'GTACAATTGTTAGAAGGTGCAAA'
aro_F <- 'GAATGTGAAATAGGATTTCCGAT'
glp_F <- 'GGATTAAGATTGCAGTTCCTAG' #GGAGATTTCTACGAGCCAAA
gmk_F <- 'GATGGTCCCGATAAAACGAT'
pta_F <- 'CCTTCAGGTAATACGATTTTAAC'
tpi_F <- 'TTCACGACGTTCAGAATGAACGAG'
yqi_F <- 'GCCAATAGGTGTCCTGTATGCTG'

# file path -------------------------------------------------------------


# ----------------------------------------------------------
# Barcode typing ----------------------------------------------------------
# Function definiction

# Function: barcode_typing
barcode_tpying <- function(seq, barcodePosition){
  
  barcode1.1  <- Biostrings::DNAString(x="CATT")
  barcode1.2  <- Biostrings::DNAString(x="TGTT")
  barcode2.3  <- Biostrings::DNAString(x="GGTC")
  barcode2.4  <- Biostrings::DNAString(x="TCCC")
  barcode3.5  <- Biostrings::DNAString(x="GGCG")
  barcode3.6  <- Biostrings::DNAString(x="AGGT")
  barcode4.7  <- Biostrings::DNAString(x="GTCT")
  barcode4.8  <- Biostrings::DNAString(x="TCTG")
  barcode5.9  <- Biostrings::DNAString(x="TAAC")
  barcode5.10 <- Biostrings::DNAString(x="TAGT")
  barcode6.11 <- Biostrings::DNAString(x="TTAG")
  barcode6.12 <- Biostrings::DNAString(x="GCGG")
  
  if(barcodePosition == "barcode1"){
    barcode.index <- c("barcode1.1", "barcode1.2")
    pairing.index <- list(barcode1.1, barcode1.2) %>%
      purrr::map(~Biostrings::pairwiseAlignment(pattern = ., subject = Biostrings::DNAString(seq))) %>% 
      purrr::map(~score(.x)) %>% unlist() %>% base::which.max()
    return(barcode.index[pairing.index])
  }
  if(barcodePosition == "barcode2"){
    barcode.index <- c("barcode2.3", "barcode2.4")
    pairing.index <- list(barcode2.3, barcode2.4) %>%
      purrr::map(~Biostrings::pairwiseAlignment(pattern = ., subject = Biostrings::DNAString(seq))) %>% 
                   purrr::map(~score(.x)) %>% unlist() %>% base::which.max()
                 return(barcode.index[pairing.index])
  }
  if(barcodePosition == "barcode3"){
    barcode.index <- c("barcode3.5", "barcode3.6")
    pairing.index <- list(barcode3.5, barcode3.6) %>%
      purrr::map(~Biostrings::pairwiseAlignment(pattern = ., subject = Biostrings::DNAString(seq))) %>% 
                   purrr::map(~score(.x)) %>% unlist() %>% base::which.max()
                 return(barcode.index[pairing.index])
  }
  if(barcodePosition == "barcode4"){
    barcode.index <- c("barcode4.7", "barcode4.8")
    pairing.index <- list(barcode1.1, barcode1.2) %>%
      purrr::map(~Biostrings::pairwiseAlignment(pattern = ., subject = Biostrings::DNAString(seq))) %>% 
                   purrr::map(~score(.x)) %>% unlist() %>% base::which.max()
                 return(barcode.index[pairing.index])
  }
  if(barcodePosition == "barcode5"){
    barcode.index <- c("barcode5.9", "barcode5.10")
    pairing.index <- list(barcode5.9, barcode5.10) %>%
      purrr::map(~Biostrings::pairwiseAlignment(pattern = ., subject = Biostrings::DNAString(seq))) %>% 
                   purrr::map(~score(.x)) %>% unlist() %>% base::which.max()
                 return(barcode.index[pairing.index])
  }
  if(barcodePosition == "barcode6"){
    barcode.index <- c("barcode6.11", "barcode6.12")
    pairing.index <- list(barcode6.11, barcode6.12) %>%
      purrr::map(~Biostrings::pairwiseAlignment(pattern = ., subject = Biostrings::DNAString(seq))) %>% 
                   purrr::map(~score(.x)) %>% unlist() %>% base::which.max()
                 return(barcode.index[pairing.index])
  }
}

# Function: barcode_typing, Calling: barcode_typing
annotation_barcode_df <- function(barcode_df){
empty.vector <- c()
for( i in 1:length(barcode_df$seqID)){
  seqID        <- barcode_df$seqID
  sequencing.0 <- barcode_df$Amplicon
  barcode.0    <- barcode_df$BarcodeType
  print(i)
  barcode.result <- barcode_tpying(seq=sequencing.0[i], barcodePosition = barcode.0[i])
  empty.vector <- c(empty.vector, barcode.result)
}
  barcode_df %>% mutate(BarcodeResult=empty.vector) %>% return()
}

# Function: barcode_sampling, Calling:barcode_typing, annotation_barcode_df
barcode_sample <- function(barcode_df){
  
  print("# ======================== check 4 barcode reads")
  complete.catch.seqID <- barcode_df %>% filter(BarcodeType %in% c('barcode1', 'barcode2', 'barcode5','barcode6')) %>% 
    group_by(seqID, BarcodeType) %>% summarise(number=n()) %>% group_by(seqID) %>% summarise(number=n()) %>% filter(number==4) %>% with(seqID)
  
  filter.barcode_df <- barcode_df %>% filter(seqID %in% complete.catch.seqID)
  print("# ======================== coding only 4 barcode reads")
  filter.barcode_df.annotation <- annotation_barcode_df(barcode_df = filter.barcode_df )
  
  two.seqID <- filter.barcode_df.annotation %>%  dplyr::select(c(1,3,5)) %>% unique() %>% group_by(seqID, BarcodeType) %>% summarise(number=n()) %>% arrange(desc(number)) %>% filter(number > 1) %>% with(seqID)
  
  print("# ======================== filter the reads with 4 barcode result")
  filter.barcode_df.annotation.corrects <- filter.barcode_df.annotation  %>% 
    dplyr::select(c(1,3,5))%>% unique() %>% filter(!seqID %in% two.seqID) %>% 
    pivot_wider(names_from = BarcodeType, values_from =BarcodeResult)
  
  print("# ======================== label for each sample")
  FAR28891_pass_766e7ef0.1.correct.barcodes.persample <-  filter.barcode_df.annotation.corrects  %>% 
    dplyr::select(c(1,2,3,6,7)) %>%
    mutate(across(2:5, ~str_remove(.x, pattern="barcode\\d\\."))) %>% 
    mutate(sample_code=paste(barcode1, barcode2, barcode5, barcode6, sep='-')) %>% 
    mutate(sample=case_when(sample_code == "1-3-9-11" ~ "sample1",
                            sample_code == "2-3-9-11" ~ "sample2",
                            sample_code == "1-4-9-11" ~ "sample3",
                            sample_code == "2-4-9-11" ~ "sample4",
                            sample_code == "1-3-10-11" ~"sample5",
                            sample_code == "2-3-10-11" ~ "sample6",
                            sample_code == "1-4-10-11" ~ "sample7",
                            sample_code == "2-4-10-11" ~ "sample8",
                            sample_code == "1-3-9-12" ~ "sample9",
                            sample_code == "2-3-9-12" ~ "sample10"))
  return(FAR28891_pass_766e7ef0.1.correct.barcodes.persample)
}
# ----------------------------------------------------------
###############################################

FAR28891_pass_766e7ef0_1_barcode <- read_delim(file='/home/weitinglin66/Documents/analysis/202205_nanopore/20220514_amplicon/FAR28891_pass_766e7ef0_1_barcode.txt',
                                               col_names = FALSE) %>% 
  dplyr::rename("seqID"=X1, "Amplicon"=X2, "BarcodeType"=X3)


FAR28891_pass_766e7ef0_2_barcode <- read_delim(file='/home/weitinglin66/Documents/analysis/202205_nanopore/20220514_amplicon/FAR28891_pass_766e7ef0_2_barcode.txt',
                                               col_names = FALSE) %>% 
  dplyr::rename("seqID"=X1, "Amplicon"=X2, "BarcodeType"=X3)


# ----------------------------------------------------------
# Data input
###############################################
# First
FAR28891_pass.barcode.filespath.list <- list.files(path='/home/weitinglin66/Documents/analysis/202205_nanopore/20220514_amplicon/', full.names = TRUE) %>%  keep(~str_detect(.x, pattern="barcode"))

FAR28891_pass.barcode.filenames.list <- list.files(path='/home/weitinglin66/Documents/analysis/202205_nanopore/20220514_amplicon/', full.names = FALSE) %>%  keep(~str_detect(.x, pattern="barcode"))

FAR28891_pass.barcode.rawdf.list <- FAR28891_pass.barcode.filespath.list %>% map(~.x %>% read_delim(file=., col_names = FALSE) %>% dplyr::rename("seqID"=X1, "Amplicon"=X2, "BarcodeType"=X3))

FAR28891_pass.barcode.rawdf.list.annotation <- list(FAR28891_pass.barcode.rawdf.list,as.list(FAR28891_pass.barcode.filenames.list)) %>% pmap(~.x %>% dplyr::mutate(filename=.y))

FAR28891_pass_766e7ef0_0_barcode.persample <- barcode_sample(barcode_df=FAR28891_pass.barcode.rawdf.list.annotation[[1]])


# Second
FAR28891_pass.barcode.filespath.list <- list.files(path='/home/weitinglin66/Documents/analysis/202205_2_nanopore/202205_2_amplicon', full.names = TRUE) %>%  keep(~str_detect(.x, pattern="barcode"))

FAR28891_pass.barcode.filenames.list <- list.files(path='/home/weitinglin66/Documents/analysis/202205_2_nanopore/202205_2_amplicon', full.names = FALSE) %>%  keep(~str_detect(.x, pattern="barcode"))

FAR28891_pass.barcode.rawdf.list <- FAR28891_pass.barcode.filespath.list %>% map(~.x %>% read_delim(file=., col_names = FALSE) %>% dplyr::rename("seqID"=X1, "Amplicon"=X2, "BarcodeType"=X3))

FAR28891_pass.barcode.rawdf.list.annotation <- list(FAR28891_pass.barcode.rawdf.list,as.list(FAR28891_pass.barcode.filenames.list)) %>% pmap(~.x %>% dplyr::mutate(filename=.y))



tmp.list <- list()
for( i in 1:53){
  print(i)
  tmp.filename <- unique(FAR28891_pass.barcode.rawdf.list.annotation[[i]]$filename)
  tmp.barcode.persample <- barcode_sample(barcode_df=FAR28891_pass.barcode.rawdf.list.annotation[[i]]) %>% mutate(filename=tmp.filename)
  tmp.list <- tmp.list %>%  purrr::prepend(., list(tmp.barcode.persample), before = NULL)
}


FAR28891.file.names <- FAR28891_pass.barcode.filenames.list %>% 
                              str_remove(., pattern = "FAR28891_pass_766e7ef0_") %>% 
                              str_remove(., pattern = "_barcode.txt")
unique(tmp.list[[1]]$filename) %>% 
  str_remove(., pattern = "FAR28891_pass_766e7ef0_") %>% 
  str_remove(., pattern = "_barcode.txt")

sample.list <- c('sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6', 'sample7', 'sample8', 'sample9', 'sample10')
# ----------------------------------------------------------
# Output sample ID files -----------------------------------
# ===========================================================
  
  for (i in 1:10){
    print(paste("Processing", j, " files"))
    file.order.name <-  unique(tmp.list[[j]]$filename) %>% 
                              str_remove(., pattern = "FAR28891_pass_766e7ef0_") %>% 
                              str_remove(., pattern = "_barcode.txt")
    file.names <- paste0('/home/weitinglin66/Documents/analysis/202205_2_nanopore/202205_2_sampleID/20220515_2_',file.order.name,'_', sample.list[i], '_SeqID.txt')
    tmp.list[[j]] %>% filter(sample == sample.list[i]) %>%
      dplyr::select(1) %>% write_delim(file = file.names, col_names = FALSE)
  }


# -----------------------------------------------------------