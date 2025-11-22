######################################
#
#  | 202204514
#  | Weitinglin66
#  | Oxford nanopore demultiplex/barcode split
#
######################################



# Loading library ---------------------------------------------------------
library(tidyverse)
library(Biostrings)
library(ggmsa)
library(msa)


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




# Data Input --------------------------------------------------------------
FAR28891_pass_766e7ef0_0_barcode <- read_delim(file='/home/weitinglin66/Documents/analysis/202205_nanopore/20220514_amplicon/FAR28891_pass_766e7ef0_0_barcode.txt',
                                               col_names = FALSE) %>% 
                                      dplyr::rename("seqID"=X1, "Amplicon"=X2, "BarcodeType"=X3)
            
FAR28891_pass_766e7ef0_1_barcode <- read_delim(file='/home/weitinglin66/Documents/analysis/202205_nanopore/20220514_amplicon/FAR28891_pass_766e7ef0_1_barcode.txt',
                                               col_names = FALSE) %>% 
                                      dplyr::rename("seqID"=X1, "Amplicon"=X2, "BarcodeType"=X3)

FAR28891_pass_766e7ef0_0_count <- read_delim(file='/home/weitinglin66/Documents/analysis/202205_nanopore/20220510_count/20220510_first_FAR28891/FAR28891_pass_766e7ef0_0_length.txt',
                                             col_names = FALSE) %>% filter(X2>3000 & X2 <4000) %>% 
                                  rename("seqID"=X1, "length"=X2)


FAR28891_pass_766e7ef0_0_locate.raw <- read_delim(file='/home/weitinglin66/Documents/analysis/202205_nanopore/20220514_locate/FAR28891_pass_766e7ef0_0_locate.txt') %>%
                                       mutate(primer=case_when(patternName == 'CGCTGATTGAAGCAGATACCT'     ~ "arc_R",
                                                                patternName == 'GATATCGTTATTAATACAACACC'  ~ "aro_R",   
                                                                patternName == 'GAATTGGACATGCGATTTTACCA'  ~ "glp_R",
                                                                patternName == 'TACGATTACGTTGTAGTTAATGA'  ~ "gmk_R", 
                                                                patternName == 'TTAAGCTTTTCAACAAAAGGGTC'  ~ "pta_R",
                                                                patternName == 'GTACAATTGTTAGAAGGTGCAAA'  ~ "tpi_R", 
                                                                patternName == 'GAATGTGAAATAGGATTTCCGAT'  ~ "aro_F",
                                                                patternName == 'GGATTAAGATTGCAGTTCCTAG'   ~ "glp_F", 
                                                                patternName == 'GATGGTCCCGATAAAACGAT'     ~ "gmk_F",
                                                                patternName == 'CCTTCAGGTAATACGATTTTAAC'  ~ "pta_F",  
                                                                patternName == 'TTCACGACGTTCAGAATGAACGAG' ~ "tpi_F",
                                                                patternName == 'GCCAATAGGTGTCCTGTATGCTG'  ~ "yqi_F")) 

FAR28891_pass_766e7ef0_0_locate.raw %>% group_by(seqID) %>% summarise(number=n())

FAR28891_pass_766e7ef0_0_locate.raw


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

###############################################

FAR28891_pass_766e7ef0_1_barcode <- read_delim(file='/home/weitinglin66/Documents/analysis/202205_nanopore/20220514_amplicon/FAR28891_pass_766e7ef0_1_barcode.txt',
                                               col_names = FALSE) %>% 
  dplyr::rename("seqID"=X1, "Amplicon"=X2, "BarcodeType"=X3)

FAR28891_pass_766e7ef0_2_barcode <- read_delim(file='/home/weitinglin66/Documents/analysis/202205_nanopore/20220514_amplicon/FAR28891_pass_766e7ef0_2_barcode.txt',
                                               col_names = FALSE) %>% 
  dplyr::rename("seqID"=X1, "Amplicon"=X2, "BarcodeType"=X3)

FAR28891_pass_766e7ef0_2_barcode <- read_delim(file='/home/weitinglin66/Documents/analysis/202205_nanopore/20220514_amplicon/FAR28891_pass_766e7ef0_2_barcode.txt',
                                               col_names = FALSE) %>% 
  dplyr::rename("seqID"=X1, "Amplicon"=X2, "BarcodeType"=X3)

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

for (j in 1:53){
  
  for (i in 1:10){
    print(paste("Processing", j, " files"))
    file.order.name <-  unique(tmp.list[[j]]$filename) %>% 
                              str_remove(., pattern = "FAR28891_pass_766e7ef0_") %>% 
                              str_remove(., pattern = "_barcode.txt")
    file.names <- paste0('/home/weitinglin66/Documents/analysis/202205_2_nanopore/202205_2_sampleID/20220515_2_',file.order.name,'_', sample.list[i], '_SeqID.txt')
    tmp.list[[j]] %>% filter(sample == sample.list[i]) %>%
      dplyr::select(1) %>% write_delim(file = file.names, col_names = FALSE)
  }
}

##############################################

identified.reads <- tmp.list %>% map(nrow) %>% unlist
id               <- tmp.list %>% map('filename') %>% map(unique) %>% map(str_remove, pattern="FAR28891_pass_766e7ef0_") %>% map(str_remove, pattern="_barcode.txt") %>% unlist %>% as.numeric()

tibble(time=id, identified_reads=identified.reads) %>% ggplot(data=.) + geom_col(aes(x=time, y=identified_reads))
##############################################

sample1.file20 <- readDNAStringSet(filepath="/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_assemlby/20220515_20_sample_1/assembly_20_sample1.fasta", format="fasta")
sample1.file21 <- readDNAStringSet(filepath="/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_assemlby/20220515_21_sample_1/assembly_21_sample1.fasta", format="fasta")
sample1.file22 <- readDNAStringSet(filepath="/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_assemlby/20220515_22_sample_1/assembly_22_sample1.fasta", format="fasta")
sample1.file23 <- readDNAStringSet(filepath="/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_assemlby/20220515_23_sample_1/assembly_23_sample1.fasta", format="fasta")
sample1.file24 <- readDNAStringSet(filepath="/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_assemlby/20220515_24_sample_1/assembly_24_sample1.fasta", format="fasta")



sample1.set <- DNAStringSet(c(sample1.file21,
                              sample1.file24))
names(sample1.set) <- c("file21", "file24")

sample1.set.msa <- msa(sample1.set,
                       method="ClustalW",
                       type="dna",
                       verbose = TRUE)

sample1.set.msa <- msa(sample1.set,
                       method="Muscle",
                       type="dna",
                       verbose = TRUE)


sample1.set.msa@unmasked %>% ggmsa(msa=., start = 3000, end=3200)

pairwiseAlignment(sample1.file21, sample1.file24) %>% consensusString()
pairwiseAlignment(pattern=arc_F, subject=sample1.file21, type="local") 
pairwiseAlignment(pattern=arc_F, subject=sample1.file21, type="overlap") 
pairwiseAlignment(pattern=arc_F, subject=sample1.file21, type="overlap") 
pairwiseAlignment(pattern=arc_R.revc, subject=sample1.file21, type="overlap") 
pairwiseAlignment(pattern=arc_F, subject=sample1.file21, type="local-global")

arc_R.revc<- arc_R %>% DNAString() %>% reverseComplement()
###################################################################################
sample1.file20 <- readDNAStringSet(filepath="/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_assemlby/20220515_20_sample_1/assembly_20_sample1.fasta", format="fasta")
sample1.file21 <- readDNAStringSet(filepath="/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_assemlby/20220515_21_sample_1/assembly_21_sample1.fasta", format="fasta")
sample1.file22 <- readDNAStringSet(filepath="/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_assemlby/20220515_22_sample_1/assembly_22_sample1.fasta", format="fasta")

arc_F <- 'TTGATTCACCAGCGCGTATTGTC' #arc up - 5' TTG ATT CAC CAG CGC GTA TTG TC -3'

arc_R <- 'AGGTATCTGCTTCAATCAGCG'

yqi_F <- 'GCCAATAGGTGTCCTGTATGCTG'
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

pairwiseAlignment(glp_F, sample1.file22)
pairwiseAlignment(glp_R, sample1.file22)


sample1.genes.filepaths <- list.files(path = "/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_sample_MLST/sample1/", full.names = TRUE)
sample1.genes.dans <- sample1.genes.filepaths %>% map(., readDNAStringSet, format="fasta")
sample1.genes.msa <- sample1.genes.dans %>% map(., 
                                                msa,
                                                method="ClustalW",
                                                type="dna",
                                                verbose = TRUE)

sample1.genes.msa[[2]]
