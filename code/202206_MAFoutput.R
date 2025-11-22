######################################
#
#  | 20220611
#  | Weitinglin66
#  | Oxford nanopore demultiplex/barcode split
#
######################################



# Loading library ---------------------------------------------------------
library(tidyverse)
library(Biostrings)
library(ggmsa)
library(msa)









sample1.file20 <- readDNAStringSet(filepath="/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_assemlby/20220515_20_sample_1/assembly_20_sample1.fasta", format="fasta")
names(sample1.file20) <- "sample1_file20"
sample1.file21 <- readDNAStringSet(filepath="/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_assemlby/20220515_21_sample_1/assembly_21_sample1.fasta", format="fasta")
names(sample1.file21) <- "sample1_file21"
sample1.file22 <- readDNAStringSet(filepath="/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_assemlby/20220515_22_sample_1/assembly_22_sample1.fasta", format="fasta")
names(sample1.file22) <- "sample1_file22"
sample1.file23 <- readDNAStringSet(filepath="/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_assemlby/20220515_23_sample_1/assembly_23_sample1.fasta", format="fasta")
names(sample1.file23) <- "sample1_file23"
sample1.file24 <- readDNAStringSet(filepath="/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_assemlby/20220515_24_sample_1/assembly_24_sample1.fasta", format="fasta")
names(sample1.file24) <- "sample1_file24"

sample1.consensus.file20 <- readDNAStringSet(filepath="/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_polish/consensus.fasta", format="fasta")
names(sample1.consensus.file20) <- "sample1_consensus"
sample1.consensus.file21 <- readDNAStringSet(filepath="/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_polish/sample1_21_consensus.fasta", format="fasta")
names(sample1.consensus.file21) <- "file21_sample1_consensus"
sample1.consensus.file22 <- readDNAStringSet(filepath="/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_polish/sample1_22_consensus.fasta", format="fasta")
names(sample1.consensus.file22) <- "file22_sample1_consensus"
sample1.consensus.file23 <- readDNAStringSet(filepath="/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_polish/sample1_23_consensus.fasta", format="fasta")
names(sample1.consensus.file23) <- "file23_sample1_consensus"
sample1.consensus.file24 <- readDNAStringSet(filepath="/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_polish/sample1_24_consensus.fasta", format="fasta")
names(sample1.consensus.file24) <- "file24_sample1_consensus"

pairwiseAlignment(sample1.file20, sample1.file21, sample1.consensus.file20)

sample1.genes.msa <- c(sample1.file20, sample1.file21, sample1.file22, sample1.file23, sample1.file24, sample1.consensus.file20,
                       sample1.consensus.file21, sample1.consensus.file22, sample1.consensus.file23, sample1.consensus.file24)%>% msa(.,
                                                method="ClustalW",
                                                type="dna",
                                                verbose = TRUE)

colM <- IRanges(start=1, end=100)
colmask(sample1.genes.msa) <- colM
sample1.genes.msa
sample1.genes.msa<- unmasked(sample1.genes.msa)

sink("/media/weitinglin66/new202205/analysis/202205_2_nanopore/202205_2_polish/file_20_sample1.txt")
print(sample1.genes.msa, showNames=TRUE, showConsensus=TRUE, show="complete")
sink()

print(sample1.genes.msa, showNames=TRUE, showConsensus=TRUE, show="complete")



# calculation the experiment result ---------------------------------------

experment.files.num.path.list <- list(list(10,15,20,25,30), list(1:100)) %>% pmap(~paste0("/media/weitinglin66/new202205/analysis/202205_2_nanopore/202206_2_experiment/sample4_assembly",
                                              "/sample4_",.x,"_assembly/sample4_",.x,"_assembly_",.y,"/sample4_",.x,"_assembly_",.y,"_polish_consensus.fasta")) 


experment.files.reads10 <- experment.files.num.path.list[1] %>% map(., readDNAStringSet, format="fasta", use.names=TRUE)

names(experment.files.reads10) <- list(1:100) %>% map(.,~paste0("sample4_reads10_num",.)) 

experment.files.reads10 %>%  map(., 
    msa,
    method="ClustalW",
    type="dna",
    verbose = TRUE)




# Calculation the figure 2 ------------------------------------------------


sample4_10to100_experiemnt_genotyping.raw <- read_delim(file='/media/weitinglin66/new202205/analysis/202205_2_nanopore/202206_2_result/sample4_genotype_result.txt', 
                                                        delim = '\t',
                                                        col_names = FALSE) %>% rename('X1'='Gene_type',
                                                                                      'X2'='SampleSource',
                                                                                      'X3'='Cigar',
                                                                                      'X4'='Tag_NM',
                                                                                      'X5'='Tag_MS',
                                                                                      'X6'='Tag_AS',
                                                                                      'X7'='Tag_NN',
                                                                                      'X8'='Tag_TP') %>% 
                                              tidyr::separate(., col=Gene_type, into=c("Gene", "Type"), sep="_") %>% 
                                              tidyr::separate(., col=SampleSource, into=c("Sample", "Reads", "Num"), sep="_") %>% 
                                              mutate(Cigar_check=!str_detect(Cigar, pattern = 'S')) 


sample4_10to100_experiemnt_genotyping.raw %>% 
  filter(Cigar_check==TRUE) %>%
  group_by(Reads, Gene, Type) %>%
  summarise(number=n()) %>% mutate(positive_rate=number/100) %>% 
  ggplot(data=.) + geom_col(aes(x=Reads, y=positive_rate, fill=Gene)) + facet_grid(rows = vars(Gene)) +
  scale_y_continuous(labels = scales::percent_format(accuracy=1)) + labs(x='Total Reads used', y="Typing Correct Rate")

sample4_10to100_experiemnt_genotyping.raw %>% 
  filter(Cigar_check==TRUE) %>%
  group_by(Reads, Gene, Type) %>%
  summarise(number=n(), .groups = 'drop') %>% mutate(positive_rate=number/100, Reads = as.integer(Reads))  %>% 
  ggplot(data=.) + geom_line(aes(x=Reads, y=positive_rate, color=Gene)) + facet_grid(rows = vars(Gene)) +
  scale_y_continuous(labels = scales::percent_format(accuracy=1)) +
  scale_x_continuous(breaks = c(5+5*(1:19))) +
  labs(x='Total Reads used', y="Typing Correct Rate")

sample4_10to100_experiemnt_genotyping.raw %>% filter(Reads==25 && Gene == "aroE")



sample4_25_experiemnt_genotyping.raw <- read_delim(file='/media/weitinglin66/new202205/analysis/202205_2_nanopore/202206_2_result/sample4_25_genotype_result.txt', 
                                                        delim = '\t',
                                                        col_names = FALSE) %>% rename('X1'='Gene_type',
                                                                                      'X2'='SampleSource',
                                                                                      'X3'='Cigar',
                                                                                      'X4'='Tag_NM',
                                                                                      'X5'='Tag_MS',
                                                                                      'X6'='Tag_AS',
                                                                                      'X7'='Tag_NN',
                                                                                      'X8'='Tag_TP') %>% 
  tidyr::separate(., col=Gene_type, into=c("Gene", "Type"), sep="_") %>% 
  tidyr::separate(., col=SampleSource, into=c("Sample", "Reads", "Num"), sep="_") %>% 
  mutate(Cigar_check=!str_detect(Cigar, pattern = 'S'))


sample4_patch_experiemnt_genotyping.raw <- read_delim(file='/media/weitinglin66/new202205/analysis/202205_2_nanopore/202206_2_result/sample4_patch_genotype_result.txt', 
                                                        delim = '\t',
                                                        col_names = FALSE) %>% rename('X1'='Gene_type',
                                                                                      'X2'='SampleSource',
                                                                                      'X3'='Cigar',
                                                                                      'X4'='Tag_NM',
                                                                                      'X5'='Tag_MS',
                                                                                      'X6'='Tag_AS',
                                                                                      'X7'='Tag_NN',
                                                                                      'X8'='Tag_TP') %>% 
                                                      tidyr::separate(., col=Gene_type, into=c("Gene", "Type"), sep="_") %>% 
                                                      tidyr::separate(., col=SampleSource, into=c("Sample", "Reads", "Num"), sep="_") %>% 
                                                      mutate(Cigar_check=!str_detect(Cigar, pattern = 'S')) 

# ===================== focus on the glpF
sample4_patch_experiemnt_genotyping.raw %>% filter(Gene == 'glpF' & Type == '1') %>% 
  filter(Cigar_check==TRUE) %>%
  group_by(Reads, Gene, Type) %>%
  summarise(number=n()) %>% mutate(positive_rate=number/100) %>% 
  ggplot(data=.) + geom_col(aes(x=Reads, y=positive_rate, fill=Gene)) + facet_grid(rows = vars(Gene)) +
  scale_y_continuous(labels = scales::percent_format(accuracy=1)) + labs(x='Total Reads used', y="Typing Correct Rate")

# ===================== on arcC, aroR, glpF, gmk, tpi, yqiL, tpi
sample4_patch_experiemnt_genotyping.raw %>% 
  filter(Cigar_check==TRUE) %>%
  group_by(Reads, Gene, Type) %>%
  summarise(number=n()) %>% mutate(GeneType=paste(Gene, Type, sep="_")) %>%
  filter(GeneType %in% c("arcC_2", "aroE_3", "glpF_1", "gmk_1", "tpi_19", "yqiL_3","tpi_4")) %>% 
  mutate(positive_rate=number/100) %>% 
  ggplot(data=.) + geom_col(aes(x=Reads, y=positive_rate, fill=Gene)) + facet_grid(rows = vars(Gene)) +
  scale_y_continuous(labels = scales::percent_format(accuracy=1)) + labs(x='Total Reads used', y="Typing Correct Rate")

sample4_patch_experiemnt_genotyping.raw %>% 
  filter(Cigar_check==TRUE) %>%
  group_by(Reads, Gene, Type) %>%
  summarise(number=n(), .groups = 'drop') %>% mutate(GeneType=paste(Gene, Type, sep="_")) %>%
  filter(GeneType %in% c("arcC_2", "aroE_3", "glpF_1", "gmk_1", "tpi_4", "yqiL_3", "pta_4")) %>% 
  mutate(positive_rate=number/100, Reads = as.integer(Reads))  %>% 
  ggplot(data=.) + 
  geom_line(aes(x=Reads, y=positive_rate, color=Gene)) + 
  geom_point(aes(x=Reads, y=positive_rate, color=Gene)) + 
  geom_text(aes(x=Reads, y=positive_rate, label=scales::percent(positive_rate)), size = 5, vjust = 1.5, check_overlap = TRUE) +
  facet_grid(rows = vars(Gene)) +
  scale_y_continuous(labels = scales::percent_format(accuracy=1), limits = c(0,1)) +
  scale_x_continuous(breaks = c(10:60), limits = c(10,60)) +
  labs(x='Total Reads used', y="Typing Correct Rate")






sample4_20.fastq.check <- read_delim(file='/media/weitinglin66/new202205/analysis/202205_2_nanopore/202206_2_experiment/sample4_check/sample4_11_fastqcheck.txt', 
                                                   delim = '\t',
                                                   col_names = FALSE) %>% rename('X1'='SeqID',
                                                                                 'X2'='SampleSource',
                                                                                 'X3'='Seqlength')



sample4_20.fastq.check %>% group_by(SampleSource) %>% summarise(mean_seq=mean(Seqlength),
                                                                median_seq=median(Seqlength))
sample4_20.fastq.check %>% group_by(SampleSource) %>% View()

set.1 <- sample4_20.fastq.check %>% filter(SampleSource == 'sample4_11_1') %>% with(SeqID)
set.2 <- sample4_20.fastq.check %>% filter(SampleSource == 'sample4_11_2') %>% with(SeqID)
setdiff(set.1, set.2)
intersect(set.1, set.2)
