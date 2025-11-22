####################################|
###
### Author: weitinglin66
### Date:2021/11/13
### Topics: Golden Gate-aided strategy
### Log:
###  - for case data:22331936
###. - update on 20251102
####################################|
# some tools related to primer design: rmelting https://github.com/aravind-j/rmelting
# BiocManager::install("openPrimeR")
# BiocManager::install("rmelting")
# NEB melting API: 


# 01:load library ------------------------------------------------------------

#install.packages("BiocManager")

library(tidyverse)
library(sangerseqR) # BiocManager::install("sangerseqR")
library(seqinr) #BiocManager::install("seqinr")
library(Biostrings) # BiocManager::install("Biostrings")
library(ggmsa) # BiocManager::install("ggmsa")
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)

# input data from publication's files
# 2020. Enabling one-pot Golden Gate assemblies of unprecedented complexity using data-optimized assembly design
#

bsaI.HF.overhange.data <- read_excel(path='/Users/weitingm4/Library/Mobile Documents/com~apple~CloudDocs/2025 Yeh/LigationRestriction_data/pone.0238592.s001.xlsx')

bsaI.HF.overhange.df <- bsaI.HF.overhange.data %>% 
                            pivot_longer(data=., cols=c(2:257), names_to="PairedSeq") %>% 
                            mutate(ComplementSeq=as.character(reverse(complement(DNAStringSet(Overhang))))) %>% 
                            mutate(Correct=(PairedSeq==ComplementSeq)) 

bsaI.HF.overhange.fidelityRate <-   bsaI.HF.overhange.df %>% 
                                        group_by(Overhang, Correct) %>% summarise(value=sum(value)) %>% 
                                        pivot_wider(names_from = "Correct", values_from="value") %>% 
                                        mutate(total=`FALSE`+`TRUE`) %>% 
                                        mutate(fidelityRate=`TRUE`/total) %>% 
                                        arrange(desc(fidelityRate))

rank <- 1:256
bsaI.HF.overhange.fidelityRate %>%
  ungroup %>% 
  mutate(rank = 257 - min_rank(fidelityRate)) 

bsaI.HF.overhange.fidelityRate %>% filter(Overhang == 'TTAG')

# Visualization of data for publication
# --- Assuming your data is loaded into a data frame named 'df' ---

# 1. Tidy the data into long format
long_data <- bsaI.HF.overhange.data %>%
  # Select all columns *except* the first one (Overhang) for pivoting
  pivot_longer(
    cols = -Overhang,
    names_to = "Following_Sequence",
    values_to = "Count"
  ) %>%
  # Ensure Count is numeric for proper visualization
  mutate(Count = as.numeric(Count))


sample.overhang<- long_data %>% with(Overhang) %>% unique() %>% sample(., size=248)

best.barcode.space <- c("AATG","AACA","GACC","GGGA","CGCC","ACCT","AGAC","CAGA","GTTA","ACTA","CTAA","CCGC")
current.barcode.MLST <- c("CATT", "TGTT", "GGTC", "GGTC", "TCCC", "GGCG", "AGGT", "GTCT", "TCTG", "TAAC", "TAGT", "TTAG", "GCGG")
good.barcode.candidate <- c(best.barcode.space,current.barcode.MLST)

filter.long_data <- long_data %>% 
                         filter(Overhang %in% good.barcode.candidate) %>% 
                         filter(Following_Sequence %in% good.barcode.candidate)

# Reordering the figure
# Assuming 'df' is your wide-format data (Overhang, AAAA, AAAG, etc.)
# 1. Prepare matrix for clustering
mat <- bsaI.HF.overhange.data %>% 
  tibble::column_to_rownames(var = "Overhang") %>% 
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

# 2. Hierarchical Clustering
row_order <- hclust(dist(mat))$order    # Cluster rows (Overhangs)
col_order <- hclust(dist(t(mat)))$order # Cluster columns (Following Sequences, need to transpose)

# 3. Extract the sequence names in the clustered order
ordered_rows <- rownames(mat)[row_order]
ordered_cols <- colnames(mat)[col_order]

# Reorder the factors in your long_data for ggplot2
long_data_ordered <- filter.long_data %>%
  mutate(
    Overhang = factor(Overhang, levels = ordered_rows),
    Following_Sequence = factor(Following_Sequence, levels = ordered_cols)
  )



# 2. Create the Heatmap
my_heatmap <- ggplot(long_data_ordered, aes(x = Following_Sequence, y = Overhang, fill = Count)) +
  # Create the tiles (the heatmap cells)
              geom_tile(color = "white", linewidth = 0.5) +
  # Optionally, add the count value as text inside each tile
  #            geom_text(aes(label = Count), color = "black", size=1) +
  # Use a perceptually uniform color scale (e.g., viridis)
              scale_fill_gradient(low = "white", high = "darkred") +
  # Set the theme and labels
              labs(
                title = "Barcode Ligase Interaction Map",
                x = "Paired Sequence",
                y = "Overhang Sequence"
                  ) +
  # Improve aesthetics
                theme_minimal() +
  # Rotate X-axis labels for better readability
                  theme(
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),  # size 3 good for larger data set
                      axis.text.y = element_text(size = 10),  # size 3 good for larger data set
                      panel.grid.major = element_blank(), # Remove grid lines
                      panel.grid.minor = element_blank()
                        )

# 2. Save the plot using ggsave()
ggsave(
  filename = "/Users/weitingm4/Library/Mobile Documents/com~apple~CloudDocs/2025 Yeh/Code/sequence_overlap_heatmap.png", # The desired filename and format
  plot = my_heatmap,                       # The plot object you want to save
  device = "png",                           # The file format (optional, often inferred from filename)
  width = 10,                               # Width of the output image in inches
  height = 8,                               # Height of the output image in inches
  units = "in",                             # Units for width/height (e.g., "in", "cm", "mm")
  dpi = 300                                 # Resolution (dots per inch) for high quality
)

# 02 ----------------------------------------------------------------------

BsaI.HF.v2.recognization.seq <- DNAString("GGTCTC")


candidate.string.list <- c()
num <- 2

for (i in 1:num){
  tmp.overhange.fidelityRate <- bsaI.HF.overhange.fidelityRate %>% mutate(rank=1:256)
  
}

bsaI.HF.overhange.df




# design primer -----------------------------------------------------------
BsaI.HF.recognized.sequence <- "GGTCTC"
arcC_R <- "AGGTATCTGCTTCAATCAGCG"
aroE_F <- "ATCGGAAATCCTATTTCACATTC"
aroE_R <- "GGTGTTGTATTAATAACGATATC"
glpF_F <- "CTAGGAACTGCAATCTTAATCC"
pta_R  <- "GACCCTTTTGTTGAAAAGCTTAA"
tpi_F  <- "TCGTTCATTCTGAACGTCGTGAA"
# use tool from NEB design golden gate
best.barcode.space <- c("AATG","AACA","GACC","GGGA","CGCC","ACCT","AGAC","CAGA","GTTA","ACTA","CTAA","CCGC")
current.barcode.MLST <- c("CATT", "TGTT", "GGTC", "GGTC", "TCCC", "GGCG", "AGGT", "GTCT", "TCTG", "TAAC", "TAGT", "TTAG", "GCGG")

nanoMLST.arc_R_seq <- "AACCTTGCAAAGTGTCAGCA"

nanoMLST.arc_R_seq %>% Biostrings::DNAString() %>% complement()




bsaI.HF.overhange.fidelityRate %>% 
  filter(!Overhang %in% current.barcode.MLST) %>%
  filter(!Overhang %in% best.barcode.space) %>% arrange(desc(fidelityRate),desc(total)) %>% View()

chose.three.barcode <- c("TAGA", "TCTA", "TATA")



# Calculation of the low interfereance ------------------------------------

Overhang.seq.vector <- c(best.barcode.space, current.barcode.MLST)
Choose.barcode.pool <- bsaI.HF.overhange.fidelityRate %>% filter(!Overhang %in% Overhang.seq.vector) %>% with(Overhang)
all.combination.from.pool 
tmp.combination <- tibble(seq1=Choose.barcode.pool) %>% 
                    tidyr::expand(seq1, seq2=seq1, seq3=seq2) %>% 
                    filter(seq1 < seq2 & seq2 < seq3)

tmp.index.combination<- tmp.combination %>% mutate(index=1:2054360)
overall.calculation.result.df <- tibble()
sample.number <- sample(c(101:2054360), 1000, replace = FALSE)
sample.number <- c(101)
for(item in sample.number){
print(item)
      tmp.three.seq <- tmp.index.combination %>% 
                            filter(index==item) %>%
                            pivot_longer(col=c(1:3),names_to = "num", values_to ="seq") %>% with(seq)
      result.df <- tibble()
      tmp.Overhang.seq.vector <- c(Overhang.seq.vector, tmp.three.seq)
      for (i in 1:length(tmp.Overhang.seq.vector)){
      
      Overhang.seq <- tmp.Overhang.seq.vector[i]
      tmp.df <- bsaI.HF.overhange.df %>% 
                    filter(Overhang==Overhang.seq & PairedSeq %in% c(best.barcode.space, current.barcode.MLST)) %>% 
                    group_by(Correct) %>% summarise(Num=sum(value)) %>% mutate(Overhang=Overhang.seq)
      result.df <-result.df %>%  bind_rows(tmp.df)
        
      }
      
      cal.result.df <- result.df %>% group_by(Correct) %>%
                          summarise(num=sum(Num)) %>% 
                          pivot_wider(names_from = Correct, values_from = num) %>% 
                          mutate(index=item)
      overall.calculation.result.df <- overall.calculation.result.df %>% bind_rows(cal.result.df)
      
}

overall.calculation.result.df %>% write_csv(., file='Documents/2020_LinkouCGMH_CP/2021_研究/202112_GoldenGateMLST/20220606_barcode.csv')
input.20220606.overall.calculation.result.df <- read_csv(file="Documents/2020_LinkouCGMH_CP/2021_研究/202112_GoldenGateMLST/20220606_barcode.csv")

#=====================
BsaI.HF.recognized.sequence <- "GGTCTC"
arcC_R <- "AGGTATCTGCTTCAATCAGCG"
aroE_F <- "ATCGGAAATCCTATTTCACATTC"
aroE_R <- "GGTGTTGTATTAATAACGATATC"
glpF_F <- "CTAGGAACTGCAATCTTAATCC"
pta_R  <- "GACCCTTTTGTTGAAAAGCTTAA"
tpi_F  <- "TCGTTCATTCTGAACGTCGTGAA"

best.bet.indexs <- overall.calculation.result.df %>% arrange(`FALSE`) %>% filter(`FALSE` == 0) %>% with(index)

barcode.set <- tmp.index.combination %>% filter(index == best.bet.indexs[4])
# ATGA, ATGC, CCTA
# set13
paste0(BsaI.HF.recognized.sequence, "A","TCAT", arcC_R)
DNAString("ATGA") %>% reverseComplement()
paste0(BsaI.HF.recognized.sequence, "T", "ATGA", aroE_F)

# set14 
paste0(BsaI.HF.recognized.sequence, "A","GCAT", aroE_R)
DNAString("ATGC") %>% reverseComplement()
paste0(BsaI.HF.recognized.sequence, "T","ATGC", glpF_F)


# set15
paste0(BsaI.HF.recognized.sequence, "A", "TAGG", pta_R)
DNAString("CCTA") %>% reverseComplement()
paste0(BsaI.HF.recognized.sequence, "T", "CCTA", tpi_F)



# design S. lugdunensis ---------------------------------------------------

# design primer -----------------------------------------------------------
BsaI.HF.recognized.sequence <- "GGTCTC"
aroE_F <- "ATCGGAGATCCGATTTCACATTC"
aroE_R <- "GGCGTTGTATTAATTATAATATC"
dat_F <- "TCGTGGTTATGTTTTTGGTGACGGT"
dat_R <- "CTATGAGAAGTAAAGCCAGGAAT"
ddl_F  <- "AGTGCGGAGCACGACGTTTCA"
ddl_R  <- "ACACTTGATCCCAAATTCGCCGGT"
gmk_F  <- "ATAGTTCTTTCCGGACCATC"
gmk_R  <- "TCATTGACTACAACGTAATCATA"
idh_F  <- "ACTTGCAGGTGCCACGTCGA"
idh_R  <- "GCTACGCATTTGCAATGGTAACGCA"
recA_F  <- "GCACGGCCACCAGGTGTTGT"
recA_R  <- "AGGCCGTCGCGTATCTAGTGT"
yqiL_F  <- "GTGCTAAACGCACACCAATTGGA"
yqiL_R  <- "CTTCAGCATCGATTAGAGGCAC"
# use tool from NEB design golden gate
best.barcode.space <- c("AATG","AACA","GACC","GGGA","CGCC","ACCT","AGAC","CAGA","GTTA","ACTA","CTAA","CCGC")
current.barcode.MLST <- c("CATT", "TGTT", "GGTC", "GGTC", "TCCC", "GGCG", "AGGT", "GTCT", "TCTG", "TAAC", "TAGT", "TTAG", "GCGG")

coding.MLST.sl.table <- tibble(name=c("aroE_R","dat_F", 
                                      "dat_R","ddl_F", 
                                      "ddl_R","gmk_F",
                                      "gmk_R", "idh_F",
                                      "idh_R","recA_F",
                                      "recA_R", "yqiL_F",
                                      "aroE_R","dat_F", 
                                      "dat_R","ddl_F", 
                                      "ddl_R","gmk_F",
                                      "gmk_R", "idh_F",
                                      "idh_R","recA_F",
                                      "recA_R", "yqiL_F"),
                               code=c(1,1,3,3,5,5,7,7,9,9,11,11,
                                      2,2,4,4,6,6,8,8,10,10,12,12))

bsaI.HF.overhange.fidelityRate %>% arrange(desc(fidelityRate), desc(total)) %>% View()
# 重新挑選
barcode.code.sl.table <- tibble(code=c(1,2,3,4,5,6,7,8,9,10,11,12),
                                F=c("CATT",# 1
                                      "TGTT",# 2
                                      "GGTC",# 3
                                      "TCCC",# 4
                                      "GGCG",# 5  =/ 8
                                      "AGGT",# 6  =/ 7
                                      "ACTA",# 7 GTCT => ACTA
                                      "AAGA",# 8 TCTG => AAGA
                                      "TAAC",# 9
                                      "TAGT",# 10
                                      "TTAG",# 11
                                      "GCGG")) # 12

barcode.code.sl.revComplement.vector <- barcode.code.sl.table %>%  as.list() %>% with(F) %>% as.list() %>% map(., DNAString) %>% map(., reverseComplement) %>% map(as.character) %>% unlist

barcode.seq.annotation.table <- barcode.code.sl.table %>% mutate(R=barcode.code.sl.revComplement.vector) %>% pivot_longer(cols = c(2,3),
                                                                                                names_to="Type",
                                                                                                values_to="Seq") %>% mutate(annotation=paste0(code,"_", Type)) %>% select(-code)

primer.MLST.sl.table <- tibble(name=c("aroE_F", "aroE_R", 
                                      "dat_F", "dat_R", 
                                      "ddl_F", "ddl_R",
                                      "gmk_F", "gmk_R",
                                      "idh_F", "idh_R",
                                      "recA_F", "recA_R",
                                      "yqiL_F", "yqiL_R"),
                               sequence=c(aroE_F, aroE_R, 
                                          dat_F, dat_R, 
                                          ddl_F, ddl_R,
                                          gmk_F, gmk_R,
                                          idh_F, idh_R,
                                          recA_F, recA_R,
                                          yqiL_F, yqiL_R))

prefix.seq.table <- tibble(type=c("F", "R"),
                           prefix=c(paste0(BsaI.HF.recognized.sequence, "T"), paste0(BsaI.HF.recognized.sequence, "A")))

primer.MLST.sl.table %>% 
  left_join(coding.MLST.sl.table, by="name") %>% 
  tidyr::separate(., col=c(1), into=c("gene", "type"), sep="_") %>% 
  mutate(type) %>% mutate(annotation=paste0(code,"_", type)) %>% 
  left_join(., barcode.seq.annotation.table, by="annotation") %>% left_join(., prefix.seq.table, by="type") %>% 
  mutate(Seq=replace_na(Seq, "")) %>% 
  mutate(prefix=case_when(is.na(Type)~"",
                          TRUE ~ prefix)) %>% 
  mutate(final_seq=paste0(prefix, Seq, sequence)) %>% 
  mutate(name=paste0(gene,"_",annotation)) %>% 
  write_excel_csv(., file="Documents/2020_LinkouCGMH_CP/2022_研究/2022_盧教授研究計畫/20220711_Sl_primers.csv")
  
primer.MLST.sl.table %>% 
  left_join(coding.MLST.sl.table, by="name") %>% 
  tidyr::separate(., col=c(1), into=c("gene", "type"), sep="_") %>% 
  mutate(type) %>% mutate(annotation=paste0(code,"_", type)) %>% 
  left_join(., barcode.seq.annotation.table, by="annotation") %>% left_join(., prefix.seq.table, by="type") %>% 
  mutate(Seq=replace_na(Seq, "")) %>% 
  mutate(prefix=case_when(is.na(Type)~"",
                          TRUE ~ prefix)) %>% 
  mutate(final_seq=paste0(prefix, Seq, sequence)) %>% 
  mutate(name=paste0(gene,"_",annotation)) %>% View()





# Barcode 108 design ------------------------------------------------------

barcode_pos1 <- c("1","2","13")
barcode_pos2 <- c("3","4","14")
barcode_pos3_4 <- c("5-7", "6-8")
barcode_pos5 <- c("9", "10", "15")
barcode_pos6 <- c("11", "12")

combaination_sample_barcodes <- list(pos1  = barcode_pos1,
                                              pos2  = barcode_pos2,
                                              pos_34= barcode_pos3_4,
                                              pos_5 = barcode_pos5,
                                              pos6  = barcode_pos6)


<- combaination_sample_barcodes %>%
      cross() %>% 
      map(unlist) %>% 
      map(t) %>% 
      map(as_tibble) %>% 
      reduce(bind_rows)
