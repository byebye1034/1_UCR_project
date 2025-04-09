# 利用SEEKR探索ncRNA UCR相关的lncRNA的功能

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(Biostrings)

# get bed file
lncRNA_with_name <- read.table(
  file = "01-data/33-SEEKR_exploration_of_lncrna_function/Homo_lncRNA_with_gene_name_1.bed", 
  sep = "\t"
) # 5757
lncRNA_without_name <- read.table(
  file = "01-data/33-SEEKR_exploration_of_lncrna_function/Homo_lncRNA_without_gene_name_1.bed", 
  sep = "\t"
) # 13109

# exchange chr id
load("D:/R_project/UCR_project/01-data/reference/chromosomes_and_corresponding_nc_numbers.Rdata")

lncRNA_with_name <- merge(GRCh38.p14_report, lncRNA_with_name, 
                          by.x = "Chromosome.name", 
                          by.y = "V1")
lncRNA_with_name <- lncRNA_with_name[, c(2:8)]

lncRNA_without_name <- merge(GRCh38.p14_report, lncRNA_without_name, 
                             by.x = "Chromosome.name", 
                             by.y = "V1")
lncRNA_without_name <- lncRNA_without_name[, c(2:7)]

write.table(lncRNA_with_name[, c(1:5)], 
            file = "02-analysis/33-SEEKR_exploration_of_lncrna_function/lncRNA_with_name_chraccid.bed", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(lncRNA_without_name[, c(1:5)], 
            file = "02-analysis/33-SEEKR_exploration_of_lncrna_function/lncRNA_without_name_chraccid.bed", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# then get fasta in linux

rm(list = ls())

# get all lncRNA sequence
lncRNA_with_name_seq <- readDNAStringSet(
  filepath = "02-analysis/33-SEEKR_exploration_of_lncrna_function/lncRNA_with_name.fasta"
)

lncRNA_without_name_seq <- readDNAStringSet(
  filepath = "02-analysis/33-SEEKR_exploration_of_lncrna_function/lncRNA_without_name.fasta"
)

names(lncRNA_with_name_seq) <- lncRNA_with_name$V5
names(lncRNA_without_name_seq) <- lncRNA_without_name$V5

writeXStringSet(lncRNA_with_name_seq, 
                "02-analysis/33-SEEKR_exploration_of_lncrna_function/clear_lncRNA_with_name.seq")
writeXStringSet(lncRNA_without_name_seq, 
                "02-analysis/33-SEEKR_exploration_of_lncrna_function/clear_lncRNA_without_name.seq")

# get ncRNA UCR related lncRNA sequence
rm(list = ls())

ncRNA_UCR <- read.table(
  file = "02-analysis/16-New_classification/02-UCR_from_Non_Coding_RNA(66).bed", 
  sep = "\t"
)

lncRNA_seq <- readDNAStringSet(
  filepath = "02-analysis/33-SEEKR_exploration_of_lncrna_function/clear_all_lncRNA(with_and_without_name).seq", 
  format = "fasta"
)

# 提取数据框第八列的名称
target_names <- ncRNA_UCR[, 8]

# 根据名称筛选出符合条件的序列名称
selected_names <- names(lncRNA_seq) %in% target_names

# 从原始 DNAStringSet 中提取符合条件的序列
ncRNA_UCR_lncRNA_seq <- lncRNA_seq[selected_names]

# 保存结果到新的 fasta 文件中
writeXStringSet(ncRNA_UCR_lncRNA_seq, 
                "02-analysis/33-SEEKR_exploration_of_lncrna_function/ncRNA_UCR_lncRNA_seq.fasta")

















