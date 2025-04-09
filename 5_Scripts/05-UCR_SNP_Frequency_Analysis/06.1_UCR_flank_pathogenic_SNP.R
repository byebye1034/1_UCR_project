# UCR左翼和右翼，以及随机片段上pathogenic SNP的数量和长度的比例

rm(list = ls())
setwd(dir = "D:/R_project/UCR_project/")

library(tidyverse)
library(ggplot2)
library(reshape2)

# 准备UCR flank和random fragments的bed文件 --------------------------------------

# 准备UCR flank的bed文件
UCR_location <- read.table(file = "data/UCR_location.txt", sep = "\t", header = T)
UCR_left <- UCR_location[, c(1, 6, 7, 4, 5)]
UCR_right <- UCR_location[, c(1, 8, 9, 4, 5)]
write.table(UCR_left, file = "data/04-UCR_flank/UCR_left.bed", quote = F, 
            sep = "\t", col.names = F, row.names = F)
write.table(UCR_right, file = "data/04-UCR_flank/UCR_right.bed", quote = F, 
            sep = "\t", col.names = F, row.names = F)

# 准备random fragments的bed文件
random_fragments <- read.table(file = "data/random_fragments_location.txt", sep = "\t", header = T)
random_fragments_bed <- random_fragments[, c(5, 2, 3, 6, 4)]
write.table(random_fragments_bed, file = "data/04-UCR_flank/random_fragments.bed", quote = F, 
            sep = "\t", col.names = F, row.names = F)

# 在nebula利用bedtools获得UCR侧翼和随机片段上的pathogenic SNPs --------------------------

UCR_left_pathogenic_SNPs <- read.table(file = "data/04-UCR_flank/UCR_left_pathogenic_SNPs.bed", sep = "\t")
UCR_right_pathogenic_SNPs <- read.table(file = "data/04-UCR_flank/UCR_right_pathogenic_SNPs.bed", sep = "\t")
random_fragments_pathogenic_SNPs <- read.table(file = "data/04-UCR_flank/random_fragments_pathogenic_SNPs.bed", sep = "\t")

# 计算每个UCR的侧翼和随机片段上的pathogenic SNPs的位点长度占UCR长度的比例
UCR_left_pathogenic_SNPs_num <- UCR_left_pathogenic_SNPs %>% 
  group_by(V4) %>% 
  summarise(pathogenic_SNPs_position_num = n_distinct(V7))
colnames(UCR_left_pathogenic_SNPs_num)[1] <- "UCR_name"
UCR_left_pathogenic_SNPs_proportion <- merge(UCR_left, UCR_left_pathogenic_SNPs_num, by = "UCR_name", all.x = T)
UCR_left_pathogenic_SNPs_proportion[is.na(UCR_left_pathogenic_SNPs_proportion)] <- 0
UCR_left_pathogenic_SNPs_proportion <- UCR_left_pathogenic_SNPs_proportion %>% 
  mutate(proportion = pathogenic_SNPs_position_num/length)

UCR_right_pathogenic_SNPs_num <- UCR_right_pathogenic_SNPs %>% 
  group_by(V4) %>% 
  summarise(pathogenic_SNPs_position_num = n_distinct(V7))
colnames(UCR_right_pathogenic_SNPs_num)[1] <- "UCR_name"
UCR_right_pathogenic_SNPs_proportion <- merge(UCR_right, UCR_right_pathogenic_SNPs_num, by = "UCR_name", all.x = T)
UCR_right_pathogenic_SNPs_proportion[is.na(UCR_right_pathogenic_SNPs_proportion)] <- 0
UCR_right_pathogenic_SNPs_proportion <- UCR_right_pathogenic_SNPs_proportion %>% 
  mutate(proportion = pathogenic_SNPs_position_num/length)

random_fragments_pathogenic_SNPs_num <- random_fragments_pathogenic_SNPs %>% 
  group_by(V4) %>% 
  summarise(pathogenic_SNPs_position_num = n_distinct(V7))
colnames(random_fragments_pathogenic_SNPs_num)[1] <- "rf_name"
random_fragments_pathogenic_SNPs_proportion <- merge(random_fragments, random_fragments_pathogenic_SNPs_num, by = "rf_name", all.x = T)
random_fragments_pathogenic_SNPs_proportion[is.na(random_fragments_pathogenic_SNPs_proportion)] <- 0
random_fragments_pathogenic_SNPs_proportion <- random_fragments_pathogenic_SNPs_proportion %>% 
  mutate(proportion = pathogenic_SNPs_position_num/Length)

# UCR上pathogenic SNPs的长度占比
load("D:/R_project/UCR_project/data/0B-SNPAnalysisData/Patho_SNP_Percentage_in_UCR.Rdata")
UCR_info <- UCR_location[, c(4:5)]
UCR_pathogenic_SNPs_proportion <- merge(UCR_info, result, by = "UCR_name", all.x = T)
UCR_pathogenic_SNPs_proportion <- UCR_pathogenic_SNPs_proportion[, c(1, 6)]
UCR_pathogenic_SNPs_proportion[is.na(UCR_pathogenic_SNPs_proportion)] <- 0

# UCR/UCR_left/UCR_right/random_fragments的pathogenic SNPs的长度占比
pathogenic_SNP_proportion <- data.frame(
  UCR = UCR_pathogenic_SNPs_proportion$SNP_proportion, 
  UCR_left = UCR_left_pathogenic_SNPs_proportion$proportion, 
  UCR_right = UCR_right_pathogenic_SNPs_proportion$proportion, 
  random_fragments = random_fragments_pathogenic_SNPs_proportion$proportion
)

# 绘图展示
ggplot(data = long_data) +
  geom_count(mapping = aes(x = value, y = variable)) + facet_wrap(~ variable, nrow = 2)














