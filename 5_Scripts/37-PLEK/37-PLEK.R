# 对于lncRNA的研究，我想应该先预测一下蛋白编码能力

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(Biostrings)
library(stringr)

# get ncRNA UCR related ncRNA bed
ncRNA_UCR <- read.table(file = "02-analysis/16-New_classification/ncRNA_UCR(66)_ensembl_id.txt", sep = "\t")
ncRNA_UCR <- unique(ncRNA_UCR)

ncRNA_bed <- read.table(file = "01-data/37-PLEK/Gene_Non_Coding_RNA.bed", 
                        sep = "\t", header = FALSE)

UCR_related_ncRNA_bed <- merge(ncRNA_UCR, 
                               ncRNA_bed, 
                               by.x = "V1", 
                               by.y = "V4")
load("D:/R_project/UCR_project/01-data/reference/chromosomes_and_corresponding_nc_numbers.Rdata")
UCR_related_ncRNA_bed <- merge(UCR_related_ncRNA_bed, 
                               GRCh38.p14_report, 
                               by.x = "V1.y", 
                               by.y = "Chromosome.name")

write.table(UCR_related_ncRNA_bed[, c(7, 3:4, 2, 5:6)], 
            file = "01-data/37-PLEK/ncRNA_UCR_related_ncRNA.bed", 
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

# 文件太大PLEK会报错 -------------------------------------------------------------

# 文件太大PLEK会报错
# 将文件分成小块再用PLEK预测
ncRNA_fasta <- readDNAStringSet(filepath = "01-data/37-PLEK/ncRNA_UCR_related_ncRNA.fasta", 
                                format = "fasta")

# Define the number of sequences per file
sequences_per_file <- 1

# Calculate the number of files needed
num_files <- ceiling(length(ncRNA_fasta) / sequences_per_file)

# Loop through the sequences and write to new FASTA files
for (i in 1:num_files) {
  # Determine the start and end indices for the sequences to include in this file
  start_index <- (i - 1) * sequences_per_file + 1
  end_index <- min(i * sequences_per_file, length(ncRNA_fasta))
  
  # Extract the subset of sequences
  subset_sequences <- ncRNA_fasta[start_index:end_index]
  
  # Create the file name
  file_name <- paste0("D:/R_project/UCR_project/01-data/37-PLEK/fastafile/", i, ".fasta")
  
  # Write the subset to a new FASTA file
  writeXStringSet(subset_sequences, file_name)
}

# 读入预测结果 ------------------------------------------------------------------

predicted <- read.table(file = "01-data/37-PLEK/predicted/1_output", sep = "\t")

# 读入预测结果
for (i in c(3, 5:49, 51:57)) {
  # 构建文件路径
  filepath <- paste0("D:/R_project/UCR_project/01-data/37-PLEK/predicted/", i, "_output")
  
  # 读取数据
  data <- read.table(file = filepath, sep = "\t")
  predicted <- rbind(predicted, data)
}

# 将uc。4和uc.50分成上下两半进行预测，都是coding
# uc.4 ENSG00000237505 coding
# uc.50 ENSG00000271860 coding

predicted <- mutate(predicted, gene_id = V3)
predicted$gene_id <- str_sub(predicted$gene_id, 2, 16)
predicted <- predicted[, c(1, 4)]
colnames(predicted) <- c("coding_potential", "gene_id")

added <- data.frame(
  coding_potential = c("Coding", "Coding"), 
  gene_id = c("ENSG00000237505", "ENSG00000271860")
)

predicted <- rbind(predicted, added)
predicted <- mutate(predicted, software = "PLEK")

save(predicted, file = "02-analysis/37-PLEK/PLEK_predicted_results.Rdata")

PLEK <- predicted
table(PLEK$coding_potential)
# coding 50 non-coding 7

# 读取CPC2的预测结果 -------------------------------------------------------------

CPC2 <- read.table(file = "02-analysis/37-PLEK/result_cpc2.txt", sep = "\t")
CPC2 <- CPC2[, c(1, 7)]
CPC2$V1 <- str_sub(CPC2$V1, 1, 15)
CPC2 <- mutate(CPC2, software = "CPC2")
colnames(CPC2)[1:2] <- c("gene_id", "coding_potential")

table(CPC2$coding_potential)
# coding 39 non-coding 18

# 合并PLEK，CPC2结果 -----------------------------------------------------------

PLEK <- PLEK[,c(2, 1, 3)]
PLEK$coding_potential <- gsub("Coding", "coding", PLEK$coding_potential)
PLEK$coding_potential <- gsub("Non-coding", "noncoding", PLEK$coding_potential)

coding_pre <- rbind(PLEK, CPC2)
coding_pre$coding_potential <- factor(x = coding_pre$coding_potential, 
                                      levels = c("noncoding", "coding"))
coding_pre$software <- factor(x = coding_pre$software, 
                              levels = c("PLEK", "CPC2"))

save(coding_pre, file = "02-analysis/37-PLEK/coding_pre.Rdata")

library(ggplot2)
library(ggview)
library(scales)
library(ggprism)
library(sysfonts)
library(showtext)

font_add(family = "Arial", regular = "C:/Windows/Fonts/arial.ttf")
font_families()

lwd_pt <- .pt*72.27/96

coding_pre_p <- ggplot(data = coding_pre) +
  stat_count(mapping = aes(x = software, fill = coding_potential), 
             position = "fill", width = 0.6) +
  
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), 
                     labels = label_percent()) +
  scale_fill_manual(values = c("#00A087FF", "#3C5488FF"), 
                    name = "", 
                    breaks = c("noncoding", "coding")) +
  
  labs(x = "", 
       y = "Percent of UCR overlapping lncRNAs(%)") +
  
  theme_prism(base_size = 7, 
              base_family = "Arial", 
              base_fontface = "plain", 
              base_line_size = 0.5/lwd_pt, 
              base_rect_size = 0.5/lwd_pt, 
              border = TRUE) +
  
  theme(
    legend.position = "top", 
    legend.key.size = unit(0.25, "cm")
  )


coding_pre_p

ggview(coding_pre_p, width = 4.8, height = 6, units = "cm", dpi = 1200)

ggsave(filename = "03-results/37-PLEK/coding_pre_p.tiff", width = 4.8, height = 4.8, units = "cm", dpi = 1200)

pdf(file = "03-results/37-PLEK/coding_pre_p.pdf", width = 480/254, height = 600/254)

# 使用showtext渲染字体
showtext_begin()

ggplot(data = coding_pre) +
  stat_count(mapping = aes(x = software, fill = coding_potential), 
             position = "fill", width = 0.6) +
  
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), 
                     labels = label_percent()) +
  scale_fill_manual(values = c("#00A087FF", "#3C5488FF"), 
                    name = "", 
                    breaks = c("noncoding", "coding")) +
  
  labs(x = "", 
       y = "Percent of UCR overlapping lncRNAs(%)") +
  
  theme_prism(base_size = 7, 
              base_family = "Arial", 
              base_fontface = "plain", 
              base_line_size = 0.5/lwd_pt, 
              base_rect_size = 0.5/lwd_pt, 
              border = TRUE) +
  
  theme(
    legend.position = "top", 
    legend.key.size = unit(0.25, "cm")
  )

# 关闭showtext渲染字体
showtext_end()

dev.off()



























