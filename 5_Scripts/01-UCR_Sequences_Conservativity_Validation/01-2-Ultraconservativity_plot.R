# 图3.1展示uc。18和uc。304在大鼠基因组里序列变化

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library(Biostrings)
library(ggmsa)
library(cowplot)

fai <- "D:/A_projects/UCR/图表/uc.304_rat.mas.fasta"
fasta <- readDNAMultipleAlignment(filepath = fai)

p304_1 <- ggmsa(msa = fasta, 
                start = 1, end = 80, 
                seq_name = TRUE)
p304_2 <- ggmsa(msa = fasta, 
                start = 81, end = 160, 
                seq_name = TRUE)
p304_3 <- ggmsa(msa = fasta, 
                start = 161, end = 240, 
                seq_name = TRUE)
p304_4 <- ggmsa(msa = fasta, 
                start = 241, end = 320, 
                seq_name = TRUE)

pdf(file = "03-results/01-UCR_Sequences_Conservativity_Validation/uc.304_81-160.pdf", 
    width = 1600/254, height = 600/254)

cowplot::plot_grid(p304_1, p304_2, p304_3, p304_4, 
                   align = c("none"), nrow = 4, ncol = 1, 
                   rel_widths = c(1, 1, 1, 1), 
                   rel_heights = c(1, 1, 1, 1))

dev.off()

# uc.18 -------------------------------------------------------------------

uc.18 <- readDNAMultipleAlignment(filepath = "D:/A_projects/UCR/图表/uc.18_rat.mas.fas", 
                                  format = "fasta")
p18_1 <- ggmsa(msa = uc.18, 
               start = 1, end = 82, 
               seq_name = TRUE)
p18_2 <- ggmsa(msa = uc.18, 
               start = 83, end = 164, 
               seq_name = TRUE)
p18_3 <- ggmsa(msa = uc.18, 
               start = 165, end = 246, 
               seq_name = TRUE)


pdf(file = "03-results/01-UCR_Sequences_Conservativity_Validation/uc.18.pdf", 
    width = 1600/254, height = 600/254)

cowplot::plot_grid(p18_1, p18_2, p18_3, 
                   align = c("none"), nrow = 3, ncol = 1)

dev.off()











