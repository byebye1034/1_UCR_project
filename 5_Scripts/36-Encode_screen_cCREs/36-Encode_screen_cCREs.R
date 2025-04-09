# 统计intergenic UCR主要和哪些cCREs重叠

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(stringr)
library(VennDiagram)

cCREs_overlap_intergenic_UCR <- read.table(
  file = "02-analysis/36-Encode_screen_cCREs/cCREs_overlap_intergenic.bed", 
  sep = "\t", header = FALSE
)
length(unique(cCREs_overlap_intergenic_UCR$V4))
cCREs_1 <- unique(cCREs_overlap_intergenic_UCR$V4)

# other UCR
other_UCR <- read.table(file = "02-analysis/16-New_classification/03-UCR_intergenic.bed")
length(other_UCR$V4)
other_UCR_1 <- other_UCR$V4

# dELS
dELS <- cCREs_overlap_intergenic_UCR[str_detect(cCREs_overlap_intergenic_UCR$V10, "dELS"), ]
length(unique(dELS$V4)) # 40个UCR
dELS_1 <- dELS$V4

# pELS
pELS <- cCREs_overlap_intergenic_UCR[str_detect(cCREs_overlap_intergenic_UCR$V10, "pELS"), ]
length(unique(pELS$V4)) # 3个UCR
pELS_1 <- pELS$V4

# CTCF-bound
CTCF_bound <- cCREs_overlap_intergenic_UCR[str_detect(cCREs_overlap_intergenic_UCR$V10, "CTCF-bound"), ]
length(unique(CTCF_bound$V4)) # 28个UCR
CTCF_bound_1 <- CTCF_bound$V4

# DNase-H3K4me3
DNase_H3K4me3 <- cCREs_overlap_intergenic_UCR[str_detect(cCREs_overlap_intergenic_UCR$V10, "DNase-H3K4me3"), ]
length(unique(DNase_H3K4me3$V4)) # 4个UCR
DNase_H3K4me3_1 <- DNase_H3K4me3$V4

lwd_pt <- .pt*72.27/96

setwd(dir = "03-results/36-Encode_screen_cCREs/")
venn1 <- venn.diagram(
  x = list(dELS_1, pELS_1, CTCF_bound_1, DNase_H3K4me3_1), 
  category.names = c("dELS", "pELS", "CTCF-bound", "DNase-H3K4me3"), 
  filename = NULL, 
  fill = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"), alpha = 0.5, 
  cex = 5/lwd_pt, cat.cex = 5/lwd_pt, cat.pos = c(-20, 20, -20, 20), unit = "px", 
  lwd = 0.5/lwd_pt
)

# venn diagram不能直接保存pdf，要设置上面：filename = NULL，再用下面的方法
pdf(file = "venn.pdf", width = 480/254, height = 480/254)
grid.draw(venn1)
dev.off()

setwd(dir = "03-results/36-Encode_screen_cCREs/")
venn2 <- venn.diagram(
  x = list(other_UCR_1, cCREs_1), 
  category.names = c("other UCR", "cCREs"), 
  filename = NULL, 
  fill = c("#ffffff", "gray20"), alpha = 0.5, 
  cex = 5/lwd_pt, cat.cex = 5/lwd_pt, cat.pos = c(-20, 20), unit = "px", 
  lwd = 0.5/lwd_pt
)

pdf(file = "venn_2.pdf", width = 480/254, height = 480/254)
grid.draw(venn2)
dev.off()























