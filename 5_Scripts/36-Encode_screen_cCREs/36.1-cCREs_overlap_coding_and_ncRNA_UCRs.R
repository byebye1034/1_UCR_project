# get and analyse how many cCREs overlap with coding UCRs and ncRNA UCRs

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library(stringr)
library(VennDiagram)

# cCREsOverlapCodingUCRs --------------------------------------------------

cCREsOverlapCodingUCRs <- read.table(
  file = "02-analysis/36-Encode_screen_cCREs/cCREsOverlapCoding.bed", 
  sep = "\t", header = FALSE)
length(unique(cCREsOverlapCodingUCRs$V4))
# 158 coding UCRs overlap with cCREs

# dELS
dELS <- cCREsOverlapCodingUCRs[str_detect(cCREsOverlapCodingUCRs$V10, "dELS"), ]
length(unique(dELS$V4)) # 117 coding UCRs
dELS_1 <- dELS$V4

# pELS
pELS <- cCREsOverlapCodingUCRs[str_detect(cCREsOverlapCodingUCRs$V10, "pELS"), ]
length(unique(pELS$V4)) # 22 coding UCRs
pELS_1 <- pELS$V4

# CTCF-bound
CTCF_bound <- cCREsOverlapCodingUCRs[str_detect(cCREsOverlapCodingUCRs$V10, "CTCF-bound"), ]
length(unique(CTCF_bound$V4)) # 79
CTCF_bound_1 <- CTCF_bound$V4

# DNase-H3K4me3
DNase_H3K4me3 <- cCREsOverlapCodingUCRs[str_detect(cCREsOverlapCodingUCRs$V10, "DNase-H3K4me3"), ]
length(unique(DNase_H3K4me3$V4)) # 3
DNase_H3K4me3_1 <- DNase_H3K4me3$V4

lwd_pt <- .pt*72.27/96

setwd(dir = "03-results/36-Encode_screen_cCREs/")
VennCodingUCRs <- venn.diagram(
  x = list(dELS_1, pELS_1, CTCF_bound_1, DNase_H3K4me3_1), 
  category.names = c("dELS", "pELS", "CTCF-bound", "DNase-H3K4me3"), 
  filename = NULL, 
  fill = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"), alpha = 0.5, 
  cex = 5/lwd_pt, cat.cex = 5/lwd_pt, cat.pos = c(-20, 20, -20, 20), unit = "px", 
  lwd = 0.5/lwd_pt
)

# venn diagram不能直接保存pdf，要设置上面：filename = NULL，再用下面的方法
pdf(file = "VennCodingUCRs.pdf", width = 480/254, height = 480/254)
grid.draw(VennCodingUCRs)
dev.off()

# cCREsOverlapNcRNAUCRs ---------------------------------------------------

rm(list = ls())
setwd(dir = "D:/R_project/UCR_project/")

cCREsOverlapNcRNAUCRs <- read.table(
  file = "02-analysis/36-Encode_screen_cCREs/cCREsOverlapNoncodingRNA.bed", 
  sep = "\t", header = FALSE)
length(unique(cCREsOverlapNcRNAUCRs$V4))
# 65 coding UCRs overlap with cCREs

# dELS
dELS <- cCREsOverlapNcRNAUCRs[str_detect(cCREsOverlapNcRNAUCRs$V10, "dELS"), ]
length(unique(dELS$V4)) # 33 nRNA UCRs
dELS_1 <- dELS$V4

# pELS
pELS <- cCREsOverlapNcRNAUCRs[str_detect(cCREsOverlapNcRNAUCRs$V10, "pELS"), ]
length(unique(pELS$V4)) # 18
pELS_1 <- pELS$V4

# CTCF-bound
CTCF_bound <- cCREsOverlapNcRNAUCRs[str_detect(cCREsOverlapNcRNAUCRs$V10, "CTCF-bound"), ]
length(unique(CTCF_bound$V4)) # 79
CTCF_bound_1 <- CTCF_bound$V4

# DNase-H3K4me3
DNase_H3K4me3 <- cCREsOverlapNcRNAUCRs[str_detect(cCREsOverlapNcRNAUCRs$V10, "DNase-H3K4me3"), ]
length(unique(DNase_H3K4me3$V4)) # 3
DNase_H3K4me3_1 <- DNase_H3K4me3$V4

lwd_pt <- .pt*72.27/96

setwd(dir = "03-results/36-Encode_screen_cCREs/")
VennNcRNAUCRs <- venn.diagram(
  x = list(dELS_1, pELS_1, CTCF_bound_1, DNase_H3K4me3_1), 
  category.names = c("dELS", "pELS", "CTCF-bound", "DNase-H3K4me3"), 
  filename = NULL, 
  fill = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"), alpha = 0.5, 
  cex = 5/lwd_pt, cat.cex = 5/lwd_pt, cat.pos = c(-20, 20, -20, 20), unit = "px", 
  lwd = 0.5/lwd_pt
)

# venn diagram不能直接保存pdf，要设置上面：filename = NULL，再用下面的方法
pdf(file = "VennNcRNAUCRs.pdf", width = 480/254, height = 480/254)
grid.draw(VennNcRNAUCRs)
dev.off()






















