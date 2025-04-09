# 用数据库的数据探索type I UCRs相关的可变剪接基因的表达模式

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(dbplyr)
library(pheatmap)
library(readxl)
library(stringr)
library(tidyverse)
library(biomaRt)
library(curl)
library(ggview)

## human brain rpkm
human_rpkm <- read.table(file = "01-data/18-Development_exp_pattern/Human_rpkm.txt", 
                         sep = " ", header = TRUE)
human_rpkm <- human_rpkm[, c(1:54)]

## 46 AS gens
typeIUCRsOverlapASGenes <- read.table(
  file = "03-results/53-TypeIUCRs_ASGenes/ASGenes.bed", sep = "\t"
)
expTypeIUCRsOverlapASGenes <- merge(
  x = typeIUCRsOverlapASGenes[, c(4:5)], y = human_rpkm, 
  by.x = "V4", by.y = "Names"
)
rownames(expTypeIUCRsOverlapASGenes) <- expTypeIUCRsOverlapASGenes[, 2]
expTypeIUCRsOverlapASGenes <- expTypeIUCRsOverlapASGenes[, c(3:55)]
expTypeIUCRsOverlapASGenes <- log10(expTypeIUCRsOverlapASGenes + 1)

lwd_pt <- .pt*72.27/96

min(expTypeIUCRsOverlapASGenes) # 0.055
max(expTypeIUCRsOverlapASGenes) # 2.46

# breaks
bk <- c(seq(0, 1.25, by = 0.01), seq(1.26, 2.5, by = 0.01))

p1 <- pheatmap(mat = expTypeIUCRsOverlapASGenes, 
               cellwidth = 5.5, cellheight = 3, 
               cluster_cols = FALSE, cluster_rows = TRUE, 
               show_rownames = FALSE, 
               angle_col = 90, main = "log10(rpkm + 1)", 
               color = c(colorRampPalette(colors = c("#083490","#f1f1f1"))(length(bk)/2), 
                         colorRampPalette(colors = c("#f1f1f1","#E60212"))(length(bk)/2)),
               legend_breaks=seq(0, 2.5, 0.5), 
               breaks = bk, fontsize = 6, 
               border = 0.5/lwd_pt, border_width = 0.5/lwd_pt)
p1

pdf(file = "03-results/18-Development_exp_pattern/HeatmapTypeIUCRsOverlapAsGenes.pdf", 
    width = 1440/254, height = 700/254)
p1
dev.off()





