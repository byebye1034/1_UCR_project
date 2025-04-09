# 查看NUCR-ncRNA在神经系统不同脑区的表达

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(pheatmap)

# get GTEx data
load("D:/R_project/UCR_project/02-analysis/22-Tissue_specificity_across_different_organs/GTEx_exp_tau_normalized.Rdata")
# get brain data
brain_expr <- GTEx_exp_tau_normalized[, c(8:21)]
brain_expr <- brain_expr %>% 
  mutate(Gene.ID = rownames(brain_expr))

# get lncRNA
ncRNA_UCR_lncRNA <- read.table(file = "02-analysis/16-New_classification/ncRNA_UCR(66)_ensembl_id.txt")
ncRNA_UCR_lncRNA <- unique(ncRNA_UCR_lncRNA)
colnames(ncRNA_UCR_lncRNA) <- "Gene.ID"

# lncRNA brain expr
brain_lncRNA_expr <- merge(ncRNA_UCR_lncRNA, brain_expr, by = "Gene.ID")
rownames(brain_lncRNA_expr) <- brain_lncRNA_expr[, 1]
brain_lncRNA_expr <- brain_lncRNA_expr[, -1]
brain_lncRNA_expr <- as.matrix(brain_lncRNA_expr)
brain_lncRNA_expr <- t(brain_lncRNA_expr)
brain_lncRNA_expr <- brain_lncRNA_expr[-14, ]

range(brain_lncRNA_expr)
# 0 1

# breaks
bk <- c(seq(0, 0.49, by = 0.01), seq(0.5, 1, by = 0.01))
pheatmap(mat = brain_lncRNA_expr, 
         cellwidth = 15, cellheight = 15, 
         cluster_rows = T, cluster_cols = T, 
         angle_col = 90, main = "log2(exp+1)", 
         color = c(colorRampPalette(colors = c("#083490","#f1f1f1"))(length(bk)/2), 
                   colorRampPalette(colors = c("#f1f1f1","#E60212"))(length(bk)/2)),
         legend_breaks=seq(0, 1, 0.2), 
         breaks = bk, 
         fontsize = 7, 
         display_numbers = TRUE)

# ?鿴һ??lncRNA????????֯???ı??? -----------------------------------------------------

load("D:/R_project/UCR_project/02-analysis/22-Tissue_specificity_across_different_organs/GTEx_exp_tau_normalized.Rdata")
GTEx_exp_tau_normalized <- GTEx_exp_tau_normalized[, c(1:54)] %>% 
  mutate(Gene.ID = rownames(GTEx_exp_tau_normalized))

# GTEx data
GTEx_exp <- read.table(file = "01-data/22-Tissue_specificity_across_different_organs/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", 
                       sep = "\t", skip = 2, header = TRUE)
GTEx_exp$Name <- str_sub(GTEx_exp$Name, 1, 15)
GTEx_exp <- GTEx_exp[, -2]

# get lncRNA
ncRNA_UCR_lncRNA <- read.table(file = "02-analysis/16-New_classification/ncRNA_UCR(66)_ensembl_id.txt")
ncRNA_UCR_lncRNA <- unique(ncRNA_UCR_lncRNA)
colnames(ncRNA_UCR_lncRNA) <- "Gene.ID"

# lncRNA expr
lncRNA_expr <- merge(ncRNA_UCR_lncRNA, GTEx_exp, by.x = "Gene.ID", by.y = "Name")
rownames(lncRNA_expr) <- lncRNA_expr[, 1]
lncRNA_expr <- lncRNA_expr[, -1]
lncRNA_expr <- as.matrix(lncRNA_expr)
lncRNA_expr <- t(lncRNA_expr)

# breaks
bk <- c(seq(0, 0.49, by = 0.01), seq(0.5, 1, by = 0.01))
pdf(file = "03-results/42-NUCR-ncRNA_expr/pheatmap.pdf", width = 1800/254, height = 2000/254)

# 99-106行

dev.off()

# 热图：去掉聚类树的同时保留聚类的顺序 ------------------------------------------------------

p1 <- pheatmap(mat = lncRNA_expr, 
               cellwidth = 5, cellheight = 5, 
               cluster_rows = T, cluster_cols = T, 
               angle_col = 90, main = "log2(exp+1)", 
               color = c(colorRampPalette(colors = c("#083490","#f1f1f1"))(length(bk)/2), 
                         colorRampPalette(colors = c("#f1f1f1","#E60212"))(length(bk)/2)),
               legend_breaks=seq(0, 1, 0.2), 
               breaks = bk, 
               fontsize = 5, 
               display_numbers = FALSE)

# 从聚类获得的结果图中获得聚类后的行和列的排序
gene_name <- colnames(lncRNA_expr)[p1$tree_col[["order"]]]
sample_name <- rownames(lncRNA_expr)[p1$tree_row[["order"]]]

# 重新排序行和列
new_lncRNA_expr<- lncRNA_expr[sample_name, gene_name]

# 绘制不使用聚类的图（在这里还变成使用默认配色）
pheatmap(mat = new_lncRNA_expr, 
         cellwidth = 5, cellheight = 5, 
         cluster_rows = F, cluster_cols = F, 
         angle_col = 90, main = "log2(exp+1)", 
         legend_breaks=seq(0, 1, 0.2), 
         breaks = bk,  
         fontsize = 5, 
         display_numbers = FALSE)
