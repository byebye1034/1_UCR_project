# coding UCR相关基因在brain当中有特殊pattern
# 利用数据库数据查看coding UCR相关基因在cerebellum、heart、kidney、liver、ovary、testis当中的表达

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(dplyr)
library(pheatmap)
library(readxl)
library(stringr)
library(tidyverse)

# human rpkm
human_rpkm <- read.table(file = "01-data/18-Development_exp_pattern/Human_rpkm.txt", 
                         sep = " ", header = TRUE)
colnames(human_rpkm)

cerebellum_human_rpkm <- human_rpkm[, c(1, 55:112)]
heart_human_rpkm <- human_rpkm[, c(1, 113:156)]
kidney_human_rpkm <- human_rpkm[, c(1, 157:192)]
liver_human_rpkm <- human_rpkm[, c(1, 193:241)]
overy_human_rpkm <- human_rpkm[, c(1, 242:259)]
testis_human_rpkm <- human_rpkm[, c(1, 260:298)]

# 获得coding UCR相关基因里富集到AS通路的基因
# 将coding UCR相关基因当中富集到alternative splicing通路的基因挑出来

AS_coding_UCR_gene <- read_xlsx(path = "01-data/10-GO_enrichment_analysis/coding_UCR_GO_results.xlsx", 
                                sheet = 2)
AS_coding_UCR_gene <- AS_coding_UCR_gene[, c(2, 6)]
AS_coding_UCR_gene$Term <- str_sub(AS_coding_UCR_gene$Term, 12, nchar(AS_coding_UCR_gene$Term))
AS_coding_UCR_gene <- AS_coding_UCR_gene[c(3, 7:10, 11), ]

# 分隔基因ID
gene_ids <- unlist(strsplit(AS_coding_UCR_gene$Genes, ","))
gene_ids <- unique(gene_ids)
gene_ids <- gsub(" ", "", gene_ids)

# 用这个coding UCR只是为了获得ensembl id和gene symbol的对应关系
load("D:/R_project/UCR_project/02-analysis/18-Development_dynamic_genes/coding_UCR_OrganDDGNum.Rdata") # 变量名叫coding UCR

AS_coding_UCR_gene <- filter(coding_UCR, V8 %in% gene_ids)
AS_coding_UCR_gene <- AS_coding_UCR_gene[!(AS_coding_UCR_gene$V9 == "."), ]
AS_coding_UCR_gene <- AS_coding_UCR_gene[, c(1:2)] # 富集到AS通路的基因

# cerebellum --------------------------------------------------------------

cerebellum_rpkm_AS_coding_UCR_gene <- merge(AS_coding_UCR_gene, cerebellum_human_rpkm, 
                                            by.x = "V8", by.y = "Names")
rownames(cerebellum_rpkm_AS_coding_UCR_gene) <- cerebellum_rpkm_AS_coding_UCR_gene$V9
cerebellum_rpkm_AS_coding_UCR_gene <- cerebellum_rpkm_AS_coding_UCR_gene[, c(-1, -2)]
heatmap_cerebellum_rpkm_AS_coding_UCR_gene <- log10(cerebellum_rpkm_AS_coding_UCR_gene + 1)
colnames(heatmap_cerebellum_rpkm_AS_coding_UCR_gene) <- 
  factor(x = colnames(heatmap_cerebellum_rpkm_AS_coding_UCR_gene), 
         levels = colnames(heatmap_cerebellum_rpkm_AS_coding_UCR_gene), 
         ordered = TRUE)

max(heatmap_cerebellum_rpkm_AS_coding_UCR_gene)
min(heatmap_cerebellum_rpkm_AS_coding_UCR_gene)

lwd_pt <- .pt*72.27/96

# breaks
bk <- c(seq(0, 1.24, by = 0.01), seq(1.25, 2.5, by = 0.01))

p1 <- pheatmap(mat = heatmap_cerebellum_rpkm_AS_coding_UCR_gene, 
               cellwidth = 5, cellheight = 5, 
               cluster_cols = FALSE, cluster_rows = TRUE, clustering_method = "average", 
               angle_col = 90, main = "log10(rpkm + 1)", 
               color = c(colorRampPalette(colors = c("#083490","#f1f1f1"))(length(bk)/2), 
                         colorRampPalette(colors = c("#f1f1f1","#E60212"))(length(bk)/2)),
               legend_breaks=seq(0, 2.5, 0.5), 
               breaks = bk, 
               fontsize = 5, 
               border = 0.5/lwd_pt)
p1

# heart -------------------------------------------------------------------

heart_rpkm_AS_coding_UCR_gene <- merge(AS_coding_UCR_gene, heart_human_rpkm, 
                                            by.x = "V8", by.y = "Names")
rownames(heart_rpkm_AS_coding_UCR_gene) <- heart_rpkm_AS_coding_UCR_gene$V9
heart_rpkm_AS_coding_UCR_gene <- heart_rpkm_AS_coding_UCR_gene[, c(-1, -2)]
heatmap_heart_rpkm_AS_coding_UCR_gene <- log10(heart_rpkm_AS_coding_UCR_gene + 1)
colnames(heatmap_heart_rpkm_AS_coding_UCR_gene) <- 
  factor(x = colnames(heatmap_heart_rpkm_AS_coding_UCR_gene), 
         levels = colnames(heatmap_heart_rpkm_AS_coding_UCR_gene), 
         ordered = TRUE)

max(heatmap_heart_rpkm_AS_coding_UCR_gene)
min(heatmap_heart_rpkm_AS_coding_UCR_gene)

lwd_pt <- .pt*72.27/96

# breaks
bk <- c(seq(0, 1.24, by = 0.01), seq(1.25, 2.5, by = 0.01))

p2 <- pheatmap(mat = heatmap_heart_rpkm_AS_coding_UCR_gene, 
               cellwidth = 5, cellheight = 5, 
               cluster_cols = FALSE, cluster_rows = TRUE, clustering_method = "average", 
               angle_col = 90, main = "log10(rpkm + 1)", 
               color = c(colorRampPalette(colors = c("#083490","#f1f1f1"))(length(bk)/2), 
                         colorRampPalette(colors = c("#f1f1f1","#E60212"))(length(bk)/2)),
               legend_breaks=seq(0, 2.5, 0.5), 
               breaks = bk, 
               fontsize = 5, 
               border = 0.5/lwd_pt)
p2

# kidney ------------------------------------------------------------------

kidney_rpkm_AS_coding_UCR_gene <- merge(AS_coding_UCR_gene, kidney_human_rpkm, 
                                       by.x = "V8", by.y = "Names")
rownames(kidney_rpkm_AS_coding_UCR_gene) <- kidney_rpkm_AS_coding_UCR_gene$V9
kidney_rpkm_AS_coding_UCR_gene <- kidney_rpkm_AS_coding_UCR_gene[, c(-1, -2)]
heatmap_kidney_rpkm_AS_coding_UCR_gene <- log10(kidney_rpkm_AS_coding_UCR_gene + 1)
colnames(heatmap_kidney_rpkm_AS_coding_UCR_gene) <- 
  factor(x = colnames(heatmap_kidney_rpkm_AS_coding_UCR_gene), 
         levels = colnames(heatmap_kidney_rpkm_AS_coding_UCR_gene), 
         ordered = TRUE)

max(heatmap_kidney_rpkm_AS_coding_UCR_gene)
min(heatmap_kidney_rpkm_AS_coding_UCR_gene)

lwd_pt <- .pt*72.27/96

# breaks
bk <- c(seq(0, 1.24, by = 0.01), seq(1.25, 2.5, by = 0.01))

p3 <- pheatmap(mat = heatmap_kidney_rpkm_AS_coding_UCR_gene, 
               cellwidth = 5, cellheight = 5, 
               cluster_cols = FALSE, cluster_rows = TRUE, clustering_method = "average", 
               angle_col = 90, main = "log10(rpkm + 1)", 
               color = c(colorRampPalette(colors = c("#083490","#f1f1f1"))(length(bk)/2), 
                         colorRampPalette(colors = c("#f1f1f1","#E60212"))(length(bk)/2)),
               legend_breaks=seq(0, 2.5, 0.5), 
               breaks = bk, 
               fontsize = 5, 
               border = 0.5/lwd_pt)
p3

# liver -------------------------------------------------------------------

liver_rpkm_AS_coding_UCR_gene <- merge(AS_coding_UCR_gene, liver_human_rpkm, 
                                        by.x = "V8", by.y = "Names")
rownames(liver_rpkm_AS_coding_UCR_gene) <- liver_rpkm_AS_coding_UCR_gene$V9
liver_rpkm_AS_coding_UCR_gene <- liver_rpkm_AS_coding_UCR_gene[, c(-1, -2)]
heatmap_liver_rpkm_AS_coding_UCR_gene <- log10(liver_rpkm_AS_coding_UCR_gene + 1)
colnames(heatmap_liver_rpkm_AS_coding_UCR_gene) <- 
  factor(x = colnames(heatmap_liver_rpkm_AS_coding_UCR_gene), 
         levels = colnames(heatmap_liver_rpkm_AS_coding_UCR_gene), 
         ordered = TRUE)

max(heatmap_liver_rpkm_AS_coding_UCR_gene)
min(heatmap_liver_rpkm_AS_coding_UCR_gene)

lwd_pt <- .pt*72.27/96

# breaks
bk <- c(seq(0, 1.24, by = 0.01), seq(1.25, 2.5, by = 0.01))

p4 <- pheatmap(mat = heatmap_liver_rpkm_AS_coding_UCR_gene, 
               cellwidth = 5, cellheight = 5, 
               cluster_cols = FALSE, cluster_rows = TRUE, clustering_method = "average", 
               angle_col = 90, main = "log10(rpkm + 1)", 
               color = c(colorRampPalette(colors = c("#083490","#f1f1f1"))(length(bk)/2), 
                         colorRampPalette(colors = c("#f1f1f1","#E60212"))(length(bk)/2)),
               legend_breaks=seq(0, 2.5, 0.5), 
               breaks = bk, 
               fontsize = 5, 
               border = 0.5/lwd_pt)
p4

# overy -------------------------------------------------------------------

overy_rpkm_AS_coding_UCR_gene <- merge(AS_coding_UCR_gene, overy_human_rpkm, 
                                        by.x = "V8", by.y = "Names")
rownames(overy_rpkm_AS_coding_UCR_gene) <- overy_rpkm_AS_coding_UCR_gene$V9
overy_rpkm_AS_coding_UCR_gene <- overy_rpkm_AS_coding_UCR_gene[, c(-1, -2)]
heatmap_overy_rpkm_AS_coding_UCR_gene <- log10(overy_rpkm_AS_coding_UCR_gene + 1)
colnames(heatmap_overy_rpkm_AS_coding_UCR_gene) <- 
  factor(x = colnames(heatmap_overy_rpkm_AS_coding_UCR_gene), 
         levels = colnames(heatmap_overy_rpkm_AS_coding_UCR_gene), 
         ordered = TRUE)

max(heatmap_overy_rpkm_AS_coding_UCR_gene)
min(heatmap_overy_rpkm_AS_coding_UCR_gene)

lwd_pt <- .pt*72.27/96

# breaks
bk <- c(seq(0, 1.24, by = 0.01), seq(1.25, 2.5, by = 0.01))

p5 <- pheatmap(mat = heatmap_overy_rpkm_AS_coding_UCR_gene, 
               cellwidth = 5, cellheight = 5, 
               cluster_cols = FALSE, cluster_rows = TRUE, clustering_method = "average", 
               angle_col = 90, main = "log10(rpkm + 1)", 
               color = c(colorRampPalette(colors = c("#083490","#f1f1f1"))(length(bk)/2), 
                         colorRampPalette(colors = c("#f1f1f1","#E60212"))(length(bk)/2)),
               legend_breaks=seq(0, 2.5, 0.5), 
               breaks = bk, 
               fontsize = 5, 
               border = 0.5/lwd_pt)
p5

# testis ------------------------------------------------------------------

testis_rpkm_AS_coding_UCR_gene <- merge(AS_coding_UCR_gene, testis_human_rpkm, 
                                        by.x = "V8", by.y = "Names")
rownames(testis_rpkm_AS_coding_UCR_gene) <- testis_rpkm_AS_coding_UCR_gene$V9
testis_rpkm_AS_coding_UCR_gene <- testis_rpkm_AS_coding_UCR_gene[, c(-1, -2)]
heatmap_testis_rpkm_AS_coding_UCR_gene <- log10(testis_rpkm_AS_coding_UCR_gene + 1)
colnames(heatmap_testis_rpkm_AS_coding_UCR_gene) <- 
  factor(x = colnames(heatmap_testis_rpkm_AS_coding_UCR_gene), 
         levels = colnames(heatmap_testis_rpkm_AS_coding_UCR_gene), 
         ordered = TRUE)

max(heatmap_testis_rpkm_AS_coding_UCR_gene)
min(heatmap_testis_rpkm_AS_coding_UCR_gene)

lwd_pt <- .pt*72.27/96

# breaks
bk <- c(seq(0, 1.24, by = 0.01), seq(1.25, 2.5, by = 0.01))

p6 <- pheatmap(mat = heatmap_testis_rpkm_AS_coding_UCR_gene, 
               cellwidth = 5, cellheight = 5, 
               cluster_cols = FALSE, cluster_rows = TRUE, clustering_method = "average", 
               angle_col = 90, main = "log10(rpkm + 1)", 
               color = c(colorRampPalette(colors = c("#083490","#f1f1f1"))(length(bk)/2), 
                         colorRampPalette(colors = c("#f1f1f1","#E60212"))(length(bk)/2)),
               legend_breaks=seq(0, 2.5, 0.5), 
               breaks = bk, 
               fontsize = 5, 
               border = 0.5/lwd_pt)
p6

# 组合图形 --------------------------------------------------------------------

# 非ggplot2绘制的图片不能用patchwork进行组合
# 攻略来自 https://www.jianshu.com/p/8fc823c39488

pdf(file = "03-results/20-Development_exp_pattern_other_organs/Development_exp_pattern_other_organs.pdf", 
    width = 12, height = 10)
cowplot::plot_grid(p1$gtable, p2$gtable, p3$gtable, p4$gtable, p5$gtable, p6$gtable, 
                   ncol = 2, labels = c("A.Cerebellum", 
                                        "B.Heart", 
                                        "C.Kidney", 
                                        "D.Liver", 
                                        "E.Overy", 
                                        "F.Testis"))
dev.off()



