# 查看UCR相关的coding gene在除人类以外的其他物种当中的exp pattern
# 绘制热图：主要关注出生前后AS基因的表达

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(pheatmap)

# rat ---------------------------------------------------------------------

# 先读入rat coding genes
rat_coding_UCR <- read.table(
  file = "02-analysis/28-Homology_rate_of_coding_ucr_related_genes/rat_coding_ucr.bed")

# get rat all genes rpkm
rat_rpkm <- read.table(file = "01-data/30-Exp_patterns_of_ucr_coding_genes_in_other_species/Rat_rpkm.txt", 
                       sep = " ", header = TRUE)
rat_brain_rpkm <- rat_rpkm[, c(1:66)]

rat_brain_ucr_coding_genes_rpkm <- merge(rat_coding_UCR[, c(8:9)], 
                                 rat_brain_rpkm, 
                                 by.x = "V8", by.y = "Names")
rat_brain_ucr_coding_genes_rpkm <- unique(rat_brain_ucr_coding_genes_rpkm)   # total 191个基因，rat_rpkm里有168个基因的数据
rat_brain_ucr_coding_genes_rpkm <- rat_brain_ucr_coding_genes_rpkm[, -1]
rownames(rat_brain_ucr_coding_genes_rpkm) <- rat_brain_ucr_coding_genes_rpkm[, 1]
rat_brain_ucr_coding_genes_rpkm <- rat_brain_ucr_coding_genes_rpkm[, -1]

# heatmap

lwd_pt <- .pt*72.27/96

heatmap_data <- log10(rat_brain_ucr_coding_genes_rpkm + 1)
max(heatmap_data) # 3.24926
min(heatmap_data) # 0

colnames(heatmap_data) <- factor(colnames(heatmap_data), levels = colnames(heatmap_data))

# breaks
bk <- c(seq(0, 1.62, by = 0.01), seq(1.63, 3.26, by = 0.01))

p1 <- pheatmap(mat = heatmap_data, 
               cellwidth = 10, cellheight = 10, 
               cluster_cols = FALSE, cluster_rows = TRUE, 
               angle_col = 90, main = "log10(rpkm + 1)", 
               color = c(colorRampPalette(colors = c("#083490","#f1f1f1"))(length(bk)/2), 
                         colorRampPalette(colors = c("#f1f1f1","#E60212"))(length(bk)/2)),
               legend_breaks=seq(0, 4, 0.4), 
               breaks = bk, 
               fontsize = 7)
p1

pdf(file = "03-results/30-Exp_patterns_of_ucr_coding_genes_in_other_species/rat_brain_ucr_coding_genes.pdf", 
    width = 10, height = 30)
p1
dev.off()

# mouse -------------------------------------------------------------------

rm(list = ls())

# get mouse coding genes
mouse_coding_genes <- read.table(
  file = "02-analysis/28-Homology_rate_of_coding_ucr_related_genes/mouse_coding_ucr.bed"
)

# get all mouse genes rpkm
mouse_rpkm <- read.table(file = "01-data/30-Exp_patterns_of_ucr_coding_genes_in_other_species/Mouse_rpkm.txt", 
                         sep = " ", header = TRUE)
mouse_brain_rpkm <- mouse_rpkm[, c(1:56)]

mouse_brain_ucr_coding_genes_rpkm <- merge(mouse_coding_genes[, c(8:9)], 
                                           mouse_brain_rpkm, 
                                           by.x = "V8", by.y = "Names")

mouse_brain_ucr_coding_genes_rpkm <- unique(mouse_brain_ucr_coding_genes_rpkm) # 204个基因，mouse rpkm里有202个数据
rownames(mouse_brain_ucr_coding_genes_rpkm) <- mouse_brain_ucr_coding_genes_rpkm[, 2]
mouse_brain_ucr_coding_genes_rpkm <- mouse_brain_ucr_coding_genes_rpkm[, c(-1, -2)]

# heatmap
lwd_pt <- .pt*72.27/96

heatmap_data_mouse <- log10(mouse_brain_ucr_coding_genes_rpkm + 1)
save(heatmap_data_mouse, file = "02-analysis/30-Exp_patterns_of_ucr_coding_genes_in_other_species/heatmap_data_mouse.Rdata")
max(heatmap_data_mouse) # 2.830871
min(heatmap_data_mouse) # 0

colnames(heatmap_data) <- factor(colnames(heatmap_data), levels = colnames(heatmap_data))

# breaks
bk <- c(seq(0, 1.42, by = 0.01), seq(1.43, 2.84, by = 0.01))

p2 <- pheatmap(mat = heatmap_data_mouse, 
               cellwidth = 10, cellheight = 10, 
               cluster_cols = FALSE, cluster_rows = TRUE, 
               angle_col = 90, main = "log10(rpkm + 1)", 
               color = c(colorRampPalette(colors = c("#083490","#f1f1f1"))(length(bk)/2), 
                         colorRampPalette(colors = c("#f1f1f1","#E60212"))(length(bk)/2)),
               legend_breaks=seq(0, 4, 0.4), 
               breaks = bk, 
               fontsize = 7)
p2

pdf(file = "03-results/30-Exp_patterns_of_ucr_coding_genes_in_other_species/mouse_brain_ucr_coding_genes.pdf", 
    width = 10, height = 30)
p2
dev.off()

































