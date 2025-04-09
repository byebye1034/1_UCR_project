# GSE173754
# U937 vs PBMC

rm(list = ls())
setwd(dir = "D:/R_project/UCR_project/")

library(tidyverse)
library(ggplot2)
library(ggrepel)
lwd_pt <- .pt*72.27/96

GSE173754 <- read_tsv(file = "02-analysis/08-PBMC_U937_DEG/GSE173754.top.table.tsv")
GSE173754 <- GSE173754 %>% 
  mutate(Type = as.factor(ifelse(log2FoldChange > 1 & padj < 0.05, "up", ifelse(log2FoldChange < -1 & padj < 0.05, "down", "no change"))))

p <- ggplot(data = GSE173754, mapping = aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(color = "#f1f1f1", size = 0.7) +        # 点的大小0.7
  
  geom_hline(yintercept = c(-log10(0.05)),           # y轴显著的线
             linewidth = 0.5/lwd_pt, 
             lty = "dashed") +
  geom_vline(xintercept = c(-1, 1),                  # x轴显著的线
             linewidth = 0.5/lwd_pt, 
             lty = "dashed") +
  
  theme(
    panel.grid = element_blank(),                    # 背景网格删除
    panel.background = element_blank(),              # 背景删除
    text = element_text(size = 7),                   # 除坐标轴字体大小
    legend.position = "none",                        # 删除标签
    line = element_line(linewidth = 0.5/lwd_pt),     # 线条宽度（坐标轴上的点）
    
    axis.text = element_text(size = 7),              # 坐标轴字体大小
    axis.line = element_line(linewidth = 0.5/lwd_pt) # 坐标轴线条粗细
  )
p

gene_of_interest <- c("PBX3", "BCL11A", "LCOR", "MECOM", "ZFHX3")
gene_of_interest <- filter(GSE173754, Symbol %in% gene_of_interest)
p + geom_point(data = gene_of_interest, 
               mapping = aes(x = log2FoldChange, y = -log10(padj), fill = Symbol), 
               size = 0.7) +
  geom_text_repel(data = filter(GSE173754, Symbol %in% gene_of_interest$Symbol), 
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 50),
                  aes(label = Symbol), size = 3)


pdf(file = "03-results/figures/PBMC_U937_DEG.pdf", height = 480/254, width = 480/254)

p + geom_point(data = gene_of_interest, 
               mapping = aes(x = log2FoldChange, y = -log10(padj), fill = Symbol), 
               size = 0.7) +
  geom_text_repel(data = filter(GSE173754, Symbol %in% gene_of_interest$Symbol), 
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 50),
                  aes(label = Symbol), size = 3)

dev.off()






















