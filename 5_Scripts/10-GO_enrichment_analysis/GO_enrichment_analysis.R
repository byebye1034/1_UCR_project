# 用intergenic/intronic/type II做GO富集分析

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(dbplyr)
library(readxl)
library(stringr)
library(ggprism)

# coding UCR GO -----------------------------------------------------------

coding_UCR <- read_xlsx(path = "01-data/10-GO_enrichment_analysis/coding_UCR_GO_results.xlsx", sheet = 2)
coding_UCR$Term <- str_sub(coding_UCR$Term, 12, nchar(coding_UCR$Term))
colnames(coding_UCR)[4] <- "GeneRatio"

# 获得-log10(FDR)
coding_UCR <- coding_UCR %>% 
  mutate(log10FDR = -log10(FDR))

# 将数据按照GeneRatio进行排序
coding_UCR <- coding_UCR[order(coding_UCR$GeneRatio, decreasing = FALSE), ]
coding_UCR$Term <- factor(x = coding_UCR$Term, levels = coding_UCR$Term, ordered = TRUE)

lwd_pt <- .pt*72.27/96

p1 <- ggplot(data = coding_UCR, mapping = aes(x = GeneRatio, y = Term)) +
  geom_point(mapping = aes(size = Count, color = log10FDR)) +
  
  scale_x_continuous(limits = c(0, 40), expand = c(0, 0)) +
  scale_color_gradient(low = "#f1f1f1", high = "#E60212", name = "-log10(FDR)") +
  
  theme(
    panel.grid = element_blank(), 
    panel.background = element_blank(), 
    text = element_text(size = 7, color = "black"), 
    line = element_line(linewidth = 0.5/lwd_pt), 
    
    axis.text = element_text(size = 7), 
    axis.line = element_line(linewidth = 0.5/lwd_pt), 
    
    panel.border = element_rect(linewidth = 0.5/lwd_pt, fill = NA)   # 如果不设置fill=NA，气泡就显示不出来
  )
p1

pdf(file = "03-results/10-GO_enrichment_analysis/coding_UCR_GO.pdf", width = 1600/254, height = 960/254)
p1
dev.off()

# lncRNA UCR GO ------------------------------------------------------------

ncRNA_UCR <- read_xlsx(path = "01-data/17-LncRNA_enrichment/GO Enrichment result download.xlsx", 
                       sheet = 2)
ncRNA_UCR <- ncRNA_UCR[c(1:20), ]

# 获得-log10(FDR)
ncRNA_UCR <- ncRNA_UCR %>% 
  mutate(log10FDR = -log10(p.adjust))

# 按照-log10(FDR)进行排序
ncRNA_UCR <- ncRNA_UCR[order(ncRNA_UCR$log10FDR, decreasing = TRUE), ]

# 根据-log10(FDR)选取前20
ncRNA_UCR <- ncRNA_UCR[order(ncRNA_UCR$log10FDR, decreasing = FALSE), ]
ncRNA_UCR$Description <- factor(x = ncRNA_UCR$Description, 
                                levels = ncRNA_UCR$Description, ordered = TRUE)

lwd_pt <- .pt*72.27/96

p1 <- ggplot(data = ncRNA_UCR, mapping = aes(x = Count, y = Description)) +
  geom_point(mapping = aes(size = Count, color = log10FDR)) +
  
  scale_x_continuous(limits = c(60, 115), expand = c(0, 0)) +
  scale_color_continuous(low = "#f1f1f1", high = "#E60212", name = "-log10(FDR)") +
  
  theme(
    panel.grid = element_blank(), 
    panel.background = element_blank(), 
    text = element_text(size = 7), 
    line = element_line(linewidth = 0.5/lwd_pt), 
    
    axis.text = element_text(size = 7, color = "black"), 
    axis.line = element_line(linewidth = 0.5/lwd_pt), 
    
    panel.border = element_rect(linewidth = 0.5/lwd_pt, fill = NA)   # 如果不设置fill=NA，气泡就显示不出来
  )
p1

pdf(file = "03-results/17-LncRNA_enrichment/ncRNA_UCR_GO.pdf", width = 1600/254, height = 960/254)
p1
dev.off()


