# exonic的go富集分析结果

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(readxl)
library(stringr)
library(dbplyr)

exonic_UCR <- read_xlsx(path = "01-data/10-GO_enrichment_analysis/exonic_GO_enrichment_analysis.xlsx", 
                        sheet = 2)
exonic_UCR$Term <- str_sub(exonic_UCR$Term, 12, nchar(exonic_UCR$Term))

# 获得-log10(FDR)
exonic_UCR <- exonic_UCR %>% 
  mutate(log10FDR = -log10(FDR))

# 按照GeneRatio进行排序
exonic_UCR <- exonic_UCR[order(exonic_UCR$GeneRatio, decreasing = FALSE), ]   # 升序后面ggplot2绘图才是降序
exonic_UCR$Term <- factor(exonic_UCR$Term, levels = exonic_UCR$Term, ordered = TRUE)

lwd_pt <- .pt*72.27/96

p1 <- ggplot(data = exonic_UCR, mapping = aes(x = GeneRatio, y = Term)) +
  geom_point(mapping = aes(size = Count, color = log10FDR)) +
  
  scale_x_continuous(limits = c(0, 35), expand = c(0, 0)) +
  scale_color_gradient(low = "#f1f1f1", high = "#E60212", name = "-log10(FDR)") +
  
  theme(
    panel.grid = element_blank(), 
    panel.background = element_blank(), 
    text = element_text(size = 7, color = "black"), 
    line = element_line(linewidth = 0.5/lwd_pt), 
    
    axis.text = element_text(size = 7, color = "black"), 
    axis.line = element_line(linewidth = 0.5/lwd_pt), 
    
    panel.border = element_rect(linewidth = 0.5/lwd_pt, fill = NA)   # 如果不设置fill=NA，气泡就显示不出来
  )
p1

pdf(file = "03-results/10-GO_enrichment_analysis/exonic_UCR_GO.pdf", width = 1600/254, height = 960/254)
p1
dev.off()






