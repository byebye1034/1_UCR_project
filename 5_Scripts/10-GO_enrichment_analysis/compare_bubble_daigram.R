# 绘制将几组基因的富集分析结果进行比较的气泡图

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(readxl)
library(ggsci)

compare_GO <- read_xlsx(path = "02-analysis/10-GO_enrichment_analysis/compare_GO_table.xlsx", sheet = 1)

lwd_pt <- .pt*72.27/96

p1 <- ggplot(data = compare_GO, mapping = aes(x = UCR_type, y = Term)) +
  geom_point(mapping = aes(color = UCR_type, size = count)) +
  
  scale_color_npg() +
  
  
  
  theme(
    panel.background = element_blank(),                       # 删除背景
    text = element_text(size = 7),                            # 除坐标轴以外的字体大小
    line = element_line(linewidth = 0.5/lwd_pt),              # 坐标轴上的点的粗细
    
    axis.text = element_text(size = 7), 
    axis.line = element_line(linewidth = 0.5/lwd_pt),         # 坐标轴的粗细
    
    panel.grid = element_line(color = "gray", linewidth = 0.5/lwd_pt),          # 添加背景网格
    panel.border = element_rect(linewidth = 0.5/lwd_pt, fill = NA),             # 如果不设置fill=NA，气泡就显示不出来
    
    legend.key.height = unit(0.25, "cm"), 
    legend.justification = c(0, 0)
  )

p1

pdf(file = "03-results/10-GO_enrichment_analysis/compare_GO.pdf", width = 720/254, height = 480/254)
p1
dev.off()












