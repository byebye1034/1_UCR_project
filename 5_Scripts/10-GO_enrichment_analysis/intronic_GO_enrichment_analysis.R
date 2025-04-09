# exonic的go富集分析结果
# barplot展示富集分析结果，网址如下
# https://mp.weixin.qq.com/s/fVIRX8ieyWRRVdArwOHQvg

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(readxl)
library(stringr)
library(dbplyr)

intronic_UCR <- read_xlsx(path = "01-data/10-GO_enrichment_analysis/intronic_GO_enrichment_analysis.xlsx", 
                        sheet = 5)
intronic_UCR$Term <- str_sub(intronic_UCR$Term, 12, nchar(intronic_UCR$Term))
intronic_UCR$Category <- str_sub(intronic_UCR$Category, 8, 9)

# 获得-log10(FDR)
intronic_UCR <- intronic_UCR %>% 
  mutate(log10FDR = -log10(FDR))

# 按照GeneRatio进行排序
intronic_UCR <- intronic_UCR[order(intronic_UCR$GeneRatio, decreasing = FALSE), ]   # 升序后面ggplot2绘图才是降序
intronic_UCR$Term <- factor(intronic_UCR$Term, levels = intronic_UCR$Term, ordered = TRUE)

lwd_pt <- .pt*72.27/96

p1 <- ggplot(data = intronic_UCR) +
  geom_point(mapping = aes(x = GeneRatio, y = Term, 
                           size = Count, color = log10FDR)) +
  
  scale_x_continuous(limits = c(0, 70), expand = c(0, 0)) +
  scale_color_gradient(low = "#f1f1f1", high = "#E60212", name = "-log10(FDR)") +
  
  theme(
    panel.grid = element_blank(), 
    panel.background = element_blank(), 
    text = element_text(size = 7, color = "black"), 
    line = element_line(linewidth = 0.5/lwd_pt), 
    
    axis.text = element_text(size = 7, color = "black"), 
    axis.line = element_line(linewidth = 0.5/lwd_pt), 
    
    panel.border = element_rect(linewidth = 0.5/lwd_pt, fill = NA)   # 如果不设置fill=NA，气泡就显示不出来
  ) +
  
  facet_grid(Category ~ ., 
             scales = "free",   # 通过添加 scales = "free"，每个分页都将根据自己的数据范围调整轴刻度，确保每个 Category 分页只包含相应的 Term。
             space = "free_y") +
  
  scale_y_discrete(labels = function(value) str_wrap(value, width = 50))
p1

pdf(file = "03-results/10-GO_enrichment_analysis/intronic_UCR_GO.pdf", width = 1450/254, height = 1080/254)
p1
dev.off()

p2 <- ggplot(data = intronic_UCR) +
  geom_bar(mapping = aes(x = GeneRatio, y = Term, fill = log10FDR), 
           stat = "identity", width = 0.8, alpha = 0.7) +
  
  geom_text(aes(x = 0.1, y = Term, label = Term), 
            size = 7, size.unit = "pt", hjust = 0) +
  
  scale_x_continuous(limits = c(0, 70), expand = c(0, 0)) +
  scale_fill_gradient(low = "#f1f1f1", high = "#E60212", name = "-log10(FDR)") +
  
  theme(
    panel.grid = element_blank(), 
    panel.background = element_blank(), 
    text = element_text(size = 7), 
    line = element_line(linewidth = 0.5/lwd_pt), 
    
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(), 
    axis.line.y = element_line(linewidth = 0.5/lwd_pt, color = "black"), 
    
    axis.text.x = element_text(size = 7, color = "black"), 
    axis.line = element_line(linewidth = 0.5/lwd_pt, color = "black")
  ) +
  
  facet_grid(
    Category ~ ., 
    scales = "free",   # 通过添加 scales = "free"，每个分页都将根据自己的数据范围调整轴刻度，确保每个 Category 分页只包含相应的 Term。
    space = "free_y"
  )
p2

pdf(file = "03-results/10-GO_enrichment_analysis/intronic_UCR_GO_barplot.pdf", 
    width = 1700/254, height = 960/254)
p2
dev.off()


























