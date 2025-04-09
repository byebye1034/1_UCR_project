# 用metascape网站做的GO富集分析结果绘图

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(readxl)

GO_result <- read_xlsx(path = "03-results/10-GO_enrichment_analysis/meatscape_coding_UCR_genes_GO/metascape_result.xlsx", 
                       sheet = 2)
GO_result <- GO_result[str_detect(GO_result$GroupID, ".*Summary"), ]

# 排序
GO_result <- GO_result[order(GO_result$LogP, decreasing = TRUE), ]
GO_result$Description <- factor(x = GO_result$Description, 
                                levels = GO_result$Description, ordered = TRUE)

lwd_pt <- .pt*72.27/96

p1 <- ggplot(data = GO_result) +
  geom_bar(mapping = aes(x = -(LogP), y = Description, fill = -(LogP)), 
           stat = "identity", width = 0.8) +
  
  geom_text(aes(x = 0.1, y = Description, label = Description), 
            size = 7, size.unit = "pt", hjust = 0) +
  
  scale_x_continuous(limits = c(0, 25), expand = c(0, 0)) +
  scale_fill_gradient(low = "#f1f1f1", high = "#E60212", name = "-log10(FDR)") +
  
  labs(x = "-log10(p-value)", 
       y = "", 
       title = "The Top 20 Enriched GO Terms(coding UCR related protein coding genes)") +
  
  theme(
    plot.title = element_text(size = 7),
    legend.key.height = unit(100/254, "cm"), # 设置图例的高度为2厘米
    legend.key.width = unit(100/254, "cm"),  # 设置图例的宽度为1厘米
    
    panel.grid = element_blank(), 
    panel.background = element_blank(), 
    text = element_text(size = 7), 
    line = element_line(linewidth = 0.5/lwd_pt), 
    
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(), 
    axis.line.y = element_line(linewidth = 0.5/lwd_pt, color = "black"), 
    
    axis.text.x = element_text(size = 7, color = "black"), 
    axis.line = element_line(linewidth = 0.5/lwd_pt, color = "black")
  )
p1

pdf(file = "03-results/10-GO_enrichment_analysis/metascape_coding_UCR_GO_barplot.pdf", 
    width = 960/254, height = 1200/254)
p1
dev.off()

ggsave(filename = "03-results/10-GO_enrichment_analysis/metascape_coding_UCR_GO_barplot.tiff", 
       width = 9.6, height = 12, units = "cm", dpi = 1200)

p1 <- ggplot(data = GO_result) +
  geom_bar(mapping = aes(x = -(LogP), y = Description), 
           stat = "identity", width = 0.8, 
           fill = "#3C5488FF", alpha = 0.5) +
  
  geom_text(aes(x = 0.1, y = Description, label = Description), 
            size = 5, size.unit = "pt", hjust = 0) +
  
  scale_x_continuous(limits = c(0, 25), expand = c(0, 0)) +
  
  labs(x = "-Log10 pvalue", 
       y = "GO biological process (coding UCR related protein coding genes)", 
       title = "") +
  
  theme(
    plot.title = element_text(size = 5),
    
    panel.grid = element_blank(), 
    panel.background = element_blank(), 
    text = element_text(size = 5), 
    line = element_line(linewidth = 0.5/lwd_pt), 
    
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(), 
    axis.line.y = element_line(linewidth = 0.5/lwd_pt, color = "black"), 
    
    axis.text.x = element_text(size = 5, color = "black"), 
    axis.line = element_line(linewidth = 0.5/lwd_pt, color = "black"), 
    
    plot.margin = unit(c(0,0.2,0,0.2), "cm")
  )
p1

ggview(p1, width = 8.25, height = 6, units = "cm", dpi = 1200)

ggsave(filename = "03-results/10-GO_enrichment_analysis/metascape_coding_UCR_GO_barplot.tiff", 
       width = 8.25, height = 6, units = "cm", dpi = 1200)

pdf(file = "03-results/10-GO_enrichment_analysis/metascape_coding_UCR_GO_barplot.pdf", 
    width = 825/254, height = 600/254)
p1
dev.off()








