# type III UCRs nearest PCGs GO
# human
# mouse

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library(stringr)
library(readxl)
library(ggplot2)
library(cowplot)
library(ggview)


# human -------------------------------------------------------------------

# human
human_GO_results <- read_xlsx(
  path = "03-results/47-Other_UCR_nearest_PCGs/humanTypeIIIUCRsNearestPCGsGO/metascape_result.xlsx", 
  sheet = 2
)
human_GO_results <- human_GO_results[str_detect(human_GO_results$GroupID, ".*Summary"), ]

# 排序
human_GO_results <- human_GO_results[order(human_GO_results$LogP, decreasing = TRUE), ]
human_GO_results$Description <- factor(
  x = human_GO_results$Description, 
  levels = human_GO_results$Description, 
  ordered = TRUE
)

lwd_pt <- .pt*72.27/96

# 绘制条形图
p1 <- ggplot(data = human_GO_results) +
  geom_bar(mapping = aes(x = -(LogP), y = Description), 
           stat = "identity", width = 0.8, 
           fill = "#E64B35FF", alpha = 0.5) +
  
  geom_text(aes(x = 0.1, y = Description, label = Description), 
            size = 5, size.unit = "pt", hjust = 0) +
  
  scale_x_continuous(limits = c(0, 15), expand = c(0, 0)) +
  
  labs(x = "-Log10 pvalue", 
       y = "GO biological process 
       (type III UCRs nearest PCGs)", 
       title = "") +
  
  theme(
    plot.title = element_text(size = 5),
    legend.key.height = unit(100/254, "cm"), # 设置图例的高度为2厘米
    legend.key.width = unit(100/254, "cm"),  # 设置图例的宽度为1厘米
    
    panel.grid = element_blank(), 
    panel.background = element_blank(), 
    text = element_text(size = 5), 
    line = element_line(linewidth = 0.5/lwd_pt), 
    
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(), 
    axis.line.y = element_line(linewidth = 0.5/lwd_pt, color = "black"), 
    
    axis.text.x = element_text(size = 5, color = "black"), 
    axis.line = element_line(linewidth = 0.5/lwd_pt, color = "black")
  )
p1

p1 + canvas(width = 8.25, height = 6, units = "cm", dpi = 1200)

# mouse -------------------------------------------------------------------

mouse_GO_results <- read_xlsx(
  path = "03-results/47-Other_UCR_nearest_PCGs/mouseTypeIIIUCRsNearestPCGsGO/metascape_result.xlsx", 
  sheet = 2
)
mouse_GO_results <- mouse_GO_results[str_detect(mouse_GO_results$GroupID, ".*Summary"), ]

mouse_GO_results <- mouse_GO_results[order(mouse_GO_results$LogP, decreasing = TRUE), ]
mouse_GO_results$Description <- factor(
  x = mouse_GO_results$Description, 
  levels = mouse_GO_results$Description, 
  ordered = TRUE
)

# get barplot
p2 <- ggplot(data = mouse_GO_results) +
  geom_bar(mapping = aes(x = -(LogP), y = Description), 
           stat = "identity", width = 0.8, 
           fill = "#4DBBD5FF", alpha = 0.5) +
  
  geom_text(aes(x = 0.1, y = Description, label = Description), 
            size = 5, size.unit = "pt", hjust = 0) +
  
  scale_x_continuous(limits = c(0, 12), expand = c(0, 0)) +
  
  labs(x = "-Log10 pvalue", 
       y = "GO biological process 
       (type III UCRs nearest PCGs)", 
       title = "") +
  
  theme(
    plot.title = element_text(size = 5),
    legend.key.height = unit(100/254, "cm"), # 设置图例的高度为2厘米
    legend.key.width = unit(100/254, "cm"),  # 设置图例的宽度为1厘米
    
    panel.grid = element_blank(), 
    panel.background = element_blank(), 
    text = element_text(size = 5), 
    line = element_line(linewidth = 0.5/lwd_pt), 
    
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(), 
    axis.line.y = element_line(linewidth = 0.5/lwd_pt, color = "black"), 
    
    axis.text.x = element_text(size = 5, color = "black"), 
    axis.line = element_line(linewidth = 0.5/lwd_pt, color = "black")
  )
p2

p2 + canvas(width = 8.25, height = 6, units = "cm", dpi = 1200)

# merge barplot
cowplot::plot_grid(p1, p2, 
                   nrow = 1, ncol = 2, 
                   align = c("h"), axis = "bt", rel_widths = c(1, 1), 
                   labels = c("c", "d"), 
                   label_size = 8)

pdf(
  file = "03-results/34-Intergenic_UCR_related_pc_GO/merged_intergenic_UCR_GO_barplot.pdf", 
  width = 1650/254, height = 600/254
)

cowplot::plot_grid(p1, p2, 
                   nrow = 1, ncol = 2, 
                   align = c("h"), axis = "bt", rel_widths = c(1, 1), 
                   labels = c("c", "d"), 
                   label_size = 8)

dev.off()




