## UCR length distribution

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(dbplyr)
library(ggplot2)

UCR_location <- read.table(file = "01-data/UCR_raw/UCR_location_refseqid.txt", 
                           sep = "\t", header = TRUE)

lwd_pt <- .pt*72.27/96

p1 <- ggplot(data = UCR_location, mapping = aes(x = length)) +
  
  geom_histogram(binwidth = 20, 
                 fill = "gray", 
                 color = "#f1f1f1", 
                 linewidth = 0.5/lwd_pt) +
  
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_x_continuous(limits = c(170, 800), expand = c(0, 0)) +
  
  labs(x = "", 
       y = "Number of UCR") +
  
  theme_prism(base_size = 7, 
              base_family = "Arial", 
              base_fontface = "plain", 
              base_line_size = 0.5/lwd_pt, 
              base_rect_size = 0.5/lwd_pt, 
              border = TRUE) +
  
  theme(
    legend.position = "top", 
    legend.key.size = unit(0.25, "cm"), 
    aspect.ratio = 1
  )
p1

pdf(file = "03-results/15-UCR_length_distribution/UCR_length_distribution.pdf", 
    width = 480/254, height = 480/254)
p1
dev.off()




# 20240602 组合UCR长度分布和GC含量云雨图
library(cowplot)

# 把长度分布图作为p4
p4 <- ggplot(data = UCR_location, mapping = aes(x = length)) +
  
  geom_histogram(binwidth = 20, 
                 fill = "gray", 
                 color = "#f1f1f1", 
                 linewidth = 0.5/lwd_pt) +
  
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_x_continuous(limits = c(170, 800), expand = c(0, 0)) +
  
  labs(x = "", 
       y = "Number of UCR") +
  
  theme_prism(base_size = 7, 
              base_family = "Arial", 
              base_fontface = "plain", 
              base_line_size = 0.5/lwd_pt, 
              base_rect_size = 0.5/lwd_pt, 
              border = TRUE) +
  
  theme(
    legend.position = "top", 
    legend.key.size = unit(0.25, "cm"), 
    aspect.ratio = 1
  )
p4

pdf(file = "04-docs/fig3-3.pdf", width = 900/254, height = 600/254)
showtext_begin()

cowplot::plot_grid(p4, GC_raincloud_plot_3, 
                   align = c("hv"), axis = "bt",
                   nrow = 1, ncol = 2, 
                   rel_widths = 1, rel_heights = 1, 
                   labels = c("a", "b"), 
                   label_size = 8)

showtext_end()
dev.off()
