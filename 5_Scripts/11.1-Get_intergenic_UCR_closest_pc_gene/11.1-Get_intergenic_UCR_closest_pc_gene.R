# 用bedtools closest查找intergenic UCR最近的coding基因，做GO富集分析
# 没有考虑非编码基因是因为非编码不能做go分析，而且非编码也是通过查找target coding gene做go

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(ggplot2)

intergenic_coding_gene <- read.table(
  file = "01-data/11.1-Get_intergenic_UCR_closest_pc_gene/Intergenic_UCR_closest_protein_coding_gene.bed", 
  sep = "\t", header = FALSE
)
length(unique(intergenic_coding_gene$V8))   # 80个基因

write.table(
  intergenic_coding_gene$V8, 
  file = "02-analysis/11.1-Get_intergenic_UCR_closest_pc_gene/Intergenic_UCR_closest_protein_coding_gene.bed", 
  sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE
)

# get bed file and then do go enrichment analysis

# get distances between UCRs and protein-coding genes
intergenic_coding_gene <- intergenic_coding_gene %>% 
  mutate(distance1 = abs(V6 - V2)) %>% 
  mutate(distance2 = abs(V6 - V3)) %>% 
  mutate(distance = ifelse(distance1 > distance2, distance2, distance1)) %>% 
  mutate(log10distance = log10(distance))

lwd_pt <- .pt*72.27/96

p1 <- ggplot(
  data = intergenic_coding_gene, 
  mapping = aes(x = log10distance)
) +
  
  geom_density(linewidth = 0.5/lwd_pt, color = "#E64B35FF") +
  geom_rug(linewidth = 0.5/lwd_pt, color = "#E64B35FF") +
  
  geom_vline(xintercept = log10(2000), 
             linetype = "dashed", 
             color = "gray", 
             linewidth = 0.5/lwd_pt) +
  
  scale_y_continuous(expand = c(0, 0)) +
  
  labs(
    x = "log10(distance between intergenic
    UCR and protein coding genes)", 
    y = "density", 
    title = "Human"
  ) +
  
  theme(
    panel.background = element_blank(), 
    panel.grid = element_blank(), 
    line = element_line(linewidth = 0.5/lwd_pt), 
    
    axis.text = element_text(size = 5, color = "#000000"), 
    axis.ticks = element_line(linewidth = 0.5/lwd_pt), 
    axis.line = element_line(linewidth = 0.5/lwd_pt), 
    axis.title = element_text(size = 5, color = "#000000"), 
    
    plot.title = element_text(size = 5),
    
    aspect.ratio = 1:1
  )
p1   # 距离主要在10000，一万bp以上

# 只有1个在2000bp以下
# 只有13个在10000bp以下

pdf(
  file = "03-results/11.1-Get_intergenic_UCR_closest_pc_gene/distance_bt_intergenic_UCR_closest_pc_genes.pdf", 
  width = 480/254, height = 480/254
)
p1
dev.off()




# 20240501 将human和mouse的distance ditribution图整合
library(cowplot)

pdf(
  file = "03-results/11.1-Get_intergenic_UCR_closest_pc_gene/distance_bt_intergenic_UCR_closest_pc_genes_human_mouse.pdf", 
  width = 1100/254, height = 600/254
)

cowplot::plot_grid(p1, p2, 
                   labels = c("A", "B"), 
                   label_size = 8)

dev.off()

ggsave(
  filename = "03-results/11.1-Get_intergenic_UCR_closest_pc_gene/distance_bt_intergenic_UCR_closest_pc_genes_human_mouse.tiff", 
  width = 11, height = 6, units = "cm", dpi = 1200
)















