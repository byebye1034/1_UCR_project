# 看human intergenic ucr在mouse和rat里是不是也是intergenic，closest protein coding gene是不是也富集在brain development

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(ggplot2)

# get intergenic ucr
ucr_mouse <- read.table(
  file = "01-data/31-Intergenic_is_intergenic_in_mouse_rat/ucr_mouse.bed", 
  sep = "\t"
)
mouse_coding_ucr <- read.table(
  file = "01-data/31-Intergenic_is_intergenic_in_mouse_rat/mouse_coding_ucr.bed", 
  sep = "\t"
)
mouse_coding_ucr_bed <- mouse_coding_ucr[, c(1:4)]

mouse_ncRNA_ucr <- read.table(
  file = "01-data/31-Intergenic_is_intergenic_in_mouse_rat/Mus_ncRNA_UCR.bed", 
  sep = "\t"
)
mouse_ncRNA_ucr_bed <- mouse_ncRNA_ucr[, c(1:4)]

mouse_intergenic_ucr <- ucr_mouse %>% 
  anti_join(mouse_coding_ucr_bed) %>% 
  anti_join(mouse_ncRNA_ucr_bed)

write.table(
  mouse_intergenic_ucr, 
  file = "02-analysis/31-Intergenic_is_intergenic_in_mouse_rat/mouse_intergenic_ucr.bed", 
  sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE
)

# homology of intergenic between human and mouse
human_intergenic_ucr <- read.table(
  file = "02-analysis/16-New_classification/03-UCR_intergenic.bed", 
  sep = "\t"
)

homo_intergenic_ucr <- merge(human_intergenic_ucr, mouse_intergenic_ucr, by = "V4")
# human(99) mouse(122) 有82个是相同的

# then get closest protein coding gene in linux

mouse_intergenic_ucr_closest_protein_coding_gene <- read.table(
  file = "01-data/31-Intergenic_is_intergenic_in_mouse_rat/mouse_intergenic_ucr_closest_protein_coding_gene.bed", 
  sep = "\t"
)

# get ensembl id 
write.table(mouse_intergenic_ucr_closest_protein_coding_gene$V8, 
            file = "01-data/31-Intergenic_is_intergenic_in_mouse_rat/mouse_intergenic_ucr_closest_protein_coding_ensemblid.bed", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# plot distance distribution
mouse_intergenic_ucr_closest_protein_coding_gene <- mouse_intergenic_ucr_closest_protein_coding_gene %>% 
  mutate(distance1 = abs(V2 - V6)) %>% 
  mutate(distance2 = abs(V3 - V6)) %>% 
  mutate(distance = ifelse(distance1 > distance2, distance2, distance1))

lwd_pt <- .pt*72.27/96

p2 <- ggplot(data = mouse_intergenic_ucr_closest_protein_coding_gene, 
             mapping = aes(x = log10(distance))) +
  
  geom_density(linewidth = 0.5/lwd_pt, color = "#4DBBD5FF") +
  geom_rug(linewidth = 0.5/lwd_pt, color = "#4DBBD5FF") +
  
  geom_vline(
    xintercept = log10(2000), 
    linetype = "dashed", 
    color = "gray", 
    linewidth = 0.5/lwd_pt
  ) +
  
  scale_y_continuous(expand = c(0, 0)) +
  
  labs(x = "log10(distance between intergenic
       UCR and protein coding genes)", 
       y = "density", 
       title = "Mouse") +
  
  theme(
    panel.background = element_blank(), 
    panel.grid = element_blank(), 
    line = element_line(linewidth = 0.5/lwd_pt), 
    
    axis.text = element_text(size = 7, color = "#000000"), 
    axis.ticks = element_line(linewidth = 0.5/lwd_pt), 
    axis.line = element_line(linewidth = 0.5/lwd_pt), 
    axis.title = element_text(size = 7, color = "#000000"), 
    
    plot.title = element_text(size = 7),
    
    aspect.ratio = 1:1
  )
p2

pdf(
  file = "03-results/31-Intergenic_is_intergenic_in_mouse_rat/distance_inter_ucr_pc_genes.pdf", 
  width = 480/254, 
  height = 480/254
)
p2
dev.off()

#  获得UCR和相关蛋白基因的中间的位置的bed文件
mouse_intergenic_ucr_closest_protein_coding_gene <- mouse_intergenic_ucr_closest_protein_coding_gene %>% 
  mutate(before = V3 - V6) %>% 
  mutate(start = ifelse(
    before > 0, V6, V2
  )) %>% 
  mutate(
    end = ifelse(
      before > 0, V3, V7
    )
  )

mouse_intergenic_ucr_closest_protein_coding_gene$V1 <- paste("chr", mouse_intergenic_ucr_closest_protein_coding_gene$V1, sep = "")

write.table(
  mouse_intergenic_ucr_closest_protein_coding_gene[, c(1, 15:16, 9, 4)], 
  file = "02-analysis/31-Intergenic_is_intergenic_in_mouse_rat/mouse_location_bt_UCR_and_related_pc_genes.bed", 
  sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE
)
