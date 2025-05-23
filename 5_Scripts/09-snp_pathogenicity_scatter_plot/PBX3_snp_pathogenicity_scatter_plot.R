setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(ggplot2)

PBX3_exon_pathogenicity <- read.table(file = "01-data/snp_pathogenicity_scatter_plot/PBX3_pathogenicity.bed")
PBX3_exon_pathogenicity <- PBX3_exon_pathogenicity[, c(8, 15:16)]
PBX3_exon_pathogenicity <- PBX3_exon_pathogenicity[!duplicated(PBX3_exon_pathogenicity), ]
# is.unsorted(PBX3_exon_pathogenicity$V8)����������

position <- table(PBX3_exon_pathogenicity$V8)
PBX3_position <- rep(c(1:1089), position)
PBX3_exon_pathogenicity <- PBX3_exon_pathogenicity %>% 
  mutate(position = PBX3_position)

lwd_pt <- .pt*72.27/96

ggplot(data = PBX3_exon_pathogenicity, mapping = aes(x = position, y = V15)) +
  geom_point(mapping = aes(color = V16), size = 0.2) +
  
  scale_color_manual(values = c("#f1f1f1", "#4DBBD5FF", "#E64B35FF")) +
  
  geom_hline(yintercept = c(0.35, 0.55), 
             linewidth = 0.5/lwd_pt, 
             lty = "dashed") +
  
  geom_vline(xintercept = c(268, 431), 
             linewidth = 0.5/lwd_pt, 
             lty = "dashed") +
  theme(
    panel.grid = element_blank(),                    # ��������ɾ��
    panel.background = element_blank(),              # ����ɾ��
    text = element_text(size = 7),                   # �������������С
    legend.position = "none",                        # ɾ����ǩ
    line = element_line(linewidth = 0.5/lwd_pt),     # �������ȣ��������ϵĵ㣩
    
    axis.text = element_text(size = 7),              # �����������С
    axis.line = element_line(linewidth = 0.5/lwd_pt) # ������������ϸ
  )

pdf(file = "03-results/figures/PBX3_pathogenicity_pointplot.pdf", width = 720/254, height = 480/254)
dev.off()





