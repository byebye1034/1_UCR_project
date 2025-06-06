setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(ggsci)

LCOR_pathogeicity <- read.table(file = "01-data/snp_pathogenicity_scatter_plot/LCOR_pathogenicity.bed")
LCOR_pathogeicity <- LCOR_pathogeicity[, c(5, 13:14)]
LCOR_pathogeicity <- LCOR_pathogeicity[!duplicated(LCOR_pathogeicity), ]
# is.unsorted(LCOR_pathogeicity$V5)

position <- table(LCOR_pathogeicity$V5)
LCOR_position <- rep(c(1:1085), position)
LCOR_pathogeicity <- LCOR_pathogeicity %>% 
  mutate(position = LCOR_position)

lwd_pt <- .pt*72.27/96

ggplot(data = LCOR_pathogeicity, mapping = aes(x = position, y = V13)) +
  geom_point(mapping = aes(color = V14), size = 0.2) +
  
  scale_color_manual(values = c("#f1f1f1", "#4DBBD5FF", "#E64B35FF"))+
  
  geom_hline(yintercept = c(0.35, 0.55), 
             linewidth = 0.5/lwd_pt, 
             lty = "dashed") +
  
  geom_vline(xintercept = c(894, 1080), 
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

pdf(file = "03-results/figures/LCOR_pathogenicity_pointplot.pdf", width = 720/254, height = 480/254)
dev.off()





