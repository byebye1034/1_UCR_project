# UCR相关的lncRNA的相互作用数目

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library(BiocManager)
library(ggplot2)
library(stringr)
library(ggbeeswarm)
library(ggsci)
library(ggpubr)
library(Hmisc)
library(ggview)
library(gginnards)
library(scales)
library(cowplot)

# get ucr related lncRNA
ncRNA_UCR_lncRNA <- read.table(
  file = "02-analysis/16-New_classification/02-UCR_from_Non_Coding_RNA.bed", 
  sep = "\t"
) # 100个ucr对应着86个lncRNA 有58个“.”，53个lncRNA有名字
ncRNA_UCR_lncRNA <- ncRNA_UCR_lncRNA[, c(8, 9)] %>% unique()

# get ncRNA RF related lncRNA
ncRNA_RF_lncRNA <- read.table(file = "02-analysis/16-New_classification/02-rf_from_Non_Coding_RNA.bed", 
                              sep = "\t")
ncRNA_RF_lncRNA <- ncRNA_RF_lncRNA[, c(8, 9)] %>% unique()

# get all lncRNA-miRNA interaction
high <- read.table(file = "01-data/25.2-LncRNA_interaction_num/mircode_highconsfamilies.txt", 
                   sep = "\t", header = TRUE)
med <- read.table(file = "01-data/25.2-LncRNA_interaction_num/mircode_medconsfamilies.txt", 
                  sep = "\t", header = TRUE)
lncRNA_miRNA_inter <- rbind(high, med)
lncRNA_miRNA_inter <- lncRNA_miRNA_inter[, c(1:4)]
lncRNA_miRNA_inter$gene_id <- str_sub(lncRNA_miRNA_inter$gene_id, 1, 15)

# ncRNA_UCR_lncRNA interaction
ncRNA_UCR_lncRNA_inter <- lncRNA_miRNA_inter %>% 
  filter(gene_id %in% ncRNA_UCR_lncRNA$V8) %>% 
  unique() %>% 
  group_by(gene_id) %>% 
  summarise(Count = n()) %>% 
  mutate(type = "UCR")

# ncRNA_RF_lncRNA interaction
ncRNA_RF_lncRNA_inter <- lncRNA_miRNA_inter %>% 
  filter(gene_id %in% ncRNA_RF_lncRNA$V8) %>% 
  unique() %>% 
  group_by(gene_id) %>% 
  summarise(Count = n()) %>% 
  mutate(type = "RF")

wilcox.test(ncRNA_UCR_lncRNA_inter$Count, ncRNA_RF_lncRNA_inter$Count)
t.test(ncRNA_UCR_lncRNA_inter$Count, ncRNA_RF_lncRNA_inter$Count)

ncRNA_UCR_RF_lncRNA_miRNA_inter <- rbind(ncRNA_UCR_lncRNA_inter, ncRNA_RF_lncRNA_inter)
save(ncRNA_UCR_RF_lncRNA_miRNA_inter, 
     file = "02-analysis/25.2-LncRNA_interaction_num/ncRNA_UCR_RF_lncRNA_miRNA_inter.Rdata")

load(file = "02-analysis/25.2-LncRNA_interaction_num/ncRNA_UCR_RF_lncRNA_miRNA_inter.Rdata")

ncRNA_UCR_RF_lncRNA_miRNA_inter$type <- factor(x = ncRNA_UCR_RF_lncRNA_miRNA_inter$type, 
                                               levels = c("UCR", "RF"), 
                                               ordered = TRUE)

p_lncRNA_inter <- ggplot(data = ncRNA_UCR_RF_lncRNA_miRNA_inter, 
                         mapping = aes(x = type, y = Count)) +
  
  geom_beeswarm(mapping = aes(color = type, fill = type), 
                size = 0.2, cex = 2) +
  
  scale_color_npg() +
  scale_fill_npg() +
  scale_x_discrete(labels = c("UCR" = "UCR (n = 62)", 
                              "RF" = "Randf (n = 47)")) +
  
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  
  labs(x = "", 
       y = "No. of lncRNA-miRNA interactions")  +
  
  theme_prism(palette = "black_and_white", 
              base_size = 7, 
              base_family = "Arial", 
              base_fontface = "plain", 
              base_line_size = 0.5/lwd_pt, 
              base_rect_size = 0.5/lwd_pt, 
              axis_text_angle = 0, 
              border = FALSE) +
  
  theme(
    legend.position = "none", 
    legend.key.size = unit(0.25, "cm"), 
    aspect.ratio = 1
  ) +
  
  # 添加平均值
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, 
               geom = "crossbar", width = 0.4, linewidth = 0.5/lwd_pt, color = "#000000") +
  # 添加误差线
  stat_summary(fun.data = function(x) median_hilow(x, 0.5), 
               geom = "errorbar", width = 0.2, linewidth = 0.5/lwd_pt, color = "#000000")

p_lncRNA_inter

compare_list <- list(c("UCR", "RF"))

p_lncRNA_inter_1 <- p_lncRNA_inter +
  
  stat_compare_means(
    comparisons = compare_list, 
    method = "t.test", 
    label = "p.signif", 
    bracket.size = 0.5/lwd_pt, 
    step.increase = 0.15
  )
p_lncRNA_inter_1

p_lncRNA_inter_1$layers[[which_layers(p_lncRNA_inter_1, "GeomSignif")]]$aes_params$textsize <- 7/lwd_pt   # 更改显著性标记的字体大小

plotncUCRGenesInterNum <- p_lncRNA_inter_1

# ggview(p_lncRNA_inter_1, width = 4.8, height = 6, units = "cm", dpi = 1200)
# 
# ggsave(filename = "03-results/25-Coding_UCR_genes_interaction_num/LncRNA_miRNA_inter_num.tiff", 
#        width = 4.8, height = 6, units = "cm", dpi = 1200)


# 合并来自脚本25的p2和p_lncRNA_inter_2 --------------------------------------------

pdf(file = "03-results/25-Coding_UCR_genes_interaction_num/merged_gene_inter.pdf", 
    width = 960/254, height = 600/254)
showtext_begin()

cowplot::plot_grid(plotcUCRGenesInterNum, 
                   plotncUCRGenesInterNum,
                   align = c("hv"), axis = c("tb"), 
                   nrow = 1, ncol = 2, rel_widths = c(1:1), 
                   labels = c("b", "c"), label_size = 8)

showtext_end()
dev.off()



# 20241001
pdf(file = "D:/C_英文论文/Fig.4a-d.pdf", width = 1800/254, height = 600/254)
showtext_begin()

cowplot::plot_grid(plotTissueSpecificity, plotcUCRGenesInterNum, 
                   plotncUCRGenesInterNum, plotHumanEssentialGene, 
                   align = c("hv"), axis = c("tb"), 
                   nrow = 1, ncol = 4, rel_widths = c(1:1), 
                   labels = c("a", "b", "c", "d"), label_size = 8)

showtext_end()
dev.off()




