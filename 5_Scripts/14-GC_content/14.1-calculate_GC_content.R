#

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(dbplyr)
library(Biostrings)

# 整理序列 --------------------------------------------------------------------

## UCR left
UCR_left <- readDNAStringSet(filepath = "02-analysis/14-GC_content/ucr_left_homemade_GRCh38p14.fa")
print(UCR_left)
names(UCR_left) <- sub("::.*", "_left", names(UCR_left))

# 保存到一个新的FASTA文件
writeXStringSet(UCR_left, "02-analysis/14-GC_content/UCR_left.fasta")

## UCR right
UCR_right <- readDNAStringSet(filepath = "02-analysis/14-GC_content/ucr_right_homemade_GRCh38p14.fa")
print(UCR_right)
names(UCR_right) <- sub("::.*", "_right", names(UCR_right))
print(UCR_right)

# 保存到一个新的FASTA文件
writeXStringSet(UCR_right, "02-analysis/14-GC_content/UCR_right.fasta")

## rf
rf <- readDNAStringSet(filepath = "02-analysis/14-GC_content/rf_homemade_GRCh38p14.fa")
print(rf)
names(rf) <- sub("::.*", "", names(rf))
print(rf)

# 保存到一个新的FASTA文件
writeXStringSet(rf, "02-analysis/14-GC_content/rf.fasta")

# 计算GC含量 ------------------------------------------------------------------

UCR_left_GC <- sapply(UCR_left, function(seq) {
  sum(letterFrequency(seq, letters = c("G", "C"))) / sum(letterFrequency(seq, letters = c("G", "C", "A", "T", "N"))) * 100
})
print(UCR_left_GC)

UCR_right_GC <- sapply(UCR_right, function(seq){
  sum(letterFrequency(seq, letters = c("G", "C"))) / sum(letterFrequency(seq, letters = c("G", "C", "A", "T", "N"))) * 100
})

rf_GC <- sapply(rf, function(seq){
  sum(letterFrequency(seq, letters = c("G", "C"))) / sum(letterFrequency(seq, letters = c("G", "C", "A", "T", "N"))) * 100
})

UCR <- readDNAStringSet(file = "01-data/UCR_raw/ucr_sequences_simple_name.fasta")
UCR_GC <- sapply(UCR, function(seq){
  sum(letterFrequency(seq, letters = c("G", "C"))) / sum(letterFrequency(seq, letters = c("G", "C", "A", "T", "N"))) * 100
})

GC_content <- data.frame(UCR = UCR_GC, UCR_left = UCR_left_GC, 
                         UCR_right = UCR_right_GC, rf = rf_GC)
save(GC_content, file = "02-analysis/14-GC_content/GC_content_sum.Rdata")

# GC含量绘图 ------------------------------------------------------------------

rm(list = ls())

library(reshape2)
library(dbplyr)
library(ggplot2)
library(ggsci) # 配色
library(ggpubr) # 添加显著性标记
library(gghalves) # 云雨图
library(gginnards) # 更改显著性标记大小
library(ggprism)
library(sysfonts)
library(showtext)

font_add(family = "Arial", regular = "C:/Windows/Fonts/arial.ttf")
font_families()

load("D:/R_project/UCR_project/02-analysis/14-GC_content/GC_content_sum.Rdata")

GC_content_L <- melt(data = GC_content, 
                     id.vars = NULL, 
                     variable.name = "Group", 
                     value.name = "GC_content")

lwd_pt <- .pt*72.27/96
compare_list <- list(c("UCR", "UCR_left"), 
                     c("UCR", "UCR_right"), 
                     c("UCR", "rf"))

pdf(file = "03-results/14-GC_content/GC_content_compare.pdf",
    width = 450/254, height = 600/254)
showtext_begin()

p1 <- ggplot(data = GC_content_L, aes(x = Group, y = GC_content, fill = Group, group = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.15, alpha = 0.4, linewidth = 0.5/lwd_pt) +
  
  stat_compare_means(
    method = "wilcox.test", 
    comparisons = compare_list, 
    aes(group = Group),  # 指定分组变量
    label = "p.format",
    method.args = list(alternative = "two.sided")  # 设置检验方向
  ) +
  
  scale_y_continuous(limits = c(0, 120), expand = c(0, 0)) +
  
  scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF"), 
                    name = "", 
                    breaks = c("UCR", "UCR_left", "UCR_right", "rf")) +
  
  guides(fill = guide_legend(nrow = 2, ncol = 2)) +
  
  labs(x = "", y = "GC content") +
  
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

p1$layers[[which_layers(p1, "GeomSignif")]]$aes_params$textsize <- 5/lwd_pt
p1

GC_raincloud_plot <- p1 +
  geom_half_violin(mapping = aes(x = Group, y = GC_content, fill = Group), 
                   linewidth = 0.5/lwd_pt, 
                   position = position_nudge(x = 0.1, y = 0),   # 位置调整，这里将其向右水平移动0.1
                   side = "R",   # 显示哪一侧， "I"代表左侧，"R"代表右侧，默认"I"
                   adjust = 1.2,   # 调整带宽，这里设为1.2使宽带略变平滑
                   trim = F,   # 小提琴图尾部的数据修整，默认为"T",表示将尾部修整到数据范围；"F"表示不修剪尾部
                   alpha = 0.4)
GC_raincloud_plot

GC_raincloud_plot_3 <- GC_raincloud_plot +
  geom_half_point(mapping = aes(x = as.numeric(Group) - 0.1, y = GC_content, color = Group),   # 散点位置向左平移0.2
                  position = position_jitter(width = 0.02),   #调整散点，使取值相同的原重合散点分散开
                  size = 0.05, 
                  side = "I") +
  
  scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF")) +
  guides(color = FALSE)
GC_raincloud_plot_3

showtext_end()
dev.off()





