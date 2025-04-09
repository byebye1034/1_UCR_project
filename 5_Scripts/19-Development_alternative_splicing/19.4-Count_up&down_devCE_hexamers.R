# human brain up&down devCE hexamers summary

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(Biostrings)
library(ggplot2)
library(ggrepel)
library(ggview)

# up ----------------------------------------------------------------------

# up
human_brain_up_devCE <- readDNAStringSet(
  filepath = "02-analysis/19-Development_alternative_splicing/human_brain_up_devCE.fasta"
)

# 计算每种六连体的出现次数
hexamers <- oligonucleotideFrequency(human_brain_up_devCE, width = 6, step = 6)

hexamers_df <- as.data.frame(hexamers)

# 创建一个空的数据框hexamers_summary
hexamers_summary <- data.frame(hexamers = character(), frequency = numeric())

# 遍历hexamers_df的列名，并计算每列的和
for (col_name in colnames(hexamers_df)) {
  # 计算每列的和
  col_sum <- sum(hexamers_df[[col_name]])
  
  # 创建一个新的行，并添加到hexamers_summary中
  new_row <- data.frame(hexamers = col_name, frequency = col_sum)
  hexamers_summary <- rbind(hexamers_summary, new_row)
}

hexamers_summary_up_devCE <- hexamers_summary

# 查看hexamers_summary
print(hexamers_summary_up_devCE)

# down --------------------------------------------------------------------

# down
human_brain_down_devCE <- readDNAStringSet(
  filepath = "02-analysis/19-Development_alternative_splicing/human_brain_down_devCE.fasta"
)

# 计算每种六连体的出现次数
hexamers <- oligonucleotideFrequency(human_brain_down_devCE, width = 6, step = 6)

hexamers_df <- as.data.frame(hexamers)

# 创建一个空的数据框hexamers_summary
hexamers_summary <- data.frame(hexamers = character(), frequency = numeric())

# 遍历hexamers_df的列名，并计算每列的和
for (col_name in colnames(hexamers_df)) {
  # 计算每列的和
  col_sum <- sum(hexamers_df[[col_name]])
  
  # 创建一个新的行，并添加到hexamers_summary中
  new_row <- data.frame(hexamers = col_name, frequency = col_sum)
  hexamers_summary <- rbind(hexamers_summary, new_row)
}

hexamers_summary_down_devCE <- hexamers_summary

# 查看hexamers_summary
print(hexamers_summary_down_devCE)

# 统计和绘图 -------------------------------------------------------------------

hexamers_summary_up_devCE$frequency <- hexamers_summary_up_devCE$frequency/sum(hexamers_summary_up_devCE$frequency)
hexamers_summary_down_devCE$frequency <- hexamers_summary_down_devCE$frequency/sum(hexamers_summary_down_devCE$frequency)

hexamers_summary_merged <- merge(hexamers_summary_up_devCE, 
                                 hexamers_summary_down_devCE, 
                                 by = "hexamers")
colnames(hexamers_summary_merged)[2] <- "frequency.up"
colnames(hexamers_summary_merged)[3] <- "frequency.down"

save(hexamers_summary_merged, 
     file = "02-analysis/19-Development_alternative_splicing/devCE_hexamers_summary_merged.Rdata")

# plot --------------------------------------------------------------------

rm(list = ls())
load("D:/R_project/UCR_project/02-analysis/19-Development_alternative_splicing/devCE_hexamers_summary_merged.Rdata")

lwd_pt <- .pt*72.27/96

hexamers_summary_merged_notop <- anti_join(hexamers_summary_merged, down_top_10)
hexamers_summary_merged_notop <- anti_join(hexamers_summary_merged_notop, up_top_10)

#   hexamers frequency.up frequency.down
# 1   AAAAAA  0.001567612    0.002387757
# 2   TTTTTT  0.001385860    0.001496987

# 绘制原始的点图
p1 <- ggplot(data = hexamers_summary_merged_notop) +
  geom_point(mapping = aes(
    x = frequency.down, y = frequency.up
  ), 
  size = 0.35, position = "jitter") +
  
  # 添加从左下到右上的虚线
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  
  labs(
    x = "human brain down devCE", 
    y = "human brain up devCE", 
    title = ""
  ) +
  
  theme(
    panel.grid = element_blank(), 
    panel.background = element_blank(),
    panel.border = element_rect(linewidth = 0.5/lwd_pt, 
                                color = "black", 
                                fill = NA), 
    
    axis.text = element_text(size = 5, color = "black"), 
    axis.line.x = element_line(linewidth = 0.5/lwd_pt), 
    axis.line.y = element_line(linewidth = 0.5/lwd_pt), 
    axis.ticks = element_line(linewidth = 0.5/lwd_pt), 
    axis.title = element_text(size = 5, color = "black"), 
    
    aspect.ratio = 1:1
  ) +
  
  scale_x_continuous(limits = c(0, 0.00150), 
                     expand = c(0.000125, 0), 
                     labels = scales::scientific) +
  scale_y_continuous(limits = c(0, 0.00150), 
                     expand = c(0.000125, 0), 
                     labels = scales::scientific)
p1

# 按照down进行排序
sorted_data <- hexamers_summary_merged %>%
  arrange(desc(frequency.down))

# 取出down排名前10的数据
down_top_10 <- head(sorted_data, 10)

# 按照up进行排序
sorted_data <- hexamers_summary_merged %>%
  arrange(desc(frequency.up))

# 取出up排名前10的数据
up_top_10 <- head(sorted_data, 10)

# 将排名前10的数据点分别展示为不同颜色
p2 <- p1 +
  geom_point(data = up_top_10, 
             mapping = aes(
               x = frequency.down, 
               y = frequency.up
             ), 
             color = "#4DBBD5FF", 
             size = 0.35) +
  
  geom_text_repel(data = up_top_10, 
                  mapping = aes(
                    x = frequency.down, 
                    y = frequency.up, 
                    label = hexamers
                  ), 
                  hjust = -0.2, vjust = 0.5, 
                  size = 1, color = "#4DBBD5FF", 
                  segment.size = 0.5/lwd_pt, 
                  max.overlaps = 20) +
  
  geom_point(data = down_top_10, 
             mapping = aes(
               x = frequency.down, 
               y = frequency.up
             ), 
             color = "#E64B35FF", 
             size = 0.35) +
  
  geom_text_repel(data = down_top_10, 
                  mapping = aes(
                    x = frequency.down, 
                    y = frequency.up, 
                    label = hexamers
                  ), 
                  hjust = 1, vjust = 0.5, 
                  size = 1, color = "#E64B35FF", 
                  segment.size = 0.5/lwd_pt, 
                  max.overlaps = 20)

p2

pdf(file = "03-results/19-Development_alternative_splicing/GGA_enrichment.pdf", 
    width = 960/254, height = 720/254)
p2
dev.off()

cowplot::plot_grid(p3, p4, p2, p3, 
                   align = c("hv"), axis = c("tb"), 
                   nrow = 1, ncol = 4, rel_widths = c(1,1,1,1), 
                   labels = c("A", "B", "C", "D"), label_size = 8)

pdf(file = "D:/C_英文论文/fig3组合图.pdf", width = 1650/254, height = 500/254)

cowplot::plot_grid(p3, p4, p2, p3, 
                   align = c("hv"), axis = c("tb"), 
                   nrow = 1, ncol = 4, rel_widths = c(1,1,1,1), 
                   labels = c("A", "B", "C", "D"), label_size = 8)

dev.off()

# HNRNPM ------------------------------------------------------------------

HNRNPM_motif <- hexamers_summary_merged[grep("GTTTGT", hexamers_summary_merged$hexamers), ]

p3 <- p1 +
  geom_point(data = HNRNPM_motif, 
             mapping = aes(
               x = frequency.down, 
               y = frequency.up
             ), 
             color = "#4DBBD5FF") +
  
  geom_text_repel(data = HNRNPM_motif, 
                  mapping = aes(
                    x = frequency.down, 
                    y = frequency.up, 
                    label = hexamers
                  ), 
                  hjust = -0.2, vjust = 0.5, 
                  size = 3, color = "#4DBBD5FF")
p3

# enrich in up

# HNRNPK ------------------------------------------------------------------

HNRNPK_motif <- hexamers_summary_merged[grep("CCCC.T", hexamers_summary_merged$hexamers), ]

p4 <- p1 +
  geom_point(data = HNRNPK_motif, 
             mapping = aes(
               x = frequency.down, 
               y = frequency.up
             ), 
             color = "#4DBBD5FF") +
  
  geom_text_repel(data = HNRNPK_motif, 
                  mapping = aes(
                    x = frequency.down, 
                    y = frequency.up, 
                    label = hexamers
                  ), 
                  hjust = -0.2, vjust = 0.5, 
                  size = 3, color = "#4DBBD5FF")
p4

# enrich in

# PCBP2 -------------------------------------------------------------------

PCBP2_motif <- hexamers_summary_merged[grep("CTTTCC", hexamers_summary_merged$hexamers), ]

p5 <- p1 +
  geom_point(data = PCBP2_motif, 
             mapping = aes(
               x = frequency.down, 
               y = frequency.up
             ), 
             color = "#4DBBD5FF") +
  
  geom_text_repel(data = PCBP2_motif, 
                  mapping = aes(
                    x = frequency.down, 
                    y = frequency.up, 
                    label = hexamers
                  ), 
                  hjust = -0.2, vjust = 0.5, 
                  size = 3, color = "#4DBBD5FF")
p5

# enrich in down(CTTTCC)

# QKI ---------------------------------------------------------------------

QKI_motif <- hexamers_summary_merged[grep("TACTAA", hexamers_summary_merged$hexamers), ]

p6 <- p1 +
  geom_point(data = QKI_motif, 
             mapping = aes(
               x = frequency.down, 
               y = frequency.up
             ), 
             color = "#4DBBD5FF") +
  
  geom_text_repel(data = QKI_motif, 
                  mapping = aes(
                    x = frequency.down, 
                    y = frequency.up, 
                    label = hexamers
                  ), 
                  hjust = -0.2, vjust = 0.5, 
                  size = 3, color = "#4DBBD5FF")
p6

# enrich in both

# GGA ---------------------------------------------------------------------

SRSF1_motif <- hexamers_summary_merged[grep("GGA", hexamers_summary_merged$hexamers), ]

p7 <- p1 +
  geom_point(data = SRSF1_motif, 
             mapping = aes(
               x = frequency.down, 
               y = frequency.up
             ), 
             color = "#4DBBD5FF") +
  
  geom_text_repel(data = SRSF1_motif, 
                  mapping = aes(
                    x = frequency.down, 
                    y = frequency.up, 
                    label = hexamers
                  ), 
                  hjust = -0.2, vjust = 0.5, 
                  size = 3, color = "#4DBBD5FF")
p7

# enrich in both



























