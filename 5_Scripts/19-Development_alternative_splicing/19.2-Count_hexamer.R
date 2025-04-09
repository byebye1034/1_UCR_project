# 在human_brain_down_devAS的序列中查找所有hexamer的数目

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(Biostrings)


# count down devAS hexamers -----------------------------------------------

# get human brain down devAS fasta
human_brain_down_devAS <- readDNAStringSet(
  filepath = "02-analysis/19-Development_alternative_splicing/human_brain_down_devAS.fasta"
)

# 计算每种六连体的出现次数
hexamers <- oligonucleotideFrequency(human_brain_down_devAS, width = 6, step = 6)

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

hexamers_summary_down <- hexamers_summary

# 查看hexamers_summary
print(hexamers_summary_down)

# count up devAS hexamers -------------------------------------------------

# get human brain up devAS fasta
human_brain_up_devAS <- readDNAStringSet(
  filepath = "02-analysis/19-Development_alternative_splicing/human_brain_up_devAS.fasta"
)

# 计算每种六连体的出现次数
hexamers <- oligonucleotideFrequency(human_brain_up_devAS, width = 6, step = 6)

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

hexamers_summary_up <- hexamers_summary

# 查看hexamers_summary
print(hexamers_summary_up)

# 统计和绘图 -------------------------------------------------------------------

hexamers_summary_up$frequency <- (hexamers_summary_up$frequency)/sum(hexamers_summary_up$frequency)
hexamers_summary_down$frequency <- (hexamers_summary_down$frequency)/sum(hexamers_summary_down$frequency)

hexamers_summary_merged <- merge(hexamers_summary_up, 
                                 hexamers_summary_down, 
                                 by = "hexamers")
colnames(hexamers_summary_merged)[2] <- "frequency.up"
colnames(hexamers_summary_merged)[3] <- "frequency.down"

save(hexamers_summary_merged, 
     file = "02-analysis/19-Development_alternative_splicing/devAS_hexamers_summary_merged.Rdata")

rm(list = ls())
load("D:/R_project/UCR_project/02-analysis/19-Development_alternative_splicing/devAS_hexamers_summary_merged.Rdata")

library(ggplot2)

lwd_pt <- .pt*72.27/96

hexamers_summary_merged <- 
  hexamers_summary_merged[!(c(hexamers_summary_merged$hexamers == "AAAAAA" | hexamers_summary_merged$hexamers == "TTTTTT")), ]

# 从数据框中选择包含"GGA"的部分
SRSF1_motif <- hexamers_summary_merged[grep("GAAGGG", hexamers_summary_merged$hexamers), ]
SRSF3_motif <- hexamers_summary_merged[grep("TGGACC", hexamers_summary_merged$hexamers), ]
SRSF6_motif <- hexamers_summary_merged[grep("GAAGAA", hexamers_summary_merged$hexamers), ]
SRSF7_motif <- hexamers_summary_merged[grep("TGGACC", hexamers_summary_merged$hexamers), ]

# 绘制原始的点图
p1 <- ggplot(data = hexamers_summary_merged) +
  geom_point(mapping = aes(
    x = frequency.down, y = frequency.up
  )) +
  
  scale_x_continuous(limits = c(0, 0.00125), 
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.00125), 
                     expand = c(0, 0)) +
  
  # 添加从左下到右上的虚线
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  
  theme(
    panel.grid = element_blank(), 
    panel.background = element_blank(),
    
    axis.line.x = element_line(linewidth = 0.5/lwd_pt), 
    axis.line.y = element_line(linewidth = 0.5/lwd_pt), 
    
    aspect.ratio = 1:1
  )

p1

# 将包含"GGA"的部分绘制成红色点并添加到图中
p2 <- p1 +
  geom_point(data = SRSF1_motif, mapping = aes(
    x = frequency.down, y = frequency.up
  ), color = "red")

p2

# 名称表示 --------------------------------------------------------------------

library(dplyr)
library(ggrepel)

# 按照第二列进行排序
sorted_data <- hexamers_summary_merged %>%
  arrange(desc(frequency.down))

# 取出排名前20的数据
top_20 <- head(sorted_data, 20)
top_40 <- head(sorted_data, 40)

# 绘制原始的点图
p1 <- ggplot(data = hexamers_summary_merged) +
  geom_point(mapping = aes(
    x = frequency.down, y = frequency.up
  )) +
  
  # 添加从左下到右上的虚线
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  
  theme(
    panel.grid = element_blank(), 
    panel.background = element_blank(),
    
    axis.line.x = element_line(linewidth = 0.5/lwd_pt), 
    axis.line.y = element_line(linewidth = 0.5/lwd_pt), 
    
    aspect.ratio = 1:1
  ) +

  scale_x_continuous(limits = c(0, 0.00125), expand = c(0, 0), 
                     labels = scales::scientific) +
  scale_y_continuous(limits = c(0, 0.00125), expand = c(0, 0), 
                     labels = scales::scientific)

print(p1)

# 将排名前20的数据的第一列数据显示在点的旁边
p2 <- p1 +
  geom_text_repel(data = top_40, mapping = aes(
    x = frequency.down, y = frequency.up, label = hexamers
  ), hjust = -0.2, vjust = 0.5, size = 3, color = "darkblue")

# 打印图
print(p2)











