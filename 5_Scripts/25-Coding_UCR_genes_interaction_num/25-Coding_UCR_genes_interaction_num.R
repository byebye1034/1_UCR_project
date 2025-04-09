# 使用cell reports的molecular interactions网络比较UCR和RF的interaction num是否有显著差异

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library(readxl)
library(biomaRt)
library(curl)
library(ggprism)

coding_UCR_genes <- 
  read.table(file = "02-analysis/16-New_classification/01-UCR_from_protein_coding_gene.bed", 
             sep = "\t")

coding_RF_genes <- 
  read.table(file = "02-analysis/16-New_classification/01-rf_form_protein_coding_gene.bed", 
             sep = "\t")

# 利用useMart链接到人类的数据库，其他数据库使用listDatasets(ensembl)查看
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
listAttributes(ensembl)

# attr为基因id目标类型，listAttributes(ensembl)展示其他id输出类型
attr <- c("ensembl_gene_id", "uniprotswissprot")

# listFilters(ensembl)展示其他id输入类型
ids <- getBM(attributes = attr, 
             filters = "ensembl_gene_id", 
             values = coding_UCR_genes$V8, 
             mart = ensembl)

ids[(ids$ensembl_gene_id == "ENSG00000248643"), ][, 2] <- "Q96PK6"
ids[(ids$ensembl_gene_id == "ENSG00000280987"), ][, 2] <- "P43243"

ids <- ids[!(ids$uniprotswissprot == ""), ]
coding_UCR_entry_id <- ids

# 获得coding RF genes的entry id
coding_RF_entry_id <- getBM(attributes = attr, 
                            filters = "ensembl_gene_id", 
                            values = coding_RF_genes$V8, 
                            mart = ensembl)
coding_RF_entry_id[(coding_RF_entry_id$ensembl_gene_id == "ENSG00000263715"), ][, 2] <- "P34998"
coding_RF_entry_id[(coding_RF_entry_id$ensembl_gene_id == "ENSG00000282301"), ][, 2] <- "P24462"

coding_RF_entry_id <- coding_RF_entry_id[!(coding_RF_entry_id$uniprotswissprot == ""), ]

# 读取cell reports里的相互作用网络
global_PPI_network <- 
  read.table(file = "01-data/13-Network_analysis/Global Protein-Protein Interaction Network.txt", 
             sep = "\t", header = TRUE)

global_PPI_network <- mutate(global_PPI_network, 
                             internum = rowSums(global_PPI_network[, c(3:75)]))
global_PPI_network <- filter(global_PPI_network, internum > 0)
save(global_PPI_network, file = "02-analysis/25-Coding_UCR_genes_interaction_num/global_PPI_network.Rdata")

# 能够发生的所有的互作都在这里
global_PPI_network_entry_id <- global_PPI_network[, c(1, 2)]
global_PPI_network_entry_id <- unique(global_PPI_network_entry_id)

# coding UCR genes --------------------------------------------------------

# 创建一个空的向量用于存储每个值出现的次数
counts <- numeric(length(coding_UCR_entry_id$uniprotswissprot))

# 循环遍历每个值并统计出现次数
for (i in seq_along(coding_UCR_entry_id$uniprotswissprot)) {
  counts[i] <- sum(global_PPI_network_entry_id$InteractorA == coding_UCR_entry_id$uniprotswissprot[i] | 
                     global_PPI_network_entry_id$InteractorB == coding_UCR_entry_id$uniprotswissprot[i])
}

# 将统计结果添加到数据框中
coding_UCR_entry_id$Count <- counts

# coding RF genes ---------------------------------------------------------

# 创建一个空的向量用于存储每个值出现的次数
counts <- numeric(length(coding_RF_entry_id$uniprotswissprot))

# 循环遍历每个值并统计出现次数
for (i in seq_along(coding_RF_entry_id$uniprotswissprot)) {
  counts[i] <- sum(global_PPI_network_entry_id$InteractorA == coding_RF_entry_id$uniprotswissprot[i] | 
                     global_PPI_network_entry_id$InteractorB == coding_RF_entry_id$uniprotswissprot[i])
}

# 将统计结果添加到数据框中
coding_RF_entry_id$Count <- counts

save(coding_UCR_entry_id, coding_RF_entry_id, 
     file = "02-analysis/25-Coding_UCR_genes_interaction_num/coding_genes_entry_id_internum.Rdata")

# 绘制统计图 -------------------------------------------------------------------

rm(list = ls())

library(ggplot2)
library(ggsci)
library(ggbeeswarm)
library(Hmisc)
library(ggpubr)
library(ggprism)
library(sysfonts)
library(showtext)
library(ggbreak)
library(scales)

load("D:/R_project/UCR_project/02-analysis/25-Coding_UCR_genes_interaction_num/coding_genes_entry_id_internum.Rdata")

result <- wilcox.test(coding_UCR_entry_id$Count, coding_RF_entry_id$Count)
print(result)
# Wilcoxon rank sum test with continuity correction
# 
# data:  coding_UCR_entry_id$Count and coding_RF_entry_id$Count
# W = 28421, p-value = 0.0009749
# alternative hypothesis: true location shift is not equal to 0

count <- c(coding_UCR_entry_id$Count, coding_RF_entry_id$Count)
interaction_num <- data.frame(type = c(rep("UCR", 202), rep("RF", 238)), 
                              num = count)
interaction_num$type <- factor(x = interaction_num$type, 
                               levels = c("UCR", "RF"))

lwd_pt <- .pt*72.27/96
font_add(family = "Arial", regular = "C:/Windows/Fonts/arial.ttf")
font_families()

BeeswarmPlot <- function(data){
  ggplot(data = data, mapping = aes(x = type, y = log10(num + 1))) +
    geom_beeswarm(mapping = aes(fill = type, color = type), 
                  size = 0.2, cex = 1.2) +
    
    scale_y_continuous(limits = c(0, 4), 
                       expand = c(0, 0)) +
    scale_color_npg() +
    
    scale_x_discrete(labels = c("UCR" = "UCR (n = 202)", "RF" = "Randf (n = 238)")) +
    
    labs(x = "", 
         y = "log10(No. of PPI)", 
         title = "") +
    
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
}

BeeswarmPlot(data = interaction_num)
BeeswarmPlot(data = interaction_num) + canvas(width = 4.8, height = 6, units = "cm", dpi = 1200)

compare_list <- list(c("UCR", "RF"))

p2 <- BeeswarmPlot(data = interaction_num) +
  stat_compare_means(
    comparisons = compare_list, 
    method = "wilcox.test", 
    label = "p.signif", 
    bracket.size = 0.5/lwd_pt, 
    step.increase = 0.25, 
    label.x = 590, 
    label.y = 3.5
  )
p2

pdf(file = "03-results/25-Coding_UCR_genes_interaction_num/cUCRGenesInterNum.pdf", 
    width = 480/254, height = 600/254)
showtext_begin()

p2

showtext_end()
dev.off()

plotcUCRGenesInterNum <- p2

# # 20250218 pie plot -------------------------------------------------------
# 
# # 修改后的分组函数
# create_groups <- function(df) {
#   # 分离num=0的组
#   zero_group <- df %>% filter(Count == 0)
#   non_zero <- df %>% filter(Count > 0)
#   
#   # 计算非零组的绝对数值范围
#   min_val <- min(non_zero$Count)
#   max_val <- max(non_zero$Count)
#   
#   # 生成基于绝对值的分界点（分成4个等宽区间）
#   breaks <- seq(min_val, max_val, length.out = 5)  # 生成5个点形成4个区间
#   
#   # 创建动态分组标签
#   non_zero <- non_zero %>%
#     mutate(
#       group = cut(Count,
#                   breaks = breaks,
#                   labels = sprintf("%.1f-%.1f", breaks[-length(breaks)], breaks[-1]),
#                   include.lowest = TRUE)
#     )
#   
#   # 合并所有组并计算比例
#   bind_rows(
#     zero_group %>% mutate(group = "0%"),
#     non_zero
#   ) %>%
#     count(group) %>%
#     mutate(percentage = n / sum(n) * 100)
# }
# 
# # 处理测试数据
# df1 <- create_groups(coding_UCR_entry_id)
# df2 <- create_groups(coding_RF_entry_id)
# 
# # 查看分组结果
# print(test_result)
# 
# # 修改后的绘图函数（保持相同）
# plot_pie <- function(data, title) {
#   ggplot(data, aes(x = "", y = percentage, fill = group)) +
#     geom_bar(width = 1, stat = "identity", color = "white") +
#     coord_polar("y", start = 0) +
#     geom_text(aes(label = paste0(round(percentage, 1), "%"))) +
#                 labs(title = title, fill = "Interaction Group") +
#                 theme_void() +
#                 scale_fill_brewer(palette = "Set3") +
#                 theme(plot.title = element_text(hjust = 0.5, face = "bold"))
# }
# 
# # 绘制测试饼图
# plot_pie(df1, "No. of protein-protein interaction UCRs")
# plot_pie(df2, "No. of protein-protein interaction RF")


# CDF 累积分布函数图 -------------------------------------------------------------

coding_UCR_entry_id <- mutate(coding_UCR_entry_id, Group = "UCR")
coding_RF_entry_id <- mutate(coding_RF_entry_id, Group = "RF")
data <- rbind(coding_UCR_entry_id, coding_RF_entry_id)

CDFPlot <- function(data){
  ggplot(data = data, mapping = aes(x = Count, color = Group)) +
    stat_ecdf(geom = "step", linewidth = 0.5/lwd_pt, alpha = 0.5) +
    scale_y_continuous(
      labels = scales::percent_format(), 
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      title = "Cumulative Distribution of Interaction Counts",
      x = "Number of Interacting Proteins (log scale)",
      y = "Cumulative Proportion",
      color = "Gene Group"
    ) +
    scale_color_manual(values = c("#3B9AB2", "#E1AF00")) +
    theme_minimal(base_size = 7) +
    theme(
      legend.position = c(0.85, 0.25),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    geom_vline(xintercept = c(10, 100), 
               linetype = "dashed", 
               color = "grey50",
               alpha = 0.6)
}

CDFPlot(data = data) + canvas(width = 4.8, height = 6, units = "cm", dpi = 1200)








