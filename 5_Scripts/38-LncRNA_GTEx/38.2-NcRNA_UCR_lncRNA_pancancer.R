# 使用easyTCGA对ncRNA相关的lncRNA进行泛癌分析

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(easyTCGA)
library(stringr)
library(ggplot2)
library(forcats)
library(ggview)

# 安装easyTCGA --------------------------------------------------------------

BiocManager::install("TCGAbiolinks")
library(devtools)
devtools::install_github("ayueme/easyTCGA")
library(easyTCGA)

# 使用TCGA ------------------------------------------------------------------

ncRNA_UCR_lncRNA <- read.table(file = "02-analysis/16-New_classification/ncRNA_UCR(66)_ensembl_id.txt")
ncRNA_UCR_lncRNA <- unique(ncRNA_UCR_lncRNA)
colnames(ncRNA_UCR_lncRNA) <- "Gene.ID"

ACC <- read.table(file = "01-data/38-LncRNA_GTEx/GEPIA2_pancaner/ACC.txt", 
                  sep = "\t", header = TRUE)
ACC$Gene.ID <- str_sub(ACC$Gene.ID, 1, 15)
ACC.res <- merge(ncRNA_UCR_lncRNA, ACC, by = "Gene.ID", all.x = TRUE)
ACC.res$Cancer.Type <- "ACC"

Allcancer.res <- ACC.res

file_list <- list.files(path = "01-data/38-LncRNA_GTEx/GEPIA2_pancaner/")
file_list <- paste("01-data/38-LncRNA_GTEx/GEPIA2_pancaner/", file_list, sep = "")

for (file in file_list) {
  # 从文件名中去掉后缀名，并作为数据框的名称
  df_name <- tools::file_path_sans_ext(file)
  df_name <- str_sub(df_name, 40, nchar(df_name))
  df_name_id <- df_name
  
  df_name <- read.table(file = file, sep = "\t", header = TRUE)
  df_name$Gene.ID <- str_sub(df_name$Gene.ID, 1, 15)
  df_name.res <- merge(ncRNA_UCR_lncRNA, df_name, by = "Gene.ID", all.x = TRUE)
  df_name.res$Cancer.Type <- df_name_id
  
  Allcancer.res <- rbind(Allcancer.res, df_name.res)
}

# 如果一个癌症当中全部都是缺失值就删除
# 假设df是你的数据框，按第七列进行分组，检查每组第四列是否全部为缺失值
Allcancer.res <- Allcancer.res %>% 
  group_by(Cancer.Type) %>% 
  filter(sum(is.na(Log2.Fold.Change.)) != 57) %>% 
  ungroup()
Allcancer.res <- filter(Allcancer.res, Cancer.Type != "ACC")

Allcancer.res$Cancer.Type <- factor(x = Allcancer.res$Cancer.Type, 
                                    levels = unique(Allcancer.res$Cancer.Type), 
                                    ordered = TRUE)
Allcancer.res$Cancer.Type <- fct_rev(Allcancer.res$Cancer.Type)

save(Allcancer.res, file = "02-analysis/38-LncRNA_GTEx/Allcancer.res.Rdata")
load("D:/R_project/UCR_project/02-analysis/38-LncRNA_GTEx/Allcancer.res.Rdata")

Allcancer.res <- Allcancer.res %>% 
  group_by(Gene.ID) %>% 
  filter(sum(is.na(Log2.Fold.Change.)) != 19) %>% 
  ungroup()

lwd_pt <- .pt*72.27/96

p_pancer <- ggplot(data = Allcancer.res) +
  geom_point(mapping = aes(x = Cancer.Type, y = Gene.ID, 
                           color = Log2.Fold.Change., size = -log10(adjp))) +
  
  scale_color_gradient2(low = "#4DBBD5FF", high = "#E64B35FF", mid = "white") +
  
  theme(
    panel.grid.major = element_line(color = "grey60"), 
    panel.background = element_blank(), 
    panel.border = element_rect(linewidth = 0.5/lwd_pt, color = "#000000", fill = NA), 
    text = element_text(size = 7, color = "#000000"), 
    line = element_line(linewidth = 0.5/lwd_pt, color = "#000000"), 
    
    axis.title = element_blank(), 
    axis.text = element_text(size = 7, color = "#000000"), 
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.ticks = element_line(linewidth = 0.5/lwd_pt, color = "#000000"), 
    axis.line = element_line(linewidth = 0.5/lwd_pt, color = "#000000"),
    
    legend.position = "top", 
    legend.key.size = unit(0.25, "cm"), 
    legend.title = element_text(size = 7, color = "#000000"), 
    legend.text = element_text(size = 7, color = "#000000"), 
    
    plot.title = element_text(size = 7, color = "#000000")
  ) +
  
  coord_flip()

p_pancer

ggview(p_pancer, width = 16, height = 12, units = "cm", dpi = 1200)

ggsave(filename = "03-results/38-LncRNA_GTEx/lncRNA_pancncer.tiff", 
       width = 16, height = 12, units = "cm", dpi = 1200)













