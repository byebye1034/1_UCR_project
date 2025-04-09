# 查看ncRNA UCR相关的lncRNA的gwRVIS
# 先看看所有的ncRNA UCR是不是都落在gwRVIS划定的范围里，如果是就不需要获得单碱基的gwRVIS了

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(ggview)
library(ggpubr)
library(ggsci)
library(ggbeeswarm)
library(Hmisc)
library(ggprism)
library(sysfonts)
library(showtext)

# gwRVIS ------------------------------------------------------------------

gwRVIS <- read.table(file = "01-data/41-JARVIS/Supplementary_File_1.tsv", 
                     header = TRUE, sep = "\t")
gwRVIS <- mutate(gwRVIS, length = end - start)
gwRVIS$chrom <- gsub("chr", "", gwRVIS$chrom)
write.table(gwRVIS[, c(1:5)], file = "02-analysis/41-JARVIS/gwRVIS.bed", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# 上传到linux，使用bedtools处理

# ncRNA UCR bed -----------------------------------------------------------

ncRNA_UCR <- read.table(file = "02-analysis/16-New_classification/02-UCR_from_Non_Coding_RNA(66).bed", 
                        sep = "\t")
ncRNA_UCR <- unique(ncRNA_UCR[, c(1:4)])
ncRNA_UCR_bed <- write.table(ncRNA_UCR, 
                             file = "02-analysis/41-JARVIS/ncRNA_UCR.bed", 
                             sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# 上传到linux，使用bedtools处理

# ncRNA UCR的gwRVIS --------------------------------------------------------

ncRNA_UCR_gw <- read.table(file = "02-analysis/41-JARVIS/ncRNA_UCR_gwRVIS.bed", 
                           sep = "\t")
ncRNA_UCR_gw <- ncRNA_UCR_gw %>% 
  mutate(start_distance = V2 - V6, 
         end_distance = V7 - V3) %>% 
  filter(start_distance > 0 & end_distance > 0)
ncRNA_UCR_gw <- ncRNA_UCR_gw[, c(4, 8, 9)]
colnames(ncRNA_UCR_gw) <- c("ncRNA_UCR_id", "gwRVIS", "JARVIS")
ncRNA_UCR_gw$type <- "ncRNA_UCR"

# ncRNA RF bed ------------------------------------------------------------

ncRNA_RF <- read.table(file = "02-analysis/16-New_classification/02-rf_from_Non_Coding_RNA.bed", 
                       sep = "\t")
ncRNA_RF <- unique(ncRNA_RF[, c(1:4)])
write.table(ncRNA_RF, file = "02-analysis/41-JARVIS/ncRNA_RF.bed", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# ncRNA RF的gwRVIS ---------------------------------------------------------

ncRNA_RF_gw <- read.table(file = "02-analysis/41-JARVIS/ncRNA_RF_gwRVIS.bed", 
                          sep = "\t")
ncRNA_RF_gw <- ncRNA_RF_gw %>% 
  mutate(start_distance = V2 - V6) %>% 
  mutate(end_distance = V7 - V3) %>% 
  filter(start_distance >= 0 & end_distance >= 0)
ncRNA_RF_gw <- ncRNA_RF_gw[, c(4, 8:9)]
colnames(ncRNA_RF_gw) <- c("ncRNA_RF_id", "gwRVIS", "JARVIS")
ncRNA_RF_gw$type <- "ncRNA_RF"

# gwRVIS ------------------------------------------------------------------

lwd_pt <- .pt*72.27/96

my_theme <- theme(
  panel.grid = element_blank(), 
  panel.background = element_blank(), 
  text = element_text(size = 5, color = "#000000"), 
  line = element_line(linewidth = 0.5/lwd_pt, color = "#000000"), 
  
  axis.title = element_text(size = 5, color = "#000000"), 
  axis.text = element_text(size = 5, color = "#000000"), 
  axis.ticks = element_line(linewidth = 0.5/lwd_pt, color = "#000000"), 
  axis.line = element_line(linewidth = 0.5/lwd_pt, color = "#000000"),
  
  legend.position = "none", 
  legend.key.size = unit(0.25, "cm"), 
  legend.title = element_text(size = 5, color = "#000000"), 
  legend.text = element_text(size = 5, color = "#000000"), 
  
  plot.title = element_text(size = 5, color = "#000000"),
  
  aspect.ratio = 1:1
)

ncRNA_UCR_RF_gwRVIS <- rbind(ncRNA_UCR_gw[, c(2, 4)], 
                             ncRNA_RF_gw[, c(2, 4)])
ncRNA_UCR_RF_gwRVIS$type <- factor(x = ncRNA_UCR_RF_gwRVIS$type, 
                                   levels = c("ncRNA_UCR", "ncRNA_RF"), 
                                   ordered = TRUE)

# 绘制基本图形
gwRVIS_p <- ggplot(data = ncRNA_UCR_RF_gwRVIS, mapping = aes(x = type, y = gwRVIS)) +
  geom_beeswarm(mapping = aes(color = type, fill = type), 
             size = 0.2, cex = 2) +
  
  scale_color_npg() +
  scale_fill_npg() +
  scale_x_discrete(labels = c("ncRNA_UCR" = "ncRNA UCR(56)", 
                              "ncRNA_RF" = "ncRNA RF(44)")) +
  scale_y_continuous(limits = c(-3, 5), expand = c(0, 0)) +
  
  labs(x = "", 
       y = "gwRVIS score") +
  
  guides(color = guide_legend(nrow = 2, ncol = 1)) +
  
  theme_prism(base_size = 7, 
              base_family = "Arial", 
              base_fontface = "plain", 
              base_line_size = 0.5/lwd_pt, 
              base_rect_size = 0.5/lwd_pt, 
              border = TRUE) +
  
  theme(
    legend.position = "top", 
    legend.key.size = unit(0.25, "cm"),
    axis.text.x = element_blank(),
    aspect.ratio = 2
  )

gwRVIS_p

ggview(gwRVIS_p, width = 4.8, height = 6, units = "cm", dpi = 1200)

# 添加中位数/均值、上、下四分位数指示线
gwRVIS_p_1 <- gwRVIS_p +
  # 添加平均值
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, 
               geom = "crossbar", width = 0.4, linewidth = 0.5/lwd_pt, color = "#000000") +
  # 添加误差线
  stat_summary(fun.data = function(x) median_hilow(x, 0.5), 
               geom = "errorbar", width = 0.2, linewidth = 0.5/lwd_pt, color = "#000000")

gwRVIS_p_1

ggview(gwRVIS_p_1, width = 4.8, height = 6, units = "cm", dpi = 1200)

compare_list <- list(c("ncRNA_UCR", "ncRNA_RF"))

gwRVIS_p_2 <- gwRVIS_p_1 +
  
  stat_compare_means(
    comparisons = compare_list, 
    method = "t.test", 
    label = "p.signif", 
    bracket.size = 0.5/lwd_pt, 
    step.increase = 0.15
  )
gwRVIS_p_2

ggview(gwRVIS_p_2, width = 4.8, height = 6, units = "cm", dpi = 1200)

ggsave(filename = "03-results/41-JARVIS/gwRVIS.tiff", 
       width = 4.8, height = 6, units = "cm", dpi = 1200)

font_add(family = "Arial", regular = "C:/Windows/Fonts/arial.ttf")
font_families()

pdf(file = "03-results/41-JARVIS/gwRVIS.pdf", width = 480/254, height = 600/254)
showtext_begin()

gwRVIS_p_2

showtext_end()
dev.off()




