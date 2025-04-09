# 统计人类七种器官中不同类别的devAS的数目的数目和占比

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(dplyr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(cowplot)
library(ggprism)
library(sysfonts)
library(showtext)
library(ggview)

human_devAS <- read.table(file = "01-data/19-Development_alternative_splicing/human.devAS", 
                          sep = ",", header = TRUE)
print(table(human_devAS$as.type))

human_devAS <- human_devAS[, c(1:15)]

# get organ devAS
human_brain_devAS <- human_devAS[c(human_devAS$pattern.brain %in% c("u", "d", "ud", "du")), ]
# > nrow(human_brain_devAS)
# [1] 9624
# > table(human_brain_devAS$as.type)
# 
#   AA      AD      CE complex      RI 
# 1445    1451    3693      91    2944
human_brain_devAS <- human_brain_devAS[, c(2, 9)]

human_cerebellum_devAS <- human_devAS[c(human_devAS$pattern.cerebellum %in% c("u", "d", "ud", "du")), ]
# > nrow(human_cerebellum_devAS)
# [1] 8889
# > table(human_cerebellum_devAS$as.type)
# 
#   AA      AD      CE complex      RI 
# 1274    1233    3598      82    2702
human_cerebellum_devAS <- human_cerebellum_devAS[, c(2, 10)]

human_heart_devAS <- human_devAS[c(human_devAS$pattern.heart %in% c("u", "d", "ud", "du")), ]
# > nrow(human_heart_devAS)
# [1] 6430
# > table(human_heart_devAS$as.type)
# 
#  AA      AD      CE complex      RI 
# 885     998    1620      58    2869
human_heart_devAS <- human_heart_devAS[, c(2, 11)]

human_kidney_devAS <- human_devAS[c(human_devAS$pattern.kidney %in% c("u", "d", "ud", "du")), ]
# > nrow(human_kidney_devAS)
# [1] 5147
# > table(human_kidney_devAS$as.type)
# 
# AA      AD      CE complex      RI 
# 702     752    1482      41    2170
human_kidney_devAS <- human_kidney_devAS[, c(2, 12)]

human_liver_devAS <- human_devAS[c(human_devAS$pattern.liver %in% c("u", "d", "ud", "du")), ]
# > nrow(human_liver_devAS)
# [1] 1978
# > table(human_liver_devAS$as.type)
# 
# AA      AD      CE complex      RI 
# 272     303     726      30     647 
human_liver_devAS <- human_liver_devAS[, c(2, 13)]

human_ovary_devAS <- human_devAS[c(human_devAS$pattern.ovary %in% c("u", "d", "ud", "du")), ]
# > nrow(human_ovary_devAS)
# [1] 1811
# > table(human_ovary_devAS$as.type)
# 
# AA      AD      CE complex      RI 
# 257     237     781      17     519 
human_ovary_devAS <- human_ovary_devAS[, c(2, 14)]

human_testis_devAS <- human_devAS[c(human_devAS$pattern.testis %in% c("u", "d", "ud", "du")), ]
# > nrow(human_testis_devAS)
# [1] 11588
# > table(human_testis_devAS$as.type)
# 
# AA      AD      CE complex      RI 
# 1702    1774    3791     100    4221
human_testis_devAS <- human_testis_devAS[, c(2, 15)]

devAS_num <- human_brain_devAS %>% 
  merge(human_cerebellum_devAS, by = "seg.id", all = TRUE) %>% 
  merge(human_heart_devAS, by = "seg.id", all = TRUE) %>% 
  merge(human_liver_devAS, by = "seg.id", all = TRUE) %>% 
  merge(human_kidney_devAS, by = "seg.id", all = TRUE) %>% 
  merge(human_ovary_devAS, by = "seg.id", all = TRUE) %>% 
  merge(human_testis_devAS, by = "seg.id", all = TRUE)

devAS_num <- merge(devAS_num, 
                   human_devAS[, c(2, 8)], 
                   by = "seg.id", 
                   all.x = TRUE)

save(devAS_num, 
     file = "02-analysis/19-Development_alternative_splicing/devAS_num.Rdata")


# 不同组织中devAS的数目占比 ---------------------------------------------------------

load("D:/R_project/UCR_project/02-analysis/19-Development_alternative_splicing/devAS_num.Rdata")

devAS_num_long <- melt(
  devAS_num, 
  id.vars = c("seg.id", "as.type"), 
  measure.vars = c(
    "pattern.brain", 
    "pattern.cerebellum", 
    "pattern.heart", 
    "pattern.liver", 
    "pattern.kidney", 
    "pattern.ovary", 
    "pattern.testis"
  ), 
  variable.name = "organ", 
  value.name = "pattern"
)

# 绘制不同器官当中不同类型可变剪接事件占比：堆叠柱形图
devAS_num_long_1 <- devAS_num_long[!(is.na(devAS_num_long$pattern)), ]

lwd_pt <- .pt*72.27/96

# 按照devAS数目进行排序
devAS_num_long_1$organ <- factor(x = devAS_num_long_1$organ, 
                                 levels = c("pattern.testis", 
                                            "pattern.brain", 
                                            "pattern.cerebellum", 
                                            "pattern.heart", 
                                            "pattern.kidney",
                                            "pattern.liver", 
                                            "pattern.ovary"))

font_add(family = "Arial", regular = "C:/Windows/Fonts/arial.ttf")
font_families()

p1 <- ggplot(data = devAS_num_long_1) +
  
  geom_bar(mapping = aes(x = organ, fill = organ), width = 0.5) +
  
  scale_y_continuous(limits = c(0, 12000), expand = c(0, 0)) +
  scale_x_discrete(labels = c("testis", "brain", 
                              "cerebellum", "heart", "kidney", 
                              "liver", "ovary")) +
  
  scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF",
                               "#8491B4FF", "#91D1C2FF")) +
  
  labs(
    x = "", 
    y = "No. of devAS"
  ) +
  
  theme(
    panel.grid = element_blank(), 
    panel.background = element_blank(), 
    text = element_text(size = 5, color = "#000000"), 
    line = element_line(linewidth = 0.5/lwd_pt), 
    
    axis.text = element_text(size = 5, color = "#000000"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(linewidth = 0.5/lwd_pt), 
    
    plot.title = element_text(size = 5, color = "#000000"), 
    
    legend.position = "none", 
    
    aspect.ratio = 1:1, # y:x
    plot.margin = unit(c(0.2, 0.1, 0, 0.1), "cm")
  )
p1

p1 <- ggplot(data = devAS_num_long_1) +
  
  geom_bar(mapping = aes(x = organ, fill = organ), width = 0.5) +
  
  scale_y_continuous(limits = c(0, 12000), expand = c(0, 0)) +
  scale_x_discrete(labels = c("testis", "brain", 
                              "cerebellum", "heart", "kidney", 
                              "liver", "ovary")) +
  
  scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF",
                               "#8491B4FF", "#91D1C2FF")) +
  
  labs(
    x = "", 
    y = "No. of devAS"
  ) +
  
  theme_prism(base_size = 7, 
              base_family = "Arial", 
              base_fontface = "plain", 
              base_line_size = 0.5/lwd_pt, 
              base_rect_size = 0.5/lwd_pt, 
              border = FALSE) +
  
  theme(
    legend.position = "none", 
    axis.text.x = element_text(angle = 45), 
    legend.key.size = unit(0.25, "cm"), 
    aspect.ratio = 1
  )
p1

# 大脑中不同类别不同pattern的AS占比 ---------------------------------------------------

# 在大脑当中CE占比最多，其中down又占比最多
human_brain_devAS_num_long <- devAS_num_long_1[c(devAS_num_long_1$organ == "pattern.brain"), ]
save(human_brain_devAS_num_long, 
     file = "02-analysis/19-Development_alternative_splicing/human_brain_devAS_num_long.Rdata")

# 按类型数目多少排序
human_brain_devAS_num_long$as.type <- factor(x = human_brain_devAS_num_long$as.type, 
                                             levels = c("CE", "RI", 
                                                        "AD", "AA", "complex"))

# p2 <- ggplot(data = human_brain_devAS_num_long) +
#   
#   geom_bar(mapping = aes(x = as.type, fill = pattern), 
#            width = 0.5) +
#   
#   scale_y_continuous(limits = c(0, 4000), expand = c(0, 0)) +
#   
#   scale_fill_manual(values = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"), 
#                     name = "") +
#   
#   labs(
#     x = "", 
#     y = "No. of different devAS in brain"
#   ) +
#   
#   theme(
#     panel.grid = element_blank(), 
#     panel.background = element_blank(), 
#     text = element_text(size = 5, color = "#000000"), 
#     line = element_line(linewidth = 0.5/lwd_pt), 
#     
#     axis.text = element_text(size = 5, color = "#000000"),
#     axis.line = element_line(linewidth = 0.5/lwd_pt), 
#     
#     plot.title = element_text(size = 5, color = "#000000"), 
#     
#     legend.key.size = unit(0.2, "cm"), 
#     legend.position = "top", 
#     
#     aspect.ratio = 1:1, 
#     plot.margin = unit(c(0.2, 0.1, 0, 0.1), "cm")
#   )
# p2

## prism style
p2 <- ggplot(data = human_brain_devAS_num_long) +
  
  geom_bar(mapping = aes(x = as.type, fill = pattern), 
           width = 0.5) +
  
  scale_y_continuous(limits = c(0, 4000), expand = c(0, 0)) +
  
  scale_fill_manual(values = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"), 
                    name = "") +
  
  labs(
    x = "", 
    y = "No. of different devAS in brain"
  ) +
  
  theme_prism(base_size = 7, 
              base_family = "Arial", 
              base_fontface = "plain", 
              base_line_size = 0.5/lwd_pt, 
              base_rect_size = 0.5/lwd_pt, 
              border = FALSE) +
  
  theme(
    legend.position = "top", 
    legend.key.size = unit(0.25, "cm"), 
    aspect.ratio = 1
  )
p2

# 组合图片 --------------------------------------------------------------------

cowplot::plot_grid(p1, p2, ncol = 2, 
                   align = c("hv"), axis = "bt", rel_widths = c(1, 1),
                   labels = c("a", "b"), 
                   label_size = 8)

pdf(file = "03-results/19-Development_alternative_splicing/pattern_type_in_brain.pdf", 
    width = 900/254, height = 600/254)
showtext_begin()

cowplot::plot_grid(p1, p2, ncol = 2, 
                   align = c("hv"), axis = "bt", rel_widths = c(1, 1),
                   labels = c("d", "e"), 
                   label_size = 8)

showtext_end()
dev.off()

































