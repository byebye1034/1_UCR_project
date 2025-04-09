# 根据UCR的基因组位置对UCR进行分类，并寻找最近的基因用于下游的富集分析

rm(list = ls())
setwd(dir = "D:/R_project/UCR_project/")

library(tidyverse)

UCR_with_overlap_ratio <- read.table(file = "data/0A-EvolutionAnalysisData/UCR classification/UCR_with_overlap_ratio.bed", 
                                     sep = "\t")
colnames(UCR_with_overlap_ratio) <- c("UCR_chr", "UCR_start", "UCR_end", "UCR_name", 
                                      "Exon_chr", "Exon_start", "Exon_end", "Gene_id", "Gene_symbol", "Overlap")
UCR_with_overlap_ratio <- UCR_with_overlap_ratio %>% 
  mutate(Overlap_in_UCR = Overlap/(UCR_end - UCR_start)) %>% 
  mutate(Overlap_in_Exon = Overlap/(Exon_end - Exon_start))
UCR_with_overlap_ratio$Type <- ifelse(UCR_with_overlap_ratio$Overlap_in_UCR == 1, "exonic", 
                                      ifelse(UCR_with_overlap_ratio$Overlap_in_Exon == 1, "exon_containing", "partly_exonic"))
UCR_classification <- UCR_with_overlap_ratio[, c(1:4, 8, 9, 13)]
UCR_classification <- unique(UCR_classification)
UCR_classification <- UCR_classification %>% 
  group_by(UCR_name) %>% 
  summarise(New_type = ifelse(n_distinct(Type) > 1, paste(Type, collapse = "/"), first(Type)))%>% 
  right_join(UCR_classification, by = "UCR_name") %>% 
  select(-Type) %>% 
  rename(Type = New_type)
UCR_classification <- unique(UCR_classification)
save(UCR_classification, file = "data/0A-EvolutionAnalysisData/UCR classification/UCR_noninter_nonintronic_classification.Rdata")

x <- UCR_classification[, c(-6:-7)]
x <- unique(x)
table(x$Type)
# noninter_nonintronic 142
# exon_containing 23
# exonic 27
# partly exonic 27
# multiple 65

# 绘制文献中不同类型UCR数目的柱形统计图 ----------------------------------------------------

rm(list = ls())
setwd(dir = "D:/R_project/UCR_project")

library(tidyverse)
library(ggplot2)
library(ggprism)
library(patchwork)

UCR_class <- data.frame(UCR_type = c("intergenic", "intronic", "exonic", "exon containing", "partly exonic", "multiple"), 
                        UCR_number = c(98, 241, 27, 23, 27, 65))
UCR_class$UCR_type <- factor(UCR_class$UCR_type, levels = c("intergenic", "intronic", "exonic", "exon containing", "partly exonic", "multiple"))

custom_palette <- c("#000000", "#E60212", "#083490", "#751384", "#007A34", "#D85F00")

lwd_pt <- .pt*72.27/96

ggplot(data = UCR_class, mapping = aes(x = UCR_type, y = UCR_number, color = UCR_type)) +
  geom_col(width = 0.8, fill = NA, size = 0.5/lwd_pt) +                         # width:柱子宽度 size:边框的粗细
  
  scale_color_manual(values = custom_palette) +
  scale_y_continuous(limits = c(0, 250), expand = c(0, 0)) +
  
  theme_prism(palette = "black_and_white", 
              base_size = 7, 
              base_family = "serif", 
              base_fontface = "bold", 
              base_line_size = 0.5/lwd_pt, 
              base_rect_size = 0.5/lwd_pt, 
              axis_text_angle = 45, 
              border = FALSE) +
  
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(), 
    legend.key.size = unit(0.25, "cm"),                                         # 调整图例大小
    legend.spacing.y = unit(1, "cm")                                            # 调整图例之间的垂直距离
  )
  

pdf(file = "D:/R_project/UCR_project/03-results/figures/重新分类后的不同类型的UCR的数目统计柱形图.pdf", width = 720/254, height = 480/254)

dev.off()

# 绘制分类饼图 ------------------------------------------------------------------

# uc.18是intronic
# uc.304是exon containing/partly exonic

# 多层饼图在excel绘制（20230306）













































