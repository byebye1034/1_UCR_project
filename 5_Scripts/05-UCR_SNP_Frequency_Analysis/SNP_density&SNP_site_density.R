# SNP density & SNP site density

setwd(dir = "D:/R_project/UCR_project")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(ggpubr)
library(patchwork)
library(ggsci)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(showtext)
library(gginnards) # ���������Ա�ǵ������С
library(export) # ������AI�����ֿ���ֱ�ӱ༭��pdf�ļ�
library(ggview)
library(ggprism)

# ��UCR������������ע�� -----------------------------------------------------------

# �����ͼƬ�???05-UCR_Classification_and_Closest���ļ�����

# ����SNP��λ��ĳ���ռUCR�ܳ��ȵı��� ---------------------------------------------------

load("02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/passed_UCR_SNPs.Rdata")
load("02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/passed_UCR_right_SNPs.Rdata")
load("02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/passed_UCR_left_SNPs.Rdata")
load("02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/passed_rf_SNPs.Rdata")

UCR_location <- read.table(file = "01-data/UCR_raw/UCR_location.txt", header = T)
random_fragments_location <- read.table(file = "01-data/UCR_raw/random_fragments_location.txt", header = T)

# UCR��SNP����Ŀ

# ʹ��group_by��summarize���м���
UCR_SNP_SUM <- passed_UCR_SNPs %>%
  group_by(UCR_name = list_type) %>%
  summarize(count = n())
UCR_SNP_SUM <- merge(UCR_SNP_SUM, UCR_location[, c(4:5)], by = "UCR_name")
UCR_SNP_SUM <- UCR_SNP_SUM %>% 
  mutate(proportion = count/length)
UCR_SNP_SUM <- UCR_SNP_SUM[-c(89, 228), ]
save(UCR_SNP_SUM, file = "02-analysis/05-UCR_SNP_Frequency_Analysis/SNP_density/UCR_SNP_density.Rdata")

## left flank
UCR_left_SNP_SUM <- passed_UCR_left_SNPs %>%
  group_by(UCR_name = list_type) %>%
  summarize(count = n())
UCR_left_SNP_SUM <- merge(UCR_left_SNP_SUM, UCR_location[, c(4:5)], by = "UCR_name")
UCR_left_SNP_SUM <- UCR_left_SNP_SUM %>% 
  mutate(proportion = count/length)
UCR_left_SNP_SUM <- UCR_left_SNP_SUM[-c(89, 228), ]
save(UCR_left_SNP_SUM, file = "02-analysis/05-UCR_SNP_Frequency_Analysis/SNP_density/UCR_left_SNP_density.Rdata")

## right flank
UCR_right_SNP_SUM <- passed_UCR_right_SNPs %>%
  group_by(UCR_name = list_type) %>%
  summarize(count = n())
UCR_right_SNP_SUM <- merge(UCR_right_SNP_SUM, UCR_location[, c(4:5)], by = "UCR_name")
UCR_right_SNP_SUM <- UCR_right_SNP_SUM %>% 
  mutate(proportion = count/length)
UCR_right_SNP_SUM <- UCR_right_SNP_SUM[-c(89, 228), ]
save(UCR_right_SNP_SUM, file = "02-analysis/05-UCR_SNP_Frequency_Analysis/SNP_density/UCR_right_SNP_density.Rdata")

## random fragments
random_fragments_SNP_SUM <- passed_rf_SNPs %>% 
  group_by(RF_name = list_type) %>% 
  summarise(count = n())
random_fragments_SNP_SUM <- merge(random_fragments_SNP_SUM, random_fragments_location[, c(4, 6)], 
                                  by.x = "RF_name", by.y = "rf_name")
random_fragments_SNP_SUM <- random_fragments_SNP_SUM %>% 
  mutate(proportion = count/Length)
random_fragments_SNP_SUM <- random_fragments_SNP_SUM[-c(89, 228), ]
save(random_fragments_SNP_SUM, file = "02-analysis/05-UCR_SNP_Frequency_Analysis/SNP_density/random_fragments_SNP_density.Rdata")

SNP_density <- cbind(UCR_SNP_SUM$proportion, UCR_left_SNP_SUM$proportion, 
                     UCR_right_SNP_SUM$proportion, random_fragments_SNP_SUM$proportion)
SNP_density <- as.data.frame(SNP_density)
SNP_density[c(437:479), 4] <- 0
colnames(SNP_density) <- c("UCR", "UCR_left", "UCR_right", "random_fragment")
save(SNP_density, file = "02-analysis/05-UCR_SNP_Frequency_Analysis/SNP_density/SNP_density.Rdata")

# ʹ��melt����������ת��Ϊ����ʽ
long_format_data_of_SNP_density <- melt(SNP_density, id.vars = NULL, variable.name = "Group", value.name = "SNP_density")
long_format_data_of_SNP_density$Group <- as.character(long_format_data_of_SNP_density$Group)
long_format_data_of_SNP_density$Group <- factor(long_format_data_of_SNP_density$Group, levels = c("UCR", "UCR_left", "UCR_right", "random_fragment"))

# SNP_density -------------------------------------------------------------

# SNP_density
lwd_pt <- .pt*72.27/96
compare_list <- list(c("UCR", "UCR_left"), 
                     c("UCR", "UCR_right"), 
                     c("UCR", "random_fragment"))

p1 <- ggplot(data = long_format_data_of_SNP_density, mapping = aes(x = Group, y = SNP_density)) +
  geom_violin(mapping = aes(color = Group), linewidth = 0.5/lwd_pt) +
  geom_boxplot(mapping = aes(color = Group),
               linewidth = 0.5/lwd_pt,
               outlier.size = 0.25, 
               width = 0.5) +
  
  scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF"), 
                     name = "", 
                     breaks = c("UCR", "UCR_left", 
                                "UCR_right", "random_fragment")) +
  
  scale_y_continuous(breaks = seq(0, 1.0, by = 0.2), 
                     limits = c(0, 1), 
                     expand = c(0, 0)) +
  
  labs(x = "", 
       y = "SNP density") +
  
  guides(color = guide_legend(nrow = 2, 
                              ncol = 2)) + # ����ͼ��Ϊ��������
  
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

ggview::ggview(last_plot(), width = 4.8, height = 6, units = "cm", dpi = 1200)

p2 <- p1 + stat_compare_means(
  comparisons = compare_list, 
  method = "wilcox.test", 
  label = "p.format", 
  bracket.size = 0.5/lwd_pt, 
  step.increase = 0.20
)
p2

p2$layers[[which_layers(p2, "GeomSignif")]]$aes_params$textsize <- 5/lwd_pt
p2

ggsave(filename = "03-results/05-SNP_density/SNP_density_violin_boxplot.tiff", 
       width = 4.8, height = 6, units = "cm", dpi = 1200)

pdf(file = "03-results/05-SNP_density/SNP_density_violin_boxplot.pdf", width = 480/254, height = 600/254)
p2
dev.off()

# SNP_site_density --------------------------------------------------------

# SNP_site_density

rm(list = ls())

load("02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/passed_UCR_SNPs.Rdata")
load("02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/passed_UCR_right_SNPs.Rdata")
load("02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/passed_UCR_left_SNPs.Rdata")
load("02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/passed_rf_SNPs.Rdata")

UCR_location <- read.table(file = "01-data/UCR_raw/UCR_location.txt", header = T)
random_fragments_location <- read.table(file = "01-data/UCR_raw/random_fragments_location.txt", header = T)

# UCR SNP site
UCR_SNP_site <- passed_UCR_SNPs[, c(1, 7)]
UCR_SNP_site <- UCR_SNP_site[!duplicated(UCR_SNP_site$chrom_start), ]
UCR_SNP_site_sum <- UCR_SNP_site %>% 
  group_by(UCR_name = list_type) %>% 
  summarise(count = n())
UCR_SNP_site_sum <- merge(UCR_SNP_site_sum, UCR_location[, c(4:5)], by = "UCR_name")
UCR_SNP_site_sum <- UCR_SNP_site_sum %>% 
  mutate(site_density = count/length)
UCR_SNP_site_sum <- UCR_SNP_site_sum[c(-89, -228), ]
save(UCR_SNP_site_sum, file = "02-analysis/05-UCR_SNP_Frequency_Analysis/SNP_site_density/UCR_SNP_site_sum.Rdata")

# UCR left SNP site
UCR_left_SNP_site <- passed_UCR_left_SNPs[, c(1, 7)]
UCR_left_SNP_site <- UCR_left_SNP_site[!duplicated(UCR_left_SNP_site$chrom_start), ]
UCR_left_SNP_site_sum <- UCR_left_SNP_site %>% 
  group_by(UCR_name = list_type) %>% 
  summarise(site_count = n())
UCR_left_SNP_site_sum <- merge(UCR_left_SNP_site_sum, UCR_location[, c(4:5)], by = "UCR_name")
UCR_left_SNP_site_sum <- UCR_left_SNP_site_sum %>% 
  mutate(site_density = site_count/length)
UCR_left_SNP_site_sum <- UCR_left_SNP_site_sum[c(-89, -228), ]
save(UCR_left_SNP_site_sum, file = "02-analysis/05-UCR_SNP_Frequency_Analysis/SNP_site_density/UCR_left_SNP_site_sum.Rdata")

# UCR right SNP site
UCR_right_SNP_site <- passed_UCR_right_SNPs[, c(1, 7)]
UCR_right_SNP_site <- UCR_right_SNP_site[!duplicated(UCR_right_SNP_site$chrom_start), ]
UCR_right_SNP_site_sum <- UCR_right_SNP_site %>% 
  group_by(UCR_name = list_type) %>% 
  summarise(site_count = n())
UCR_right_SNP_site_sum <- merge(UCR_right_SNP_site_sum, UCR_location[, c(4:5)], by = "UCR_name")
UCR_right_SNP_site_sum <- UCR_right_SNP_site_sum %>% 
  mutate(site_density = site_count/length)
UCR_right_SNP_site_sum <- UCR_right_SNP_site_sum[c(-89, -228), ]
save(UCR_right_SNP_site_sum, file = "02-analysis/05-UCR_SNP_Frequency_Analysis/SNP_site_density/UCR_right_SNP_site_sum.Rdata")

# random fragments
rf_SNP_site <- passed_rf_SNPs[, c(1, 7)]
rf_SNP_site <- rf_SNP_site[!duplicated(rf_SNP_site$chrom_start), ]
rf_SNP_site_sum <- rf_SNP_site %>% 
  group_by(rf_name = list_type) %>% 
  summarise(site_count = n())
rf_SNP_site_sum <- merge(rf_SNP_site_sum, random_fragments_location[, c(6, 4)], by = "rf_name")
rf_SNP_site_sum <- rf_SNP_site_sum %>% 
  mutate(site_density = site_count/Length)
rf_SNP_site_sum <- rf_SNP_site_sum[c(-89, -228), ]
save(rf_SNP_site_sum, file = "02-analysis/05-UCR_SNP_Frequency_Analysis/SNP_site_density/rf_SNP_site_sum.Rdata")

SNP_site_density <- cbind(UCR_SNP_site_sum$site_density, UCR_left_SNP_site_sum$site_density, 
                          UCR_right_SNP_site_sum$site_density, rf_SNP_site_sum$site_density)
SNP_site_density <- as.data.frame(SNP_site_density)
SNP_site_density[c(437:479), 4] <- 0
colnames(SNP_site_density) <- c("UCR", "UCR_left", "UCR_right", "random_fragment")
save(SNP_site_density, file = "02-analysis/05-UCR_SNP_Frequency_Analysis/SNP_site_density/SNP_site_density.Rdata")

# ʹ��melt����ת���ɳ�����
long_format_data_of_SNP_site_density <- melt(SNP_site_density, id.vars = NULL, variable.name = "Group", value.name = "SNP_site_density")
long_format_data_of_SNP_site_density$Group <- as.character(long_format_data_of_SNP_site_density$Group)
long_format_data_of_SNP_site_density$Group <- factor(long_format_data_of_SNP_site_density$Group, levels = c("UCR", "UCR_left", "UCR_right", "random_fragment"))

# SNP site density --------------------------------------------------------

lwd_pt <- .pt*72.27/96
compare_list <- list(c("UCR", "UCR_left"), 
                     c("UCR", "UCR_right"), 
                     c("UCR", "random_fragment"))

p3 <- ggplot(data = long_format_data_of_SNP_site_density, mapping = aes(x = Group, y = SNP_site_density)) +
  geom_violin(mapping = aes(color = Group), linewidth = 0.5/lwd_pt) +
  geom_boxplot(mapping = aes(color = Group), 
               linewidth = 0.5/lwd_pt, 
               outlier.size = 0.25, 
               width = 0.5) +
  
  scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF"), 
                     name = "", 
                     breaks = c("UCR", "UCR_left", 
                                "UCR_right", "random_fragment")) +
  
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.3), 
                     limits = c(0, 0.9), 
                     expand = c(0, 0)) +
  
  labs(x = "", 
       y = "SNP site density") +
  
  guides(color = guide_legend(nrow = 2, 
                              ncol = 2)) +
  
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
p3

p4 <- p3 + stat_compare_means(
  comparisons = compare_list, 
  method = "wilcox.test", 
  label = "p.format", 
  bracket.size = 0.5/lwd_pt, 
  step.increase = 0.2
)
p4

p4$layers[[which_layers(p4, "GeomSignif")]]$aes_params$textsize <- 5/lwd_pt
p4

ggview(p4, width = 4.8, height = 6, units = "cm", dpi = 1200)

ggsave(filename = "03-results/05-SNP_density/SNP_site_density_violin_boxplot.tiff", 
       width = 4.8, height = 6, units = "cm", dpi = 1200)

pdf(file = "03-results/05-SNP_density/SNP_site_density_violin_boxplot.pdf", 
    width = 480/254, height = 600/254)
p4
dev.off()

# 20230604 ��SNP density��SNP site densityͼ��???

library(cowplot)
library(sysfonts)
library(showtext)

font_add(family = "Arial", regular = "C:/Windows/Fonts/arial.ttf")
font_families()

cowplot::plot_grid(p2, p4, 
                   align = c("hv"), axis = c("tb"), 
                   nrow = 1, ncol = 2, rel_widths = c(1, 1), 
                   labels = c("c", "d"), label_size = 8)
ggview(plot = ggplot2::last_plot(), width = 9, height = 6, units = "cm", dpi = 1200)

pdf(file = "03-results/05-SNP_density/SNP&SNP_site_density.pdf", width = 900/254, height = 600/254)
showtext_begin()

cowplot::plot_grid(p2, p4, 
                   align = c("hv"), axis = c("tb"), 
                   nrow = 1, ncol = 2, rel_widths = c(1, 1), 
                   labels = c("a", "b"), label_size = 8)

showtext_end()
dev.off()

# С����ͼ + ����ͼ --------------------------------------------------------------

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(ggprism)
library(reshape2)
library(ggpubr) # stat_compare_means����Pֵ
library(gginnards) # ���������Ա�ǵ������С
library(export) # ������AI�����ֿ���ֱ�ӱ༭��pdf�ļ�
library(ggsci) # ʹ��npg��ɫ
library(ggview)

# ɾ��uc.18(89��)��uc.304(228��)
load("D:/R_project/UCR_project/02-analysis/05-UCR_SNP_Frequency_Analysis/SNP_proportion_SUM.Rdata")
SNP_proportion_SUM <- SNP_proportion_SUM[c(-89, -228), ]

# ת���ɳ�����
SNP_proportion_SUM_long <- melt(SNP_proportion_SUM, id.vars = NULL, variable.name = "Group", value.name = "SNP_Ratio")
SNP_proportion_SUM_long$Group <- factor(SNP_proportion_SUM_long$Group, levels = c("UCR", "UCR_left", "UCR_right", "random_fragment"))

lwd_pt <- .pt*72.27/96
custom_palette <- c("#000000", "#E60212", "#083490", "#751384")
compare_list <- list(c("UCR", "UCR_left"), 
                     c("UCR", "UCR_right"), 
                     c("UCR", "random_fragment"))

p1 <- ggplot(data = SNP_proportion_SUM_long, mapping = aes(x = Group, y = SNP_Ratio)) +
  geom_violin(mapping = aes(color = Group), linewidth = 0.5/lwd_pt) +                                           # �����colorָ���Ǹ���˭������ɫ
  geom_boxplot(mapping = aes(color = Group), 
               linewidth = 0.5/lwd_pt,                                          # ��������ͼ��������ϸ
               outlier.size = 0.25, 
               width = 0.5) +                                           # ������Ⱥֵ��Ĵ�???
  
  scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF"), 
                     name = "", 
                     breaks = c("UCR", "UCR_left", 
                                "UCR_right", "random_fragment")) +
  
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.3), limits = c(0, 0.9), 
                     expand = c(0, 0)) +
  
  guides(color = guide_legend(nrow = 2, 
                              ncol = 2)) +
  
  labs(x = "", 
       y = "SNP density") +
  
  theme(
    panel.background = element_blank(), 
    panel.grid = element_blank(), 
    text = element_text(size = 7), 
    line = element_line(linewidth = 0.5/lwd_pt), 
    
    axis.text = element_text(size = 7, color = "#000000"), 
    axis.line = element_line(linewidth = 0.5/lwd_pt), 
    axis.text.x = element_blank(), 
    
    aspect.ratio = 1,                                                           # ���ֻ�ͼ����Ϊ1��1
    
    legend.position = "top", 
    legend.key.size = unit(0.2, "cm")
    
    # legend.key.width = unit(0.25, "cm"),  # ����ͼ�����Ŀ���
    # legend.key.height = unit(0.25, "cm"),  # ����ͼ�����ĸ߶�
    # legend.spacing.y = unit(1, "cm"), 
    # legend.text = element_text(size = 7),  # ����ͼ���ı��Ĵ�С
    # legend.title = element_text(size = 7)  # ����ͼ������Ĵ�???
  )

p1

p2 <- p1 + stat_compare_means(
  comparisons = compare_list, 
  method = "wilcox.test", 
  label = "p.signif", 
  bracket.size = 0.5/lwd_pt, 
  step.increase = 0.15
)
p2

p2$layers[[which_layers(p2, "GeomSignif")]]$aes_params$textsize <- 5/lwd_pt

ggview(p2, width = 4.8, height = 6, units = "cm", dpi = 1200)

pdf(file = "03-results/figures/UCR������SNP��ռ����С����ͼ+����ͼ.pdf", width = 480/254, height = 480/254)
p2
dev.off()


