# 查看UCR和RF的基因的组织特异性
# tau计算公式来源：
# https://academic.oup.com/bioinformatics/article/21/5/650/220059

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(stringr)
library(ggview)

GTEx_exp <- read.table(file = "01-data/22-Tissue_specificity_across_different_organs/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", 
                       sep = "\t", skip = 2, header = TRUE)
GTEx_exp$Name <- str_sub(GTEx_exp$Name, 1, 15)
GTEx_exp <- GTEx_exp[, -2]

# UCR genes
coding_UCR_genes <- read.table(file = "02-analysis/16-New_classification/coding_UCR_ensembl_id.txt", 
                               sep = "\t")
ncRNA_UCR_genes <- read.table(file = "02-analysis/16-New_classification/ncRNA_UCR_ensembl_id.txt", 
                              sep = "\t")
UCR_genes <- rbind(coding_UCR_genes, ncRNA_UCR_genes)
UCR_genes <- merge(UCR_genes, GTEx_exp, by.x = "V1", by.y = "Name")

UCR_genes_tau <- UCR_genes
UCR_genes_tau <- unique(UCR_genes_tau)
rownames(UCR_genes_tau) <- UCR_genes_tau[, 1]
UCR_genes_tau <- UCR_genes_tau[, -1]

# 将每行数据除以每行的最大值
UCR_genes_tau_normalized <- t(apply(UCR_genes_tau, 1, function(x) x / max(x)))

# 转置结果以保持数据框形式
UCR_genes_tau_normalized <- as.data.frame(UCR_genes_tau_normalized)
UCR_genes_tau_normalized <- 1 - UCR_genes_tau_normalized

UCR_genes_tau_normalized <- mutate(UCR_genes_tau_normalized, 
                                   tau = rowSums(UCR_genes_tau_normalized)/53)

# RF genes
coding_rf_genes <- read.table(file = "02-analysis/16-New_classification/01-rf_form_protein_coding_gene.bed", 
                              sep = "\t", header = FALSE)
coding_rf_genes <- coding_rf_genes[, c(8, 10)]

ncRNA_rf_genes <- read.table(file = "02-analysis/16-New_classification/02-rf_from_Non_Coding_RNA.bed", 
                             sep = "\t", header = FALSE)
ncRNA_rf_genes <- ncRNA_rf_genes[, c(8, 10)]
rf_genes <- rbind(coding_rf_genes, ncRNA_rf_genes)
rf_genes <- merge(rf_genes, GTEx_exp, by.x = "V8", by.y = "Name")
rf_genes <- rf_genes[, -2]

rf_genes_tau <- rf_genes
rf_genes_tau <- unique(rf_genes_tau)
rf_genes_tau <- rf_genes_tau[!(c(rf_genes_tau$V8 == "ENSG00000223511" & rf_genes_tau$Adipose...Subcutaneous == 0)), ]
rownames(rf_genes_tau) <- rf_genes_tau[, 1]
rf_genes_tau <- rf_genes_tau[, -1]

# 将每行数据除以每行的最大值
rf_genes_tau_normalized <- t(apply(rf_genes_tau, 1, function(x) x / max(x)))

# 转置结果以保持数据框形式
rf_genes_tau_normalized <- as.data.frame(rf_genes_tau_normalized)
rf_genes_tau_normalized <- 1 - rf_genes_tau_normalized

rf_genes_tau_normalized <- mutate(rf_genes_tau_normalized, 
                                   tau = rowSums(rf_genes_tau_normalized)/53)

save(UCR_genes_tau_normalized, 
     file = "02-analysis/22-Tissue_specificity_across_different_organs/UCR_genes_tau_normalized.Rdata")
save(rf_genes_tau_normalized, 
     file = "02-analysis/22-Tissue_specificity_across_different_organs/rf_genes_tau_normalized.Rdata")

# 绘制统计图 -------------------------------------------------------------------

library(ggplot2)
library(sysfonts)
library(showtext)
library(ggprism)

# UCR基因和RF基因比较
UCR_genes_tau_normalized <- UCR_genes_tau_normalized[!(UCR_genes_tau_normalized$tau == "NaN"), ]
rf_genes_tau_normalized <- rf_genes_tau_normalized[!(rf_genes_tau_normalized$tau == "NaN"), ]

p1 <- ggplot(mapping = aes(x = tau)) +
  geom_density(data = UCR_genes_tau_normalized) +
  geom_density(data = rf_genes_tau_normalized)
p1

wilcox.test(UCR_genes_tau_normalized$tau, 
            rf_genes_tau_normalized$tau)

# UCR基因和all基因比较
GTEx_exp_tau <- GTEx_exp

# 计算第2-55列数值的和
row_sums <- rowSums(GTEx_exp_tau[, 2:55])

# 删除和为0的行
GTEx_exp_tau_filtered <- GTEx_exp_tau[row_sums != 0, ]

rownames(GTEx_exp_tau_filtered) <- GTEx_exp_tau_filtered[, 1]
GTEx_exp_tau_filtered <- GTEx_exp_tau_filtered[, -1]

# 将每行数据除以每行的最大值
GTEx_exp_tau_normalized <- t(apply(GTEx_exp_tau_filtered, 1, function(x) x / max(x)))

# 转置结果以保持数据框形式
GTEx_exp_tau_normalized <- as.data.frame(GTEx_exp_tau_normalized)
GTEx_exp_tau_normalized <- 1 - GTEx_exp_tau_normalized

GTEx_exp_tau_normalized <- mutate(GTEx_exp_tau_normalized, 
                                  tau = rowSums(GTEx_exp_tau_normalized)/53)
save(GTEx_exp_tau_normalized, 
     file = "02-analysis/22-Tissue_specificity_across_different_organs/GTEx_exp_tau_normalized.Rdata")

lwd_pt <- .pt*72.27/96

font_add(family = "Arial", regular = "C:/Windows/Fonts/arial.ttf")
font_families()

pdf(file = "03-results/22-Tissue_specificity_across_different_organs/Tissue specificity(τ) across different organs_2.pdf", 
    width = 480/254, height = 600/254)
showtext_begin()

p2 <- ggplot() +
  geom_density(data = UCR_genes_tau_normalized, 
               aes(x = tau, fill = "UCR Genes"),  # Mapping fill to a label
               linewidth = 0.5/lwd_pt, 
               alpha = 0.7, 
               position = position_nudge(x = 0, y = 1)) +
  geom_density(data = GTEx_exp_tau_normalized, 
               aes(x = tau, fill = "GTEx Exp"),  # Mapping fill to a label
               linewidth = 0.5/lwd_pt,
               alpha = 0.7) +
  
  scale_fill_manual(values = c("UCR Genes" = "#E64B35FF", "GTEx Exp" = "#4DBBD5FF")) +  # Define colors for the legend
  scale_x_continuous(limits = c(0.25, 1), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  
  labs(x = "Tissue specificity (τ)
  across different organs", 
       y = "", 
       title = "", 
       fill = "Dataset") +  # Add a label for the legend
  
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

p3 <- p2 +
  annotate("text", x = 0.5, y = 6.5, label = expression(paste("p value < 2.2 \u00d7 10"^{-16})), size = 5/lwd_pt)
p3

plotTissueSpecificity <- p3

showtext_end()
dev.off()

# p4 <- p3 +
#   annotate("text", x = 0.5, y = 5.5, size = 5/lwd_pt, 
#            label = expression(paste(τ == sum((1-x[i]), i==1, n) / (N-1))))
# p4
# 
# ggview(p4, width = 4.8, height = 6, units = "cm", dpi = 1200)
# 
# ggsave(filename = "03-results/22-Tissue_specificity_across_different_organs/Tissue specificity(τ) across different organs.tiff", 
#        width = 4.8, height = 6, units = "cm", dpi = 1200)
# 
# pdf(file = "03-results/22-Tissue_specificity_across_different_organs/Tissue specificity(τ) across different organs_1.pdf", 
#     width = 480/254, height = 600/254)
# p4
# dev.off()

wilcox.test(UCR_genes_tau_normalized$tau, 
            GTEx_exp_tau_normalized$tau)
# Wilcoxon rank sum test with continuity correction
# 
# data:  UCR_genes_tau_normalized$tau and GTEx_exp_tau_normalized$tau
# W = 4001854, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0











