# coding UCR相关基因参与的phase separated condensates

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(readxl)
library(biomaRt)
library(curl)

coding_UCR_ensembl_id <- read.table(file = "02-analysis/16-New_classification/coding_UCR_ensembl_id.txt")

# 利用useMart链接到人类的数据库，其他数据库使用listDatasets(ensembl)查看
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# attr为基因id目标类型，listAttributes(ensembl)展示其他id输出类型
attr <- c("ensembl_gene_id", "uniprotswissprot")

# listFilters(ensembl)展示其他id输入类型
ids <- getBM(attributes = attr, 
             filters = "ensembl_gene_id", 
             values = coding_UCR_ensembl_id$V1, 
             mart = ensembl)
coding_UCR_uniprot_id <- ids

# Human_phase_separated_condensates
Human_phase_separated_condensates <- 
  read_xlsx(path = "01-data/23-Phase_separated_condensates/1-s2.0-S2211124723008227-mmc3.xlsx", 
            sheet = 6, skip = 3)
colnames(Human_phase_separated_condensates) <- c("Uniprot_accession", 
                                                 "Condensates", 
                                                 "Major_biological_process")

coding_UCR_genes_PSC <- merge(coding_UCR_uniprot_id, Human_phase_separated_condensates, 
                              by.x = "uniprotswissprot", by.y = "Uniprot_accession")
save(coding_UCR_genes_PSC, 
     file = "02-analysis/23-Phase_separated_condensates/coding_UCR_genes_PSC.Rdata")

# coding RF
coding_RF_ensembl_id <- read.table(file = "02-analysis/16-New_classification/01-rf_form_protein_coding_gene.bed")
coding_RF_ensembl_id <- coding_RF_ensembl_id[, c(8:9)]

ids <- getBM(attributes = attr, 
             filters = "ensembl_gene_id", 
             values = coding_RF_ensembl_id$V8, 
             mart = ensembl)
coding_UCR_uniprot_id <- ids
coding_RF_ensembl_id <- ids

coding_RF_genes_PSC <- merge(coding_RF_ensembl_id, Human_phase_separated_condensates, 
                             by.x = "uniprotswissprot", by.y = "Uniprot_accession")
save(coding_RF_genes_PSC, 
     file = "02-analysis/23-Phase_separated_condensates/coding_RF_genes_PSC.Rdata")

# 绘制结果图 -------------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(UpSetR)
library(ggplot2)
library(forcats) # 颠倒因子顺序
library(ggsci)

load("D:/R_project/UCR_project/02-analysis/23-Phase_separated_condensates/coding_UCR_genes_PSC.Rdata")

coding_UCR_genes_PSC <- mutate(coding_UCR_genes_PSC, Type = "UCR")

coding_UCR_genes_PSC$Condensates <- factor(x = coding_UCR_genes_PSC$Condensates, 
                      levels = c("Stress granule", "P-body", 
                                 "Nucleolus", "Postsynaptic density", 
                                 "Nuclear speckle", "Paraspeckle", 
                                 "Droplet", "Centrosome/Spindle pole body", 
                                 "Others", "PML nuclear body", 
                                 "Nuclear stress body", "RNP granules", 
                                 "Transcription factories", "Nuclear protein granule", 
                                 "Sam68 nuclear body", "Neuronal granule", 
                                 "Nuclear body", "Cajal body"), 
                      ordered = TRUE)
coding_UCR_genes_PSC$Condensates <- fct_rev(coding_UCR_genes_PSC$Condensates)

lwd_pt <- .pt*72.27/96

p1 <- ggplot(data = coding_UCR_genes_PSC) +
  geom_histogram(mapping = aes(x = Condensates), stat = "count", 
                 binwidth = 7, alpha = 0.7, fill = "#E64B35FF") +
  
  scale_y_continuous(limits = c(0, 35), expand = c(0, 0)) +
  
  labs(x = "", 
       y = "", 
       title = "") +
  
  theme(
    panel.grid = element_blank(), 
    panel.background = element_blank(), 
    text = element_text(size = 7, color = "black"), 
    line = element_line(linewidth = 0.5/lwd_pt), 
    
    axis.text = element_text(size = 7, color = "black"), 
    axis.line = element_line(linewidth = 0.5/lwd_pt),
    axis.text.y = element_text(hjust = 0),  # 设置 Y 轴标签左对齐
    
    aspect.ratio = 1:1.2
  ) +
  coord_flip()
p1

pdf(file = "03-results/23-Phase_separated_condensates/histogram.pdf", width = 1080/254, height = 720/254)
p1
dev.off()

# sankey diagram 桑基图 ------------------------------------------------------

library(ggalluvial)






