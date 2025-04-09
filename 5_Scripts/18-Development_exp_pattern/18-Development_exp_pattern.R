# explore exp pattern of coding UCR genes using datasets data

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library(dbplyr)
library(pheatmap)
library(readxl)
library(stringr)
library(ggview)
library(ComplexHeatmap)

# human rpkm
human_rpkm <- read.table(file = "01-data/18-Development_exp_pattern/Human_rpkm.txt", 
                         sep = " ", header = TRUE)
human_rpkm <- human_rpkm[, c(1:54)]

# 获得在brain当中被识别为DDG的基因
load("D:/R_project/UCR_project/02-analysis/18-Development_dynamic_genes/coding_UCR_OrganDDGNum.Rdata")
coding_UCR_BrainDDG <- coding_UCR[, c(1:3)]
coding_UCR_BrainDDG <- coding_UCR_BrainDDG %>% 
  dplyr::filter(BrainDDG == 1)

coding_UCR_BrainDDG <- coding_UCR_BrainDDG[, c(1:2)]
coding_UCR_BrainDDG <- merge(coding_UCR_BrainDDG, human_rpkm, by.x = "V8", by.y = "Names")
colnames(coding_UCR_BrainDDG)[1:2] <- c("Ensembl_ID", "Gene_symbol")
coding_UCR_BrainDDG <- coding_UCR_BrainDDG[, -1]
rownames(coding_UCR_BrainDDG) <- coding_UCR_BrainDDG$Gene_symbol
coding_UCR_BrainDDG <- coding_UCR_BrainDDG[, -1]

heatmap_data <- log10(coding_UCR_BrainDDG + 1)
colnames(heatmap_data) <- factor(colnames(heatmap_data), levels = colnames(heatmap_data))

lwd_pt <- .pt*72.27/96

# breaks
bk <- c(seq(0, 1.9, by = 0.01), seq(2, 4, by = 0.01))

p1 <- pheatmap(mat = heatmap_data, 
               cellwidth = 10, cellheight = 10, 
               cluster_cols = FALSE, cluster_rows = TRUE, 
               angle_col = 90, main = "log10(rpkm + 1)", 
               color = c(colorRampPalette(colors = c("#083490","#f1f1f1"))(length(bk)/2), 
                         colorRampPalette(colors = c("#f1f1f1","#E60212"))(length(bk)/2)),
               legend_breaks=seq(0, 4, 0.4), 
               breaks = bk, 
               fontsize = 7)
p1

pdf(file = "03-results/18-Development_exp_pattern/coding_UCR_BrainDDG_heatmap.pdf", width = 10, height = 30)
p1
dev.off()

# 将coding UCR相关基因当中富集到alternative splicing通路的基因挑出来

AS_coding_UCR_gene <- read_xlsx(path = "01-data/10-GO_enrichment_analysis/coding_UCR_GO_results.xlsx", 
                                sheet = 2)
AS_coding_UCR_gene <- AS_coding_UCR_gene[, c(2, 6)]
AS_coding_UCR_gene$Term <- str_sub(AS_coding_UCR_gene$Term, 12, nchar(AS_coding_UCR_gene$Term))
AS_coding_UCR_gene <- AS_coding_UCR_gene[c(3, 7:10, 11), ]

# 分隔基因ID
gene_ids <- unlist(strsplit(AS_coding_UCR_gene$Genes, ","))
gene_ids <- unique(gene_ids)
gene_ids <- gsub(" ", "", gene_ids)

# 获得gene symbol
coding_UCR_in_AS_pathway <- filter(coding_UCR, V8 %in% gene_ids)
coding_UCR_in_AS_pathway <- coding_UCR_in_AS_pathway[!(coding_UCR_in_AS_pathway$V9 == "."), ]
coding_UCR_in_AS_pathway <- coding_UCR_in_AS_pathway[, c(1:2)]

coding_UCR_in_AS_pathway <- merge(coding_UCR_in_AS_pathway, human_rpkm, by.x = "V8", by.y = "Names")
coding_UCR_in_AS_pathway <- coding_UCR_in_AS_pathway[, -1]
rownames(coding_UCR_in_AS_pathway) <- coding_UCR_in_AS_pathway$V9
coding_UCR_in_AS_pathway <- coding_UCR_in_AS_pathway[, -1]

heatmap_coding_UCR_in_AS_pathway <- log10(coding_UCR_in_AS_pathway + 1)
colnames(heatmap_coding_UCR_in_AS_pathway) <- factor(x = colnames(heatmap_coding_UCR_in_AS_pathway), 
                                                     levels = colnames(heatmap_coding_UCR_in_AS_pathway), 
                                                     ordered = TRUE)
max(heatmap_coding_UCR_in_AS_pathway)
min(heatmap_coding_UCR_in_AS_pathway)

# breaks
bk <- c(seq(0, 1.24, by = 0.01), seq(1.25, 2.5, by = 0.01))

p2 <- pheatmap(mat = heatmap_coding_UCR_in_AS_pathway, 
               cellwidth = 5, cellheight = 5, 
               cluster_cols = FALSE, cluster_rows = TRUE, clustering_method = "average", 
               angle_col = 90, main = "log10(rpkm + 1)", 
               color = c(colorRampPalette(colors = c("#083490","#f1f1f1"))(length(bk)/2), 
                         colorRampPalette(colors = c("#f1f1f1","#E60212"))(length(bk)/2)),
               legend_breaks=seq(0, 2.5, 0.5), 
               breaks = bk, 
               fontsize = 5, 
               border = 0.5/lwd_pt, border_width = 0.5/lwd_pt)
p2

ggview(p2, width = 14, height = 10, unit = "cm", dpi = 1200)

pdf(file = "03-results/18-Development_exp_pattern/coding_UCR_in_AS_pathway_heatmap.pdf", width = 1400/254, height = 1000/254)
p2
dev.off()

gene_name <- rownames(heatmap_coding_UCR_in_AS_pathway)[p2$tree_row[["order"]]]
heatmap_coding_UCR_in_AS_pathway <- heatmap_coding_UCR_in_AS_pathway[gene_name, ]

p3 <- pheatmap(mat = heatmap_coding_UCR_in_AS_pathway, 
         cellwidth = 4, cellheight = 5, 
         cluster_cols = FALSE, cluster_rows = FALSE, clustering_method = "average", 
         angle_col = 90, main = "log10(rpkm + 1)", 
         color = c(colorRampPalette(colors = c("#083490","#f1f1f1"))(length(bk)/2), 
                   colorRampPalette(colors = c("#f1f1f1","#E60212"))(length(bk)/2)),
         legend_breaks=seq(0, 2.5, 0.5), 
         breaks = bk, 
         fontsize = 5, 
         border = 0.5/lwd_pt, border_width = 0.5/lwd_pt)

pdf(file = "03-results/18-Development_exp_pattern/notree_coding_UCR_in_AS_pathway_heatmap.pdf", width = 1400/254, height = 1000/254)
p3
dev.off()

# heatmap of a significant cluster
genesWithDecreasedExpression <- c("MATR3", "RBFOX2", "SFPQ", "HNRNPH1", "SRSF1", 
                                  "HNRNPU", "SRSF3", "SRSF6", "NR2F1", "ATP5MC2", 
                                  "TRA2A", "HNRNPM", "SRSF7", "PCBP2", "SRSF11", 
                                  "RBM39", "DDX3X", "PUM2", "OGT", "ZBTB18", 
                                  "CSDE1", "DDX5", "HNRNPDL", "HNRNPK", "ZFR", "G3BP2")

library(biomaRt)
library(curl)

# 利用useMart链接到人类的数据库，其他数据库使用listDatasets(ensembl)查看
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# attr为基因id目标类型，listAttributes(ensembl)展示其他id输出类型
attr <- c("ensembl_gene_id", "hgnc_symbol")

# listFilters(ensembl)展示其他id输入类型
ids <- getBM(attributes = attr, 
             filters = "hgnc_symbol", 
             values = genesWithDecreasedExpression, 
             mart = ensembl)

ids <- ids[!(ids$ensembl_gene_id %in% c("ENSG00000284254", "ENSG00000277564")), ]

# 将表达矩阵与ids合并
expDownregulatedGenes <- merge(ids, human_rpkm, by.x = "ensembl_gene_id", by.y = "Names")
expDownregulatedGenes <- expDownregulatedGenes[, -1]
rownames(expDownregulatedGenes) <- expDownregulatedGenes[, 1]
expDownregulatedGenes <- expDownregulatedGenes[, -1]

expDownregulatedGenes <- log10(expDownregulatedGenes + 1)
colnames(expDownregulatedGenes) <- factor(x = colnames(expDownregulatedGenes), 
                                          levels = colnames(expDownregulatedGenes), 
                                          ordered = TRUE)
max(expDownregulatedGenes)
min(expDownregulatedGenes)

# breaks
bk <- c(seq(0.8, 1.80, by = 0.01), seq(1.90, 2.8, by = 0.01))

p4 <- pheatmap(mat = expDownregulatedGenes, 
               cellwidth = 5, cellheight = 5, 
               cluster_cols = FALSE, cluster_rows = TRUE, clustering_method = "average", 
               angle_col = 90, main = "log10(rpkm + 1)", 
               color = c(colorRampPalette(colors = c("#083490","#f1f1f1"))(length(bk)/2), 
                         colorRampPalette(colors = c("#f1f1f1","#E60212"))(length(bk)/2)),
               legend_breaks=seq(0, 2.5, 0.5), 
               breaks = bk, 
               fontsize = 5, 
               border = 0.5/lwd_pt, border_width = 0.5/lwd_pt)
p4

ggview(p4, width = 14, height = 10, unit = "cm", dpi = 1200)

pdf(file = "03-results/18-Development_exp_pattern/heatmapOfDownRegulatedGene.pdf", width = 1400/254, height = 1000/254)
p4
dev.off()
