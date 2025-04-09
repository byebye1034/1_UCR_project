# 用数据库的数据探索coding UCR相关基因的表达模式

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(dbplyr)
library(pheatmap)
library(readxl)
library(stringr)
library(tidyverse)
library(biomaRt)
library(curl)

# human brain rpkm
human_rpkm <- read.table(file = "01-data/18-Development_exp_pattern/Human_rpkm.txt", 
                         sep = " ", header = TRUE)
human_rpkm <- human_rpkm[, c(1:54)]

# human brain ucr related coding genes rpkm
coding_UCR_related_genes <- read.table(
  file = "02-analysis/16-New_classification/coding_UCR_ensembl_id.txt", 
  sep = "\t"
  )
coding_UCR_related_genes <- unique(coding_UCR_related_genes)
coding_rpkm <- merge(coding_UCR_related_genes, 
                     human_rpkm, 
                     by.x = "V1", 
                     by.y = "Names")

# data prepare
heatmap_data <- coding_rpkm
rownames(heatmap_data) <- heatmap_data[, 1]
heatmap_data <- heatmap_data[, -1]
heatmap_data <- log10(heatmap_data + 1)

colnames(heatmap_data) <- factor(x = colnames(heatmap_data), levels = colnames(heatmap_data))

lwd_pt <- .pt*72.27/96

max(heatmap_data) # 3.659032
min(heatmap_data) # 0

# breaks
bk <- c(seq(0, 1.9, by = 0.01), seq(2, 4, by = 0.01))

p1 <- pheatmap(mat = heatmap_data, 
               cellwidth = 6, cellheight = 2, 
               cluster_cols = FALSE, cluster_rows = TRUE, 
               show_rownames = FALSE, 
               angle_col = 90, main = "log10(rpkm + 1)", 
               color = c(colorRampPalette(colors = c("#083490","#f1f1f1"))(length(bk)/2), 
                         colorRampPalette(colors = c("#f1f1f1","#E60212"))(length(bk)/2)),
               legend_breaks=seq(0, 4, 0.4), 
               breaks = bk, 
               fontsize = 7)
p1

pdf(file = "03-results/18-Development_exp_pattern/All_coding_UCR_related_coding_exp_heatmap.pdf", 
    width = 1440/254, height = 1920/254)
p1
dev.off()

# identity module
geneid <- rownames(heatmap_data)

# 利用useMart链接到人类的数据库，其他数据库使用listDatasets(ensembl)查看
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# attr为基因id目标类型，listAttributes(ensembl)展示其他id输出类型
attr <- c("ensembl_gene_id", "hgnc_symbol")

# listFilters(ensembl)展示其他id输入类型
ids <- getBM(attributes = attr, 
             filters = "ensembl_gene_id", 
             values = geneid, 
             mart = ensembl)

# 将表达矩阵与ids合并
heatmap_data <- heatmap_data %>% 
  mutate(Names = rownames(heatmap_data))
heatmap_data <- merge(ids, heatmap_data, 
                 by.x = "ensembl_gene_id", 
                 by.y = "Names")
heatmap_data <- heatmap_data[!c(heatmap_data$hgnc_symbol == ""), ]
row.names(heatmap_data) <- heatmap_data[, 2]
heatmap_data <- heatmap_data[, c(-1, -2)]

p2 <- pheatmap(mat = heatmap_data, 
               cellwidth = 6, cellheight = 6, 
               cluster_cols = FALSE, cluster_rows = TRUE, 
               angle_col = 90, main = "log10(rpkm + 1)", 
               color = c(colorRampPalette(colors = c("#083490","#f1f1f1"))(length(bk)/2), 
                         colorRampPalette(colors = c("#f1f1f1","#E60212"))(length(bk)/2)),
               legend_breaks=seq(0, 4, 0.4), 
               breaks = bk, 
               fontsize = 7)
p2

pdf(file = "03-results/18-Development_exp_pattern/All_coding_UCR_related_coding_exp_heatmap_detail.pdf", 
    width = 1920/254, height = 5760/254)
p2
dev.off()

# classify gene to go pathway
GO_id <- getBM(attributes = c('ensembl_gene_id', 'go_id'), 
               filters = "ensembl_gene_id", 
               values = geneid, 
               mart = ensembl)

coding_UCR_related_genes_enriched_GO_id <- read_xlsx(
  path = "03-results/10-GO_enrichment_analysis/meatscape_coding_UCR_genes_GO/metascape_result.xlsx", 
  sheet = 2
)

GO_id <- GO_id[GO_id$go_id %in% coding_UCR_related_genes_enriched_GO_id$Term, ]
GO_id <- merge(GO_id, 
               coding_UCR_related_genes_enriched_GO_id[, c(3:4)], 
               by.x = "go_id", 
               by.y = "Term")
GO_id <- merge(GO_id, ids, 
               by = "ensembl_gene_id")

# add bar label
mRNA_metabolic_process <- c(
  "GO:0016071", "GO:0008380", "GO:0043484", "GO:0006397", 
  "GO:0000377", "GO:0000398", "GO:0000375", "GO:0048024", 
  "GO:0050684", "GO:0000000", "GO:1903311", "GO:0000381"
)

brain_development <- c(
  "GO:0007420", "GO:0060322", "GO:0045664", "GO:0045665", 
  "GO:0045596", "GO:0061564", "GO:0007409", "GO:0048667", 
  "GO:0045666", "GO:0000902", "GO:0048812", "GO:0120039", 
  "GO:0048858", "GO:0001764", "GO:0031175", "GO:0032989", 
  "GO:0007411", "GO:0097485"
)

alternative_mRNA_splicing <- c(
  "GO:0000380", "GO:0033119", "GO:1903312", "GO:0048025", 
  "GO:0050686", "GO:0006376", "GO:0000245", "GO:0022613", 
  "GO:0022618", "GO:0071826", "GO:0045165"
)

mRNA_metabolic_process_go <- GO_id[GO_id$go_id %in% mRNA_metabolic_process, ]
brain_development_go <- GO_id[GO_id$go_id %in% brain_development, ]
alternative_mRNA_splicing_go <- GO_id[GO_id$go_id %in% alternative_mRNA_splicing, ]

# brain development go
brain_development_go <- brain_development_go %>% 
  mutate(Term = "brain_development")
annotation_row_brain_dev <- brain_development_go[, c(4:5)]
annotation_row_brain_dev <- unique(annotation_row_brain_dev)
row.names(annotation_row_brain_dev) <- annotation_row_brain_dev[, 1]
# annotation_row_brain_dev <- (annotation_row_brain_dev[, -1])

# mRNA metabolic process
mRNA_metabolic_process_go <- mRNA_metabolic_process_go %>% 
  mutate(Term = "mRNA_metabolic_process")
annotation_row_mRNA_meta <- mRNA_metabolic_process_go[, c(4:5)]
annotation_row_mRNA_meta <- unique(annotation_row_mRNA_meta)
row.names(annotation_row_mRNA_meta) <- annotation_row_mRNA_meta[, 1]
# annotation_row_mRNA_meta <- (annotation_row_mRNA_meta[, -1])

# AS
alternative_mRNA_splicing_go <- alternative_mRNA_splicing_go %>% 
  mutate(Term = "alternative_mRNA_splicing")
annotation_row_AS <- alternative_mRNA_splicing_go[, c(4:5)]
annotation_row_AS <- unique(annotation_row_AS)
row.names(annotation_row_AS) <- annotation_row_AS[, 1]

# get annotation row
annotation_row <- merge(annotation_row_mRNA_meta, 
                        annotation_row_brain_dev, 
                        by = "hgnc_symbol", 
                        all = TRUE)
annotation_row <- merge(annotation_row, 
                        annotation_row_AS, 
                        by = "hgnc_symbol", 
                        all = TRUE)
row.names(annotation_row) <- annotation_row[, 1]
annotation_row <- annotation_row[, -1]
colnames(annotation_row) <- c("mRNA_metabolic_process", 
                              "brain_development", 
                              "alternative_splicing")

p3 <- pheatmap(mat = heatmap_data, 
               cellwidth = 6, cellheight = 6, 
               cluster_cols = FALSE, cluster_rows = TRUE, 
               annotation_row = annotation_row, annotation_legend = FALSE, 
               angle_col = 90, main = "log10(rpkm + 1)", 
               color = c(colorRampPalette(colors = c("#083490","#f1f1f1"))(length(bk)/2), 
                         colorRampPalette(colors = c("#f1f1f1","#E60212"))(length(bk)/2)),
               legend_breaks=seq(0, 4, 0.4), 
               breaks = bk, 
               fontsize = 7)
p3

pdf(file = "03-results/18-Development_exp_pattern/detailed_heatmap_of_all_coding_UCR_related_genes_containing_enriched_pathway_markers.pdf", 
    width = 1920/254, height = 5760/254)
p3
dev.off()

# 不展示基因名称
p4 <- pheatmap(mat = heatmap_data, 
               cellwidth = 6, cellheight = 2, 
               cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = FALSE, # 不展示行名
               annotation_row = annotation_row, annotation_legend = FALSE, # 不展示注释标签
               angle_col = 90, main = "log10(rpkm + 1)", 
               color = c(colorRampPalette(colors = c("#083490","#f1f1f1"))(length(bk)/2), 
                         colorRampPalette(colors = c("#f1f1f1","#E60212"))(length(bk)/2)),
               legend_breaks=seq(0, 4, 0.4), 
               breaks = bk, 
               fontsize = 7)
p4

pdf(file = "03-results/18-Development_exp_pattern/heatmap_of_all_coding_UCR_related_genes_containing_enriched_pathway_markers.pdf", 
    width = 1700/254, height = 1700/254)
p4
dev.off()














