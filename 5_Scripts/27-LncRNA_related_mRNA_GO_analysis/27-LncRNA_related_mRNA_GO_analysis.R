# 用lncRNA相关的蛋白编码基因做富集分析

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(biomaRt)

ncRNA_UCR <- read.table(file = "02-analysis/16-New_classification/ncRNA_UCR(66)_ensembl_id.txt", sep = "\t")
ncRNA_UCR <- unique(ncRNA_UCR)

# 利用useMart链接到人类的数据库，其他数据库使用listDatasets(ensembl)查看
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# attr为基因id目标类型，listAttributes(ensembl)展示其他id输出类型
attr <- c("ensembl_gene_id", "hgnc_symbol")

# listFilters(ensembl)展示其他id输入类型
ids <- getBM(attributes = attr, 
             filters = "ensembl_gene_id", 
             values = ncRNA_UCR$V1, 
             mart = ensembl)

# 将表达矩阵与ids合并
write.table(ids$hgnc_symbol, file = "02-analysis/27-LncRNA_related_mRNA_GO_analysis/ncRNA_UCR_symbol.txt", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)







