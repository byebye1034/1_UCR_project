# lncRNA的富集分析显示和gland发育，氧气应答等有关
# 在各组织器官看看这些lncRNA的表达pattern，推测其功能

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(biomaRt)
library(curl)

ncRNA_UCR_ensembl_id <- read.table(file = "02-analysis/16-New_classification/ncRNA_UCR_ensembl_id.txt")
ncRNA_UCR_ensembl_id <- unique(ncRNA_UCR_ensembl_id)

# 利用useMart链接到人类的数据库，其他数据库使用listDatasets(ensembl)查看
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# attr为基因id目标类型，listAttributes(ensembl)展示其他id输出类型
attr <- c("ensembl_gene_id", "hgnc_symbol")

# listFilters(ensembl)展示其他id输入类型
ids <- getBM(attributes = attr, 
             filters = "ensembl_gene_id", 
             values = ncRNA_UCR_ensembl_id$V1, 
             mart = ensembl)

ids <- ids[!(ids$hgnc_symbol == ""), ]

write.table(ids$hgnc_symbol, 
            file = "02-analysis/24-Development_lncRNA_exp_pattern/ncRNA_UCR_hgnc_symbol.bed", 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

RNAenrich_converted_lncRNA_id <- 
  read.table(file = "01-data/24-Development_lncRNA_exp_pattern/RNAenrich_converted_lncRNA_id.csv", 
             sep = ",", header = TRUE)

































