# 从devCE里筛选和大脑发育相关的基因（用devCE所有的基因做GO，挑选神经发育pathway里的），检测其CE是否随发育PSI下降
# 用interproscan检测CE影响的domain
# 筛选devCE里包含GGA的

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(readxl)
library(biomaRt)
library(curl)

# get human brain down devCE ----------------------------------------------

human_devAS <- read.table(file = "01-data/19-Development_alternative_splicing/human.devAS", 
                          sep = ",", header = TRUE)
human_brain_down_devCE <- human_devAS[, c(1:9, 16, 23)] %>% 
  filter(pattern.brain == "d") %>% 
  filter(as.type == "CE")
# 1412个human brain down devCE

length(unique(human_brain_down_devCE$ens_id))
# 981个基因

# 先用GO富集分析找到和神经发育有关的基因，再筛选包含GGA的基因
write.table(unique(human_brain_down_devCE$ens_id), 
            file = "02-analysis/40-Choose_experimentally_validated_devCE/human_brain_down_devCE_ensembl_id.bed", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# 读入GO富集分析结果里神经发育相关的基因
metascape_result <- read_xlsx(path = "02-analysis/40-Choose_experimentally_validated_devCE/all.tq0uubhcs/metascape_result.xlsx", 
                              sheet = 3)
neu_dev_gene <- metascape_result$Symbols[1]
neu_dev_gene <- str_split_1(neu_dev_gene, ",")
print(neu_dev_gene)

# 利用useMart链接到人类的数据库，其他数据库使用listDatasets(ensembl)查看
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# attr为基因id目标类型，listAttributes(ensembl)展示其他id输出类型
attr <- c("ensembl_gene_id", "hgnc_symbol")

# listFilters(ensembl)展示其他id输入类型
ids <- getBM(attributes = attr, 
             filters = "hgnc_symbol", 
             values = neu_dev_gene, 
             mart = ensembl)

neu_dev_gene_id <- ids

human_brain_down_devCE <- filter(human_brain_down_devCE, 
                                 human_brain_down_devCE$ens_id %in% neu_dev_gene_id$ensembl_gene_id)

























