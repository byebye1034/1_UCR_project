# https://home.kaessmannlab.org/resources#apps
# 根据数据框的数据做一些发育相关分析

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(dbplyr)
library(readxl)

DDG_lncRNA <- read_xlsx(path = "01-data/18-Development/development_dynamic_lncRNA.xlsx", 
                        sheet = 1)
DDG_lncRNA <- DDG_lncRNA[, c(3, 5, 6)]

ncRNA_UCR_ensembl_id <- read.table(file = "02-analysis/16-New_classification/ncRNA_UCR_ensembl_id.txt", 
                                   sep = "\t")

ncRNA_UCR_dynamic <- merge(ncRNA_UCR_ensembl_id, DDG_lncRNA, by.x = "V1", by.y = "ENSEMBL75 id")

























