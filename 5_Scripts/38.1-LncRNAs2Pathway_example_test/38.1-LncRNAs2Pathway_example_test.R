# 用LncRNAsPathway文章里给的数据测试一下我使用的方法对不对

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(LncPath)
library(stringr)

diff_exp_lnc <- read.table(file = "01-data/38.1-LncRNAs2Pathway_example_test/Table_S3_DELncRNA_prostate_cancer.txt", 
                           sep = "\t", header = FALSE)

# get lncRNA-mRNA interaction network
NetLncPath <- getNet()
print(head(NetLncPath), row.names = FALSE)
length(unique(NetLncPath$V1))
length(unique(NetLncPath$V2))

# get lncRNA set
SigLncs <- diff_exp_lnc$V1

Result <- lncPath(LncRNAList = SigLncs, 
                  Network = NetLncPath, 
                  Weighted = TRUE, 
                  PathwayDataSet = "KEGG", 
                  nperm = 1000, 
                  minPathSize = 15, 
                  maxPathSize = 500)

save(Result, file = "02-analysis/38.1-LncRNAs2Pathway_example_test/Result.Rdata")

# Generate a table of the summary of each pathway
PathwaySummaryTable <- lncPath2Table(Result)
print(head(PathwaySummaryTable), row.names = FALSE)

PathwaySummaryTable <- filter(PathwaySummaryTable, PathwaySummaryTable$`False Discovery Rate` < 0.01)

PathwaySummaryTable$`Gene Set Name` <- gsub("KEGG_", "", PathwaySummaryTable$`Gene Set Name`)
PathwaySummaryTable$`Gene Set Name` <- str_to_lower(PathwaySummaryTable$`Gene Set Name`)

print(head(PathwaySummaryTable, 15), row.names = FALSE)
























