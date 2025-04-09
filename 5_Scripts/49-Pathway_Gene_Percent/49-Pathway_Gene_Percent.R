# get the proportion of genes overlap with type I UCRs from different pathways

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library(readxl)

GOResults <- read_xlsx(
  path = "03-results/10-GO_enrichment_analysis/meatscape_coding_UCR_genes_GO/metascape_result.xlsx", 
  sheet = 2)
GOResults <- GOResults[str_detect(GOResults$GroupID, "Summary"), ]
























