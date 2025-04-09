# How many statistically coding UCR-related genes are associated with AS

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library()

UCRs_from_PCGs <- read.table(
  file = "02-analysis/16-New_classification/01-UCR_from_protein_coding_gene.bed")
length(unique(UCRs_from_PCGs$V4))
length(unique(UCRs_from_PCGs$V9))





