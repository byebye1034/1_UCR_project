# organize the classification of UCR

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library(ggplot2)

codingUCR <- read.table(
  file = "02-analysis/16-New_classification/01-UCR_from_protein_coding_gene.bed", 
  sep = "\t"
)
length(unique(codingUCR$V4)) # 314
codingUCR <- codingUCR %>% 
  mutate(type = "codingUCR")

ncRNAUCR <- read.table(
  file = "02-analysis/16-New_classification/02-UCR_from_Non_Coding_RNA(66).bed", 
  sep = "\t"
)
length(unique(ncRNAUCR$V4)) # 66
ncRNAUCR <- ncRNAUCR %>% 
  mutate(type = "ncRNAUCR")

otherUCR <- read.table(
  file = "02-analysis/16-New_classification/03-UCR_intergenic.bed", 
  sep = "\t"
)
length(unique(otherUCR$V4)) # 99
otherUCR <- otherUCR %>% 
  mutate(type = "otherUCR")

UCRType <- rbind(unique(codingUCR[, c(1:4, 11)]), unique(ncRNAUCR[, c(1:4, 11)]))
UCRType <- rbind(UCRType, otherUCR)
UCRType <- UCRType %>% 
  arrange(V1, V2, V3, V4)
save(UCRType, file = "02-analysis/16-New_classification/UCRType.rdata")












