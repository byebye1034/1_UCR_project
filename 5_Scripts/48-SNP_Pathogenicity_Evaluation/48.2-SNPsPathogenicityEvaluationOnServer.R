# ClinPred Score

setwd(dir = "/home/wangxiangting/baiyun0216/projects/UCR/data/41-SNPsPathogenicityEvaluation")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library(stringr)

## load all SNPs info
# load("passed_UCR_SNPs.Rdata")
# load("passed_UCR_left_SNPs.Rdata")
# load("passed_UCR_right_SNPs.Rdata")
# load("passed_rf_SNPs.Rdata")

## load ClinPred
load("ClinPred_hg38.Rdata")

## import UCRs SNPs
snpsUCRs <- read.delim(
  file = "snpsUCRsVCF4FATH.vcf", 
  header = TRUE
)

snpsUCRsClinPred <- merge(snpsUCRs[, c(1:5)], 
                          ClinPred, 
                          by.x = c("X.CHROM", "POS", "REF", "ALT"), 
                          by.y = c("Chr", "Start", "Ref", "Alt"))

write.table(snpsUCRsClinPred, 
            file = "snpsUCRsClinPred.txt", 
            sep = "\t", col.names = TRUE, 
            row.names = FALSE, quote = FALSE)

snpsUCRsClinPred <- read.table(
  file = "01-data/48-SNP_Pathogenicity_Evaluation/snpsUCRsClinPred.txt", 
  sep = "\t", header = TRUE
)





























