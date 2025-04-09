# Pathogenicity evaluation of SNPs within UCRs and control fragments
# organize SNPs within UCRs and control regions to VCF files
# obtain all pathogenicity score respectively

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(tidyr)
library(data.table)
library(stringr)

# oragnize SNPs within UCRs and control regions to VCF files --------------

## load all SNPs info
# UCR
load("D:/R_project/UCR_project/02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/passed_UCR_SNPs.Rdata")
write.table(passed_UCR_SNPs$refsnp_id, 
            file = "01-data/48-SNP_Pathogenicity_Evaluation/snpPassedWithinUCRs.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE) # for VEP

# oragnize to VCF format
snpsUCRsVCF <- passed_UCR_SNPs %>% 
  separate(allele, into = c("REF", "ALT"), sep = "/", extra = "merge", fill = "right") %>% 
  mutate(ALT = str_replace_all(ALT, "/", ","))
snpsUCRsVCF <- snpsUCRsVCF[, c(5, 8, 2, 6, 7)]
snpsUCRsVCF <- unique(snpsUCRsVCF)
colnames(snpsUCRsVCF) <- c("#CHROM", "POS", "ID", "REF", "ALT")
snpsUCRsVCF <- snpsUCRsVCF %>% 
  mutate(QUAL = ".") %>% 
  mutate(FITTER = ".") %>% 
  mutate(INFO = ".") %>% 
  mutate(POS = as.numeric(POS)) %>% 
  arrange(`#CHROM`, POS)

# set chunk size: 2000 rows
chunkSize <- 2000

# calculate num of splits needed
nChunk <- ceiling(nrow(snpsUCRsVCF) / chunkSize)

# save all splited files
for (i in 1:nChunk) {
  # calculate start row and end row for every chunk
  startRow <- (i - 1)*chunkSize + 1
  endRow <- min(i*2000, nrow(snpsUCRsVCF))
  
  chunk <- snpsUCRsVCF[startRow:endRow, ]
  
  # save to VCF files
  write.table(
    chunk, 
    file = paste0("02-analysis/48-SNP_Pathogenicity_Evaluation/SIFT4G/SNPs/snpsUCRsVCF-", startRow, "-", endRow, ".vcf"), 
    sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE
  )
}

# set chunk size: 1000 rows for CADD
chunkSize <- 1000

# calculate num of splits needed
nChunk <- ceiling(nrow(snpsUCRsVCF) / chunkSize)

# save all splited files
for (i in 1:nChunk) {
  # calculate start row and end row for every chunk
  startRow <- (i - 1)*chunkSize + 1
  endRow <- min(i*1000, nrow(snpsUCRsVCF))
  
  chunk <- snpsUCRsVCF[startRow:endRow, ]
  
  # save to VCF files
  write.table(
    chunk, 
    file = paste0("02-analysis/48-SNP_Pathogenicity_Evaluation/CADD/SNPs/snpsUCRsVCF-", startRow, "-", endRow, ".vcf"), 
    sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE
  )
}

# FATHMM-XF: one variant in every line
snpsUCRsVCF4FATH <- snpsUCRsVCF %>% 
  separate_rows(5, sep = ",") %>% 
  arrange(`#CHROM`, POS) # 24269

write.table(
  snpsUCRsVCF4FATH, 
  file = "02-analysis/48-SNP_Pathogenicity_Evaluation/FATHMM-XF/SNPs/snpsUCRsVCF4FATH.vcf", 
  sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE
)

chunkSize <- 2000
nChunk <- ceiling(nrow(snpsUCRsVCF4FATH) / chunkSize)

# save all splited files
for (i in 1:nChunk) {
  # calculate start row and end row for every chunk
  startRow <- (i - 1)*chunkSize + 1
  endRow <- min(i*2000, nrow(snpsUCRsVCF4FATH))
  
  chunk <- snpsUCRsVCF4FATH[startRow:endRow, ]
  
  # save to VCF files
  write.table(
    chunk, 
    file = paste0("02-analysis/48-SNP_Pathogenicity_Evaluation/FATHMM-XF/SNPs/snpsUCRsVCF-", startRow, "-", endRow, ".vcf"), 
    sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE
  )
}


# UCR left
load("D:/R_project/UCR_project/02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/passed_UCR_left_SNPs.Rdata")
write.table(passed_UCR_left_SNPs$refsnp_id, 
            file = "01-data/48-SNP_Pathogenicity_Evaluation/snpPassedWithinUCRsleft.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

snpsUCRsLeftVCF <- passed_UCR_left_SNPs %>% 
  separate(allele, into = c("REF", "ALT"), sep = "/", extra = "merge", fill = "right") %>% 
  mutate(ALT = str_replace_all(ALT, "/", ","))
snpsUCRsLeftVCF <- snpsUCRsLeftVCF[, c(5, 8, 2, 6, 7)]
snpsUCRsLeftVCF <- unique(snpsUCRsLeftVCF)
colnames(snpsUCRsLeftVCF) <- c("#CHROM", "POS", "ID", "REF", "ALT")
snpsUCRsLeftVCF <- snpsUCRsLeftVCF %>% 
  mutate(QUAL = ".") %>% 
  mutate(FITTER = ".") %>% 
  mutate(INFO = ".") %>% 
  mutate(POS = as.numeric(POS)) %>% 
  arrange(`#CHROM`, POS)

# all
write.table(snpsUCRsLeftVCF, 
            file = "02-analysis/48-SNP_Pathogenicity_Evaluation/SNPsUCRsLeft/snpsUCRsLeftVCF-1.vcf", 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# set chunk size: 1000 rows
chunkSize <- 1000

# calculate num of splits needed
nChunk <- ceiling(nrow(snpsUCRsLeftVCF) / chunkSize)

# save all splited files
for (i in 1:nChunk) {
  # calculate start row and end row for every chunk
  startRow <- (i - 1)*chunkSize + 1
  endRow <- min(i*1000, nrow(snpsUCRsLeftVCF))
  
  chunk <- snpsUCRsLeftVCF[startRow:endRow, ]
  
  # save to VCF files
  write.table(
    chunk, 
    file = paste0("02-analysis/48-SNP_Pathogenicity_Evaluation/SNPsUCRsLeft/snpsUCRsLeftVCF-1-", startRow, "-", endRow, ".vcf"), 
    sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE
  )
}

# set chunk size: 2000 rows
chunkSize <- 2000

# calculate num of splits needed
nChunk <- ceiling(nrow(snpsUCRsLeftVCF) / chunkSize)

# save all splited files
for (i in 1:nChunk) {
  # calculate start row and end row for every chunk
  startRow <- (i - 1)*chunkSize + 1
  endRow <- min(i*2000, nrow(snpsUCRsLeftVCF))
  
  chunk <- snpsUCRsLeftVCF[startRow:endRow, ]
  
  # save to VCF files
  write.table(
    chunk, 
    file = paste0("02-analysis/48-SNP_Pathogenicity_Evaluation/SNPsUCRsLeft/snpsUCRsLeftVCF-1-", startRow, "-", endRow, ".vcf"), 
    sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE
  )
}

# type 2
snpsUCRsLeftVCF2 <- snpsUCRsLeftVCF %>% 
  separate_rows(5, sep = ",") %>% 
  arrange(`#CHROM`, POS)

write.table(
  snpsUCRsLeftVCF2, 
  file = "02-analysis/48-SNP_Pathogenicity_Evaluation/SNPsUCRsLeft/snpsUCRsLeftVCF-2.vcf", 
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)

# set chunk size: 1000 rows
chunkSize <- 1000

# calculate num of splits needed
nChunk <- ceiling(nrow(snpsUCRsLeftVCF2) / chunkSize)

# save all splited files
for (i in 1:nChunk) {
  # calculate start row and end row for every chunk
  startRow <- (i - 1)*chunkSize + 1
  endRow <- min(i*1000, nrow(snpsUCRsLeftVCF2))
  
  chunk <- snpsUCRsLeftVCF2[startRow:endRow, ]
  
  # save to VCF files
  write.table(
    chunk, 
    file = paste0("02-analysis/48-SNP_Pathogenicity_Evaluation/SNPsUCRsLeft/snpsUCRsLeftVCF-2-", startRow, "-", endRow, ".vcf"), 
    sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE
  )
}

# set chunk size: 2000 rows
chunkSize <- 2000

# calculate num of splits needed
nChunk <- ceiling(nrow(snpsUCRsLeftVCF2) / chunkSize)

# save all splited files
for (i in 1:nChunk) {
  # calculate start row and end row for every chunk
  startRow <- (i - 1)*chunkSize + 1
  endRow <- min(i*2000, nrow(snpsUCRsLeftVCF2))
  
  chunk <- snpsUCRsLeftVCF2[startRow:endRow, ]
  
  # save to VCF files
  write.table(
    chunk, 
    file = paste0("02-analysis/48-SNP_Pathogenicity_Evaluation/SNPsUCRsLeft/snpsUCRsLeftVCF-2-", startRow, "-", endRow, ".vcf"), 
    sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE
  )
}



# UCR right
load("D:/R_project/UCR_project/02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/passed_UCR_right_SNPs.Rdata")
write.table(passed_UCR_right_SNPs$refsnp_id, 
            file = "01-data/48-SNP_Pathogenicity_Evaluation/snpPassedWithinUCRsright.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

snpsUCRsRightVCF <- passed_UCR_right_SNPs %>% 
  separate(allele, into = c("REF", "ALT"), sep = "/", extra = "merge", fill = "right") %>% 
  mutate(ALT = str_replace_all(ALT, "/", ","))
snpsUCRsRightVCF <- snpsUCRsRightVCF[, c(5, 8, 2, 6, 7)]
snpsUCRsRightVCF <- unique(snpsUCRsRightVCF)
colnames(snpsUCRsRightVCF) <- c("#CHROM", "POS", "ID", "REF", "ALT")
snpsUCRsRightVCF <- snpsUCRsRightVCF %>% 
  mutate(QUAL = ".") %>% 
  mutate(FITTER = ".") %>% 
  mutate(INFO = ".") %>% 
  mutate(POS = as.numeric(POS)) %>% 
  arrange(`#CHROM`, POS)

# all
write.table(snpsUCRsRightVCF, 
            file = "02-analysis/48-SNP_Pathogenicity_Evaluation/SNPsUCRsRight/snpsUCRsRightVCF-1.vcf", 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# set chunk size: 1000 rows
chunkSize <- 1000

# calculate num of splits needed
nChunk <- ceiling(nrow(snpsUCRsRightVCF) / chunkSize)

# save all splited files
for (i in 1:nChunk) {
  # calculate start row and end row for every chunk
  startRow <- (i - 1)*chunkSize + 1
  endRow <- min(i*1000, nrow(snpsUCRsRightVCF))
  
  chunk <- snpsUCRsRightVCF[startRow:endRow, ]
  
  # save to VCF files
  write.table(
    chunk, 
    file = paste0("02-analysis/48-SNP_Pathogenicity_Evaluation/SNPsUCRsRight/snpsUCRsRightVCF-1-", startRow, "-", endRow, ".vcf"), 
    sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE
  )
}

# set chunk size: 2000 rows
chunkSize <- 2000

# calculate num of splits needed
nChunk <- ceiling(nrow(snpsUCRsRightVCF) / chunkSize)

# save all splited files
for (i in 1:nChunk) {
  # calculate start row and end row for every chunk
  startRow <- (i - 1)*chunkSize + 1
  endRow <- min(i*2000, nrow(snpsUCRsRightVCF))
  
  chunk <- snpsUCRsRightVCF[startRow:endRow, ]
  
  # save to VCF files
  write.table(
    chunk, 
    file = paste0("02-analysis/48-SNP_Pathogenicity_Evaluation/SNPsUCRsRight/snpsUCRsRightVCF-1-", startRow, "-", endRow, ".vcf"), 
    sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE
  )
}

# type 2
snpsUCRsRightVCF2 <- snpsUCRsRightVCF %>% 
  separate_rows(5, sep = ",") %>% 
  arrange(`#CHROM`, POS)

write.table(
  snpsUCRsRightVCF2, 
  file = "02-analysis/48-SNP_Pathogenicity_Evaluation/SNPsUCRsRight/snpsUCRsRightVCF-2.vcf", 
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)

# set chunk size: 1000 rows
chunkSize <- 1000

# calculate num of splits needed
nChunk <- ceiling(nrow(snpsUCRsRightVCF2) / chunkSize)

# save all splited files
for (i in 1:nChunk) {
  # calculate start row and end row for every chunk
  startRow <- (i - 1)*chunkSize + 1
  endRow <- min(i*1000, nrow(snpsUCRsRightVCF2))
  
  chunk <- snpsUCRsRightVCF2[startRow:endRow, ]
  
  # save to VCF files
  write.table(
    chunk, 
    file = paste0("02-analysis/48-SNP_Pathogenicity_Evaluation/SNPsUCRsRight/snpsUCRsRightVCF-2-", startRow, "-", endRow, ".vcf"), 
    sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE
  )
}

# set chunk size: 2000 rows
chunkSize <- 2000

# calculate num of splits needed
nChunk <- ceiling(nrow(snpsUCRsRightVCF2) / chunkSize)

# save all splited files
for (i in 1:nChunk) {
  # calculate start row and end row for every chunk
  startRow <- (i - 1)*chunkSize + 1
  endRow <- min(i*2000, nrow(snpsUCRsRightVCF2))
  
  chunk <- snpsUCRsRightVCF2[startRow:endRow, ]
  
  # save to VCF files
  write.table(
    chunk, 
    file = paste0("02-analysis/48-SNP_Pathogenicity_Evaluation/SNPsUCRsRight/snpsUCRsRightVCF-2-", startRow, "-", endRow, ".vcf"), 
    sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE
  )
}



# Random Fragments
load("D:/R_project/UCR_project/02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/passed_rf_SNPs.Rdata")
write.table(passed_rf_SNPs$refsnp_id, 
            file = "01-data/48-SNP_Pathogenicity_Evaluation/snpPassedWithinRandomFragments.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

snpsRandfVCF <- passed_rf_SNPs %>% 
  separate(col = allele, into = c("REF", "ALT"), sep = "/", extra = "merge", fill = "right") %>% 
  mutate(ALT = str_replace_all(ALT, "/", ","))
snpsRandfVCF <- snpsRandfVCF[, c(5, 8, 2, 6, 7)]
snpsRandfVCF <- unique(snpsRandfVCF)
colnames(snpsRandfVCF) <- c("#CHROM", "POS", "ID", "REF", "ALT")
snpsRandfVCF <- snpsRandfVCF %>% 
  mutate(QUAL = ".") %>% 
  mutate(FITTER = ".") %>% 
  mutate(INFO = ".") %>% 
  mutate(POS = as.numeric(POS)) %>% 
  arrange(`#CHROM`, POS)

# all
write.table(snpsRandfVCF, 
            file = "02-analysis/48-SNP_Pathogenicity_Evaluation/SNPsRandf/snpsRandfVCF-1.vcf", 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# set chunk size: 1000 rows
chunkSize <- 1000

# calculate num of splits needed
nChunk <- ceiling(nrow(snpsRandfVCF) / chunkSize)

# save all splited files
for (i in 1:nChunk) {
  # calculate start row and end row for every chunk
  startRow <- (i - 1)*chunkSize + 1
  endRow <- min(i*1000, nrow(snpsRandfVCF))
  
  chunk <- snpsRandfVCF[startRow:endRow, ]
  
  # save to VCF files
  write.table(
    chunk, 
    file = paste0("02-analysis/48-SNP_Pathogenicity_Evaluation/SNPsRandf/snpsRandfVCF-1-", startRow, "-", endRow, ".vcf"), 
    sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE
  )
}

# set chunk size: 2000 rows
chunkSize <- 2000

# calculate num of splits needed
nChunk <- ceiling(nrow(snpsRandfVCF) / chunkSize)

# save all splited files
for (i in 1:nChunk) {
  # calculate start row and end row for every chunk
  startRow <- (i - 1)*chunkSize + 1
  endRow <- min(i*2000, nrow(snpsRandfVCF))
  
  chunk <- snpsRandfVCF[startRow:endRow, ]
  
  # save to VCF files
  write.table(
    chunk, 
    file = paste0("02-analysis/48-SNP_Pathogenicity_Evaluation/SNPsRandf/snpsRandfVCF-1-", startRow, "-", endRow, ".vcf"), 
    sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE
  )
}

# type 2
snpsRandfVCF2 <- snpsRandfVCF %>% 
  separate_rows(5, sep = ",") %>% 
  arrange(`#CHROM`, POS)

write.table(
  snpsRandfVCF2, 
  file = "02-analysis/48-SNP_Pathogenicity_Evaluation/SNPsRandf/snpsRandfVCF-2.vcf", 
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)

# set chunk size: 1000 rows
chunkSize <- 1000

# calculate num of splits needed
nChunk <- ceiling(nrow(snpsRandfVCF2) / chunkSize)

# save all splited files
for (i in 1:nChunk) {
  # calculate start row and end row for every chunk
  startRow <- (i - 1)*chunkSize + 1
  endRow <- min(i*1000, nrow(snpsRandfVCF2))
  
  chunk <- snpsRandfVCF2[startRow:endRow, ]
  
  # save to VCF files
  write.table(
    chunk, 
    file = paste0("02-analysis/48-SNP_Pathogenicity_Evaluation/SNPsRandf/snpsRandfVCF-2-", startRow, "-", endRow, ".vcf"), 
    sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE
  )
}

# set chunk size: 2000 rows
chunkSize <- 2000

# calculate num of splits needed
nChunk <- ceiling(nrow(snpsRandfVCF2) / chunkSize)

# save all splited files
for (i in 1:nChunk) {
  # calculate start row and end row for every chunk
  startRow <- (i - 1)*chunkSize + 1
  endRow <- min(i*2000, nrow(snpsRandfVCF2))
  
  chunk <- snpsRandfVCF2[startRow:endRow, ]
  
  # save to VCF files
  write.table(
    chunk, 
    file = paste0("02-analysis/48-SNP_Pathogenicity_Evaluation/SNPsRandf/snpsRandfVCF-2-", startRow, "-", endRow, ".vcf"), 
    sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE
  )
}




# FATHMM-XF ---------------------------------------------------------------

UCRsFATHMMXF <- read.table(
  file = "02-analysis/48-SNP_Pathogenicity_Evaluation/FATHMM-XF/FATHMMXFScores/UCRs-FATHMMXF.txt", 
  sep = "\t"
)
colnames(UCRsFATHMMXF) <- c("#CHROM", "POS", "REF", "ALT", 
                            "CodingScore", "NonCodingScore", "Warning")
nrow(UCRsFATHMMXF[UCRsFATHMMXF$CodingScore == "--", ]) # 21737
nrow(UCRsFATHMMXF[UCRsFATHMMXF$NonCodingScore == "--", ]) # 4543
nrow(UCRsFATHMMXF[UCRsFATHMMXF$CodingScore != "--" & UCRsFATHMMXF$NonCodingScore != "--", ])

UCRsFATHMMXF <- UCRsFATHMMXF %>% 
  mutate(FATHMMXF = ifelse(CodingScore == "--", NonCodingScore, CodingScore))
nrow(UCRsFATHMMXF[UCRsFATHMMXF$FATHMMXF == "--", ]) # no codingScore and noncodingScore: 2011

# ClinPred ----------------------------------------------------------------

ClinPred <- fread(input = "01-data/48-SNP_Pathogenicity_Evaluation/ClinPred_hg38.txt", 
                  sep = "\t", header = TRUE)
ClinPred <- as.tibble(ClinPred)
head(ClinPred)
save(ClinPred, file = "01-data/48-SNP_Pathogenicity_Evaluation/ClinPred_hg38.Rdata")

#       Chr  Start    Ref    Alt     ClinPred_Score
#    <char> <char> <char> <char>             <char>
# 1:      1  69091      A      C 0.0722438414657657
# 2:      1  69091      A      G 0.0408461190417061
# 3:      1  69091      A      T 0.0598572379421533
# 4:      1  69092      T      A   0.38129135966301
# 5:      1  69092      T      C  0.121829316020012
# 6:      1  69092      T      G  0.373158395290375

head(passed_UCR_SNPs)
subset(ClinPred, ClinPred$Start == 10537645)
subset(ClinPred, ClinPred$Start == 69091)

ClinPred <- ClinPred %>% 
  mutate(Allele = str_c(Ref, Alt, sep = "/"))
head(ClinPred)

passed_UCR_SNPs_1000 <- head(passed_UCR_SNPs, 1000)
ClinPred_1000 <- head(ClinPred, 1000)

UCRsClinPred <- passed_UCR_SNPs_1000 %>% 
  select(chr_name, allele, chrom_start) %>% 
  left_join(ClinPred_1000 %>% 
              select(Chr, Start, Allele), 
            by = c("chr_name" = "Chr", "allele" = "Allele", "chrom_start" = "Start"))

UCRsClinPred <- ClinPred_1000 %>% 
  merge(passed)





















