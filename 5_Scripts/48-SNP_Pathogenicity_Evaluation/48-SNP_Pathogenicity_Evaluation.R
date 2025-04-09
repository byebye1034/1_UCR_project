# Pathogenicity evaluation of SNPs within UCRs and control fragments

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(stringr)

# obtainAllSNPsWithinUCRs -------------------------------------------------

# UCR
load("D:/R_project/UCR_project/02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/passed_UCR_SNPs.Rdata")
write.table(passed_UCR_SNPs$refsnp_id, 
            file = "01-data/48-SNP_Pathogenicity_Evaluation/snpPassedWithinUCRs.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# UCR left
load("D:/R_project/UCR_project/02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/passed_UCR_left_SNPs.Rdata")
write.table(passed_UCR_left_SNPs$refsnp_id, 
            file = "01-data/48-SNP_Pathogenicity_Evaluation/snpPassedWithinUCRsleft.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# UCR right
load("D:/R_project/UCR_project/02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/passed_UCR_right_SNPs.Rdata")
write.table(passed_UCR_right_SNPs$refsnp_id, 
            file = "01-data/48-SNP_Pathogenicity_Evaluation/snpPassedWithinUCRsright.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Random Fragments
load("D:/R_project/UCR_project/02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/passed_rf_SNPs.Rdata")
write.table(passed_rf_SNPs$refsnp_id, 
            file = "01-data/48-SNP_Pathogenicity_Evaluation/snpPassedWithinRandomFragments.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# organizeSNPPathogenicityScoresObtainedFromVEP ---------------------------

SNPPathogenicityScores <- read.delim(
  file = "01-data/48-SNP_Pathogenicity_Evaluation/vepUCRsResults.txt", 
  sep = "\t", header = TRUE, 
)
SNPPathogenicityScores <- SNPPathogenicityScores[c(5:299691), ]
colnames(SNPPathogenicityScores)
headerSNPPathogenicityScores <- head(SNPPathogenicityScores, 15)

# get Scores
SNPPathogenicityScores <- SNPPathogenicityScores[, c("X.Uploaded_variation", 
                                                     "Existing_variation", 
                                                     "REF_ALLELE", "UPLOADED_ALLELE", 
                                                     "SIFT", "PolyPhen", "AF", 
                                                     "CLIN_SIG", "SOMATIC", 
                                                     "PHENO", "PUBMED", "VAR_SYNONYMS", 
                                                     "NMD", "MaxEntScan_alt", 
                                                     "MaxEntScan_diff", "MaxEntScan_ref", 
                                                     "DisGeNET", "PHENOTYPES", 
                                                     "ada_score", "rf_score", 
                                                     "am_class", "am_pathogenicity")]
headerSNPPathogenicityScores <- head(SNPPathogenicityScores, 15)

length(unique(SNPPathogenicityScores$Existing_variation)) # 21708

# SIFT --------------------------------------------------------------------

SIFTScores <- SNPPathogenicityScores[, c(2, 5)]
SIFTScores <- unique(SIFTScores) # total: 24714
table(SIFTScores$SIFT) # no SIFTScores: 21708

SIFTScores <- subset(SIFTScores, !(SIFTScores$SIFT == "-"))
SIFTScores$SIFTScore <- str_extract(SIFTScores$SIFT, 
                                    "\\((.*?)\\)")
SIFTScores$SIFTScore <- gsub("\\(", "", SIFTScores$SIFTScore)
SIFTScores$SIFTScore <- gsub("\\)", "", SIFTScores$SIFTScore)
SIFTScores$SIFTScore <- as.numeric(SIFTScores$SIFTScore)
SIFTScores <- unique(SIFTScores)

lwd_pt <- .pt*72.27/96

myTheme <- theme(
panel.grid = element_blank(), 
panel.background = element_blank(), 
text = element_text(size = 7, color = "#000000"), 
line = element_line(linewidth = 0.5/lwd_pt, color = "#000000"), 

axis.title = element_text(size = 7, color = "#000000"), 
axis.text = element_text(size = 7, color = "#000000"), 
axis.ticks = element_line(linewidth = 0.5/lwd_pt, color = "#000000"), 
axis.line = element_line(linewidth = 0.5/lwd_pt, color = "#000000"),

legend.position = "top", 
legend.key.size = unit(0.25, "cm"), 
legend.title = element_text(size = 7, color = "#000000"), 
legend.text = element_text(size = 7, color = "#000000"), 

plot.title = element_text(size = 7, color = "#000000"),

aspect.ratio = 1:1
)

plotSIFTDistribution <- ggplot(
  data = SIFTScores, 
  mapping = aes(x = SIFTScore)
) +
  geom_density(color = "#E64B35FF", 
               adjust = 1.5) +
  
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  
  myTheme
  
plotSIFTDistribution

# NMD ---------------------------------------------------------------------

# check the association of coding UCR genes and NMD
# obtain all UCRs and type info 
load("D:/R_project/UCR_project/02-analysis/16-New_classification/UCRType.rdata")

# obtain all SNPs info
load("D:/R_project/UCR_project/02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/passed_UCR_SNPs.Rdata")

UCRSNPsINFO <- merge(UCRType[, c(4:5)], passed_UCR_SNPs[,c(1:2)], 
                     by.x = "V4", by.y = "list_type")
unique(UCRSNPsINFO)

UCRSNPsNMD <- merge(UCRSNPsINFO, SNPPathogenicityScores[, c(1, 46)], 
                    by.x = "refsnp_id", by.y = "X.Uploaded_variation")
UCRSNPsNMD <- unique(UCRSNPsNMD)
table(UCRSNPsNMD$NMD)

# Ada_score ---------------------------------------------------------------

AdaScore <- SNPPathogenicityScores[, c(2:4, 19)]
AdaScore <- unique(AdaScore) # 22052
length(unique(AdaScore$Existing_variation)) # 21708
table(AdaScore$ada_score) # not avail

notAvail <- subset(AdaScore, AdaScore$ada_score == "-")
length(unique(notAvail$Existing_variation)) # 21484

# 21484/21708 is not avail

# RF_score ----------------------------------------------------------------

rfScore <- SNPPathogenicityScores[, c(2:4, 20)]
rfScore <- unique(rfScore)
length(unique(rfScore$Existing_variation)) # 21708

rfScoreNA <- subset(rfScore, rfScore$rf_score == "-")
length(unique(rfScoreNA$Existing_variation)) # 21488

# "-" to "0" and character to numeric
rfScore <- rfScore %>% 
  mutate(rf_score = ifelse(rf_score == "-", "0", rf_score)) %>% 
  mutate(rf_score = as.numeric(rf_score))

summary(rfScore$rf_score)

# the num of SNPs predicted to be pathogenic by rf score: 75
sum(rfScore$rf_score > 0.6, na.rm = TRUE)

plotRfScoreDistribution <- 
  ggplot(data = rfScore, mapping = aes(x = rf_score)) +
  geom_histogram()

plotRfScoreDistribution

# CADD --------------------------------------------------------------------

CADD <- read.table(file = "D:/1-UCRÏîÄ¿/CADD/UCR.ucname.#rs.CADD.loc", 
                   sep = "\t")
colnames(CADD) <- c("variation", "chr", "start", "end", 
                    "Allele", "CADD_PHRED", "CADD_RAW", "UCR")
summary(CADD$CADD_PHRED)

plotCADDDistribution <- 
  ggplot(data = CADD, mapping = aes(x = CADD_PHRED)) +
  geom_density() +
  
  myTheme

plotCADDDistribution

ggview::ggview(plot = plotCADDDistribution, 
               width = 8, height = 6, units = c("cm"), dpi = 1200)

# organize VCF files for CADD
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
write.table(snpsUCRsVCF, file = "02-analysis/48-SNP_Pathogenicity_Evaluation/snpsUCRsVCF.vcf", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(snpsUCRsVCF[c(1:1000), ], file = "02-analysis/48-SNP_Pathogenicity_Evaluation/snpsUCRsVCF-1-1000.vcf", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(snpsUCRsVCF[c(1001:2000), ], file = "02-analysis/48-SNP_Pathogenicity_Evaluation/snpsUCRsVCF-1001-2000.vcf", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(snpsUCRsVCF[c(2001:3000), ], file = "02-analysis/48-SNP_Pathogenicity_Evaluation/snpsUCRsVCF-2001-3000.vcf", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(snpsUCRsVCF[c(3001:4000), ], file = "02-analysis/48-SNP_Pathogenicity_Evaluation/snpsUCRsVCF-3001-4000.vcf", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(snpsUCRsVCF[c(4001:5000), ], file = "02-analysis/48-SNP_Pathogenicity_Evaluation/snpsUCRsVCF-4001-5000.vcf", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(snpsUCRsVCF[c(5001:6000), ], file = "02-analysis/48-SNP_Pathogenicity_Evaluation/snpsUCRsVCF-5001-6000.vcf", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(snpsUCRsVCF[c(6001:7000), ], file = "02-analysis/48-SNP_Pathogenicity_Evaluation/snpsUCRsVCF-6001-7000.vcf", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(snpsUCRsVCF[c(7001:8000), ], file = "02-analysis/48-SNP_Pathogenicity_Evaluation/snpsUCRsVCF-7001-8000.vcf", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(snpsUCRsVCF[c(8001:9000), ], file = "02-analysis/48-SNP_Pathogenicity_Evaluation/snpsUCRsVCF-8001-9000.vcf", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(snpsUCRsVCF[c(9001:10000), ], file = "02-analysis/48-SNP_Pathogenicity_Evaluation/snpsUCRsVCF-9001-10000.vcf", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(snpsUCRsVCF[c(10001:11000), ], file = "02-analysis/48-SNP_Pathogenicity_Evaluation/snpsUCRsVCF-10001-11000.vcf", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(snpsUCRsVCF[c(11001:12000), ], file = "02-analysis/48-SNP_Pathogenicity_Evaluation/snpsUCRsVCF-11001-12000.vcf", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(snpsUCRsVCF[c(12001:13000), ], file = "02-analysis/48-SNP_Pathogenicity_Evaluation/snpsUCRsVCF-12001-13000.vcf", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(snpsUCRsVCF[c(13001:14000), ], file = "02-analysis/48-SNP_Pathogenicity_Evaluation/snpsUCRsVCF-13001-14000.vcf", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(snpsUCRsVCF[c(14001:15000), ], file = "02-analysis/48-SNP_Pathogenicity_Evaluation/snpsUCRsVCF-14001-15000.vcf", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(snpsUCRsVCF[c(15001:16000), ], file = "02-analysis/48-SNP_Pathogenicity_Evaluation/snpsUCRsVCF-15001-16000.vcf", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(snpsUCRsVCF[c(16001:17000), ], file = "02-analysis/48-SNP_Pathogenicity_Evaluation/snpsUCRsVCF-16001-17000.vcf", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(snpsUCRsVCF[c(17001:18000), ], file = "02-analysis/48-SNP_Pathogenicity_Evaluation/snpsUCRsVCF-17001-18000.vcf", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(snpsUCRsVCF[c(18001:19000), ], file = "02-analysis/48-SNP_Pathogenicity_Evaluation/snpsUCRsVCF-18001-19000.vcf", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(snpsUCRsVCF[c(19001:20000), ], file = "02-analysis/48-SNP_Pathogenicity_Evaluation/snpsUCRsVCF-19001-20000.vcf", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(snpsUCRsVCF[c(20001:21000), ], file = "02-analysis/48-SNP_Pathogenicity_Evaluation/snpsUCRsVCF-20001-21000.vcf", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(snpsUCRsVCF[c(21001:22000), ], file = "02-analysis/48-SNP_Pathogenicity_Evaluation/snpsUCRsVCF-21001-22000.vcf", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(snpsUCRsVCF[c(22001:22299), ], file = "02-analysis/48-SNP_Pathogenicity_Evaluation/snpsUCRsVCF-22001-22299.vcf", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# organize vcf file for dbNSFP
# requirements
# 1. Please make sure you don't create a new line after your last query.
# 2. Please make sure use space between input parameters, don't use tab.
# 3. Upload and sumbit it at the box below.
write.table(snpsUCRsVCF[c(1:5000), ], 
            file = "02-analysis/48-SNP_Pathogenicity_Evaluation/dbNSFPSNPsUCRs-1-5000.vcf", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(snpsUCRsVCF[c(2001:4000), ], 
            file = "02-analysis/48-SNP_Pathogenicity_Evaluation/dbNSFPSNPsUCRs-2001-4000.vcf", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


snpsUCRsCADD <- read.table(
  file = gzfile
)
