# SNPs Pathogenicity Scores Distribution Visualization

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(stringr)
library(ggview)
library(ggprism)
library(sysfonts)
library(showtext)

library(data.table)

lwd_pt <- .pt*72.27/96

font_add(family = "Arial", regular = "C:/Windows/Fonts/arial.ttf")
font_families()

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

# FATHMM-XF ---------------------------------------------------------------

## UCRs
UCRsFATHMMXF <- read.delim(
  file = "02-analysis/48-SNP_Pathogenicity_Evaluation/FATHMM-XF/FATHMMXFScores/UCRs-FATHMMXF.txt", 
  sep = "\t", header = TRUE
)
colnames(UCRsFATHMMXF) <- c("CHROM", "POS", "REF", "ALT", 
                            "CodingScore", "NonCodingScore", "Warning")
table(UCRsFATHMMXF$Warning)
UCRsFATHMMXF <- UCRsFATHMMXF %>% 
  filter(grepl("benign|pathogenic", UCRsFATHMMXF$Warning)) %>% 
  mutate(Region = "UCRs")


## UCRs Left
UCRsLeftFATHMMXF <- read.delim(
  file = "02-analysis/48-SNP_Pathogenicity_Evaluation/SNPsUCRsLeft/snpsUCRsLeftFATHMMXF.txt", 
  sep = "\t", header = TRUE
)
colnames(UCRsLeftFATHMMXF) <- c("CHROM", "POS", "REF", "ALT", 
                                "CodingScore", "NonCodingScore", "Warning")
table(UCRsLeftFATHMMXF$Warning)
UCRsLeftFATHMMXF <- UCRsLeftFATHMMXF %>% 
  filter(grepl("benign|pathogenic", UCRsLeftFATHMMXF$Warning)) %>% 
  mutate(Region = "UCRsLeft")


## UCRs Right
UCRsRightFATHMMXF <- read.delim(
  file = "02-analysis/48-SNP_Pathogenicity_Evaluation/SNPsUCRsRight/snpsUCRsRightFATHMMXF.txt", 
  sep = "\t", header = TRUE
)
colnames(UCRsRightFATHMMXF) <- c("CHROM", "POS", "REF", "ALT", 
                                 "CodingScore", "NonCodingScore", "Warning")
table(UCRsRightFATHMMXF$Warning)
UCRsRightFATHMMXF <- UCRsRightFATHMMXF %>% 
  filter(grepl("benign|pathogenic", UCRsRightFATHMMXF$Warning)) %>% 
  mutate(Region = "UCRsRight")


## Randf
RandfFATHMMXF <- read.delim(
  file = "02-analysis/48-SNP_Pathogenicity_Evaluation/SNPsRandf/snpsRandfFATHMMXF.txt", 
  sep = "\t", header = TRUE
)
colnames(RandfFATHMMXF) <- c("CHROM", "POS", "REF", "ALT", 
                             "CodingScore", "NonCodingScore", "Warning")
table(RandfFATHMMXF$Warning)
RandfFATHMMXF <- RandfFATHMMXF %>% 
  filter(grepl("benign|pathogenic", RandfFATHMMXF$Warning)) %>% 
  mutate(Region = "Randf")

FATHMMXF <- rbind(
  UCRsFATHMMXF, UCRsLeftFATHMMXF, 
  UCRsRightFATHMMXF, RandfFATHMMXF
)
nrow(FATHMMXF[FATHMMXF$CodingScore == "--", ]) # 110646
nrow(FATHMMXF[FATHMMXF$NonCodingScore == "--", ]) # 6420
FATHMMXF <- FATHMMXF %>% 
  mutate(FATHMMXFScore = ifelse(
    CodingScore == "--", NonCodingScore, CodingScore
  ))
FATHMMXF$FATHMMXFScore <- as.numeric(FATHMMXF$FATHMMXFScore)

FATHMMXF$Region <- factor(x = FATHMMXF$Region, 
                          levels = c("UCRs", "UCRsLeft", "UCRsRight", "Randf"), 
                          ordered = TRUE)

save(FATHMMXF, file = "03-results/48-SNP_Pathogenicity_Evaluation/FATHMMXF.Rdata")

plotFATHMMXFDistribution <- ggplot(
  data = FATHMMXF, mapping = aes(x = FATHMMXFScore, 
                                 group = Region, color = Region)
) +
  geom_density(adjust = 1.5) +
  
  scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF", 
                                "#00A087FF", "#3C5488FF")) +
  
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  
  theme_prism(
    base_size = 7, 
    base_family = "Arial", 
    base_fontface = "plain", 
    base_line_size = 0.5/lwd_pt, 
    base_rect_size = 0.5/lwd_pt, 
    border = TRUE
  ) +
  
  myTheme

plotFATHMMXFDistribution

plotFATHMMXFDistribution + canvas(width = 8, height = 6, units = "cm", dpi = 1200)

pdf(file = "03-results/48-SNP_Pathogenicity_Evaluation/FATHMMXFScoreDistribution.pdf", 
    width = 800/254, height = 600/254)

showtext_begin()

plotFATHMMXFDistribution

showtext_end()

dev.off()

# ClinPred ----------------------------------------------------------------


# REVEL -------------------------------------------------------------------

REVEL <- fread(
  file = "01-data/48-SNP_Pathogenicity_Evaluation/revel-v1.3_all_chromosomes/revel_with_transcript_ids"
)
head(REVEL)
# > head(REVEL)
#   chr hg19_pos grch38_pos ref alt aaref aaalt REVEL Ensembl_transcriptid
# 1   1    35142      35142   G   A     T     M 0.027      ENST00000417324
# 2   1    35142      35142   G   C     T     R 0.035      ENST00000417324
# 3   1    35142      35142   G   T     T     K 0.043      ENST00000417324
# 4   1    35143      35143   T   A     T     S 0.018      ENST00000417324
# 5   1    35143      35143   T   C     T     A 0.034      ENST00000417324
# 6   1    35143      35143   T   G     T     P 0.039      ENST00000417324
class(REVEL)
REVEL <- REVEL[, c(1, 3:5, 8, 9)]
save(REVEL, file = "01-data/48-SNP_Pathogenicity_Evaluation/REVEL.Rdata")

load(file = "01-data/48-SNP_Pathogenicity_Evaluation/REVEL.Rdata")
head(REVEL)
REVEL$grch38_pos <- as.numeric(REVEL$grch38_pos)

## SNPsUCRs
SNPsUCRs <- read.delim(
  file = "02-analysis/48-SNP_Pathogenicity_Evaluation/FATHMM-XF/SNPs/snpsUCRsVCF4FATH.vcf", 
  sep = "\t", header = TRUE
)
head(SNPsUCRs)
# > head(SNPsUCRs)
# X.CHROM      POS           ID REF ALT QUAL FITTER INFO
# 1       1 10537642  rs989713876   C   T    .      .    .
# 2       1 10537645  rs909245667   C   T    .      .    .
# 3       1 10537646 rs1218943258   G   A    .      .    .
# 4       1 10537654 rs1396443709   C   A    .      .    .
# 5       1 10537655  rs921922938   C   T    .      .    .
# 6       1 10537657  rs533943014   G   T    .      .    .

SNPsUCRsREVELScore <- SNPsUCRs[, c(1:5)] %>%
  inner_join(REVEL, 
             by = c("X.CHROM" = "chr", 
                    "POS" = "grch38_pos", 
                    "REF" = "ref", 
                    "ALT" = "alt"))
save(SNPsUCRsREVELScore, 
     file = "02-analysis/48-SNP_Pathogenicity_Evaluation/SNPsUCRsREVELScore.Rdata")
