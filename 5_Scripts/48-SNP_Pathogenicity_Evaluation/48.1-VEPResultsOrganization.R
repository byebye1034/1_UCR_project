# organize the vep results

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(stringr)
library(cowplot)
library(ggview)
library(ggprism)
library(sysfonts)
library(showtext)

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

font_add(family = "Arial", regular = "C:/Windows/Fonts/arial.ttf")
font_families()

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

UCRsSNPsPathogenicityScores <- read.delim(
  file = "01-data/48-SNP_Pathogenicity_Evaluation/vepUCRs.txt", 
  sep = "\t", header = TRUE, 
)
UCRsSNPsPathogenicityScores <- UCRsSNPsPathogenicityScores[c(5:299691), ] # remove NA caused by X chr
UCRsLeftSNPsPathogenicityScores <- read.delim(
  file = "01-data/48-SNP_Pathogenicity_Evaluation/vepUCRsLeft.txt", 
  sep = "\t", header = TRUE
)
UCRsLeftSNPsPathogenicityScores <- UCRsLeftSNPsPathogenicityScores[c(2:350440), ]

UCRsRightSNPsPathogenicityScores <- read.delim(
  file = "01-data/48-SNP_Pathogenicity_Evaluation/vepUCRsRight.txt", 
  sep = "\t", header = TRUE
)
UCRsRightSNPsPathogenicityScores <- UCRsRightSNPsPathogenicityScores[c(3:364271), ]

RandfSNPsPathogenicityScores <- read.delim(
  file = "01-data/48-SNP_Pathogenicity_Evaluation/vepRandf.txt", 
  sep = "\t", header = TRUE
)

# get Scores
commonColnames <- intersect(colnames(UCRsSNPsPathogenicityScores), 
                            colnames(UCRsLeftSNPsPathogenicityScores))
commonColnames <- intersect(commonColnames, 
                            colnames(UCRsRightSNPsPathogenicityScores))
commonColnames <- intersect(commonColnames, 
                            colnames(RandfSNPsPathogenicityScores))

# CADD --------------------------------------------------------------------

UCRsCADD <- UCRsSNPsPathogenicityScores[, c(1, 44)]
UCRsLeftCADD <- UCRsLeftSNPsPathogenicityScores[, c(1, 47)]
UCRsRightCADD <- UCRsRightSNPsPathogenicityScores[, c(1, 47)]
RandfCADD <- RandfSNPsPathogenicityScores[, c(1, 46)]

length(unique(UCRsCADD$X.Uploaded_variation)) # 22231
length(unique(UCRsLeftCADD$X.Uploaded_variation)) # 25172
length(unique(UCRsRightCADD$X.Uploaded_variation)) # 25445
length(unique(RandfCADD$X.Uploaded_variation)) # 42869

UCRsCADD <- unique(UCRsCADD) # 26053
UCRsLeftCADD <- unique(UCRsLeftCADD) # 31050
UCRsRightCADD <- unique(UCRsRightCADD) # 31180
RandfCADD <- unique(RandfCADD) # 53741

UCRsCADD <- mutate(UCRsCADD, Region = "UCRs")
UCRsLeftCADD <- mutate(UCRsLeftCADD, Region = "UCRsLeft")
UCRsRightCADD <- mutate(UCRsRightCADD, Region = "UCRsRight")
RandfCADD <- mutate(RandfCADD, Region = "Randf")

CADDScore <- rbind(
  UCRsCADD, UCRsLeftCADD, 
  UCRsRightCADD, RandfCADD
)
CADDScore$CADD_PHRED <- as.numeric(CADDScore$CADD_PHRED)
CADDScore$Region <- factor(
  x = CADDScore$Region, 
  levels = c("UCRs", "UCRsLeft", "UCRsRight", "Randf"), 
  ordered = TRUE
)
save(CADDScore, file = "03-results/48-SNP_Pathogenicity_Evaluation/CADDScore.Rdata")

plotCADDDistribution <- ggplot(data = CADDScore) +
  
  geom_density(mapping = aes(x = CADD_PHRED, 
                             color = Region)) +
  
  scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF", 
                                "#00A087FF", "#3C5488FF")) +
  
  scale_x_continuous(limits = c(0, 60), expand = c(0, 0)) +
  
  guides(color = guide_legend(nrow = 2, ncol = 2)) +
  
  theme_prism(base_size = 7, 
              base_family = "Arial", 
              base_fontface = "plain", 
              base_line_size = 0.5/lwd_pt, 
              base_rect_size = 0.5/lwd_pt, 
              border = TRUE) +
  
  theme(
    panel.border = element_rect(linewidth = 0.5/lwd_pt, color = "#000000", fill = NA),
    legend.position = "top", 
    legend.key.size = unit(0.25, "cm"), 
    aspect.ratio = 1
  )

plotCADDDistribution

plotCADDDistribution + canvas(width = 4.8, height = 6, units = "cm", dpi = 1200)

pdf(file = "03-results/48-SNP_Pathogenicity_Evaluation/CADD_PHRED.pdf", 
    width = 480/254, height = 600/254)

showtext_begin()

plotCADDDistribution

showtext_end()

dev.off()











