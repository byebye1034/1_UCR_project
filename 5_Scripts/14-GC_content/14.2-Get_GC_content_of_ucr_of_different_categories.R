# 比较不同类别的UCR上的GC含量

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(Biostrings)
library(ggplot2)
library(ggpubr)

# get ucr type
load("D:/R_project/UCR_project/02-analysis/12-karyoplote/UCR_type_gene.Rdata")

# get ucr sequence
ucr_sequence <- readDNAStringSet(filepath = "01-data/UCR_raw/ucr_sequences_simple_name.fasta")
coding_ucr <- UCR_type_gene[(UCR_type_gene$UCR_type == "coding_UCR"), ]
ncRNA_ucr <- UCR_type_gene[(UCR_type_gene$UCR_type == "ncRNA_UCR"), ]
intergenic_ucr <- UCR_type_gene[(UCR_type_gene$UCR_type == "intergenic"), ]

# classify
coding_ucr_seq <- ucr_sequence[names(ucr_sequence) %in% coding_ucr$UCR_name, ]
ncRNA_ucr_seq <- ucr_sequence[names(ucr_sequence) %in% ncRNA_ucr$UCR_name, ]
intergenic_ucr_seq <- ucr_sequence[names(ucr_sequence) %in% intergenic_ucr$UCR_name, ]

# calculate GC content
coding_ucr_GC <- sapply(coding_ucr_seq, function(seq){
  sum(letterFrequency(seq, letters = c("G", "C"))) / sum(letterFrequency(seq, letters = c("G", "C", "A", "T", "N"))) * 100
})

ncRNA_ucr_GC <- sapply(ncRNA_ucr_seq, function(seq){
  sum(letterFrequency(seq, letters = c("G", "C"))) / sum(letterFrequency(seq, letters = c("G", "C", "A", "T", "N"))) * 100
})

intergenic_ucr_GC <- sapply(intergenic_ucr_seq, function(seq){
  sum(letterFrequency(seq, letters = c("G", "C"))) / sum(letterFrequency(seq, letters = c("G", "C", "A", "T", "N"))) * 100
})

# get GC content in different type ucr
GC_content <- data.frame(
  ucr_type = c(rep("coding", 314), rep("ncRNA", 66), rep("intergenic", 99)), 
  GC_content = c(coding_ucr_GC, ncRNA_ucr_GC, intergenic_ucr_GC)
)
save(GC_content, 
     file = "02-analysis/14-GC_content/GC_content_in_different_type_ucr/GC_content_in_different_type_ucr.Rdata")

# plot
rm(list = ls())

load("D:/R_project/UCR_project/02-analysis/14-GC_content/GC_content_in_different_type_ucr/GC_content_in_different_type_ucr.Rdata")
GC_content$ucr_type <- factor(x = GC_content$ucr_type, 
                              levels = c("coding", "ncRNA", "intergenic"), 
                              ordered = TRUE)

lwd_pt <- .pt*72.27/96

p1 <- ggplot(data = GC_content) +
  
  geom_violin(mapping = aes(x = ucr_type, 
                            y = GC_content, 
                            fill = ucr_type), 
              linewidth = 0.5/lwd_pt) +
  
  scale_fill_manual(values = c("#E64B35FF", 
                               "#4DBBD5FF", 
                               "#00A087FF"), 
                    name = "UCR type") +
  
  labs(x = "", 
       y = "GC content") +
  
  theme(
    panel.background = element_blank(), 
    panel.grid = element_blank(), 
    text = element_text(size = 7, color = "#000000"), 
    line = element_line(linewidth = 0.5/lwd_pt), 
    
    axis.text = element_text(size = 7, color = "#000000"), 
    axis.line = element_line(linewidth = 0.5/lwd_pt), 
    axis.ticks = element_line(linewidth = 0.5/lwd_pt), 
    
    legend.text = element_text(size = 7), 
    legend.position = "top",
  #  legend.key.width = unit(0.4, units = "cm"), 
    legend.key.size = unit(0.2, units = "cm"), 
    
    margin(r = 0, b = 0, l = 0), 
    aspect.ratio = 1:1
  )
p1

# 添加显著性标记
compare_list <- list(c("coding", "ncRNA"), 
                     c("ncRNA", "intergenic"), 
                     c("coding", "intergenic"))

p2 <- p1 + stat_compare_means(
  comparisons = compare_list, 
  method = "wilcox.test", 
  label = "p.signif", 
  bracket.size = 0.5/lwd_pt, 
  step.increase = 0.25
)
p2

ggsave(filename = "03-results/14-GC_content/GC_content_in_different_type_UCR.tiff", 
       width = 4.8, height = 5.8, units = "cm", dpi = 1200)

pdf(file = "03-results/14-GC_content/GC_content_in_different_type_UCR.pdf", 
    width = 480/254, height = 600/254)
p1
dev.off()











