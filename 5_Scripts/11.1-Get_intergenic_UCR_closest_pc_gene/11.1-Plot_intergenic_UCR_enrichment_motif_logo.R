# MEME找到一个富集在intergenic UCR上的motif，下面绘制一个logo

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(ggseqlogo)

motif <- read.table(
  file = "01-data/11.1-Get_intergenic_UCR_closest_pc_gene/CCMMKCTGVCWGCYK_counts.txt", 
  sep = " "
)















