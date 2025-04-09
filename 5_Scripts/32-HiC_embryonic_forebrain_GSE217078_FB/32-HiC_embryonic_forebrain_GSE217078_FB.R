# 用GEO数据库查到的HiC数据可视化，看看uc.157和Otp的promoter之间有无相互作用

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)

FB_hic <- read.table(
  file = "01-data/32-HiC_embryonic_forebrain_GSE217078_FB/GSE217078_FB.ibed/GSE217078_FB.ibed", 
  sep = "\t", header = TRUE
)

FB_hic <- FB_hic[, c(1:3, 5:7, 10)]

FB_hic <- FB_hic %>% 
  mutate(temp = paste(otherEnd_chr, otherEnd_start, sep = ":")) %>% 
  mutate(temp_1 = paste(temp, otherEnd_end, sep = "-")) %>% 
  mutate(temp_2 = paste(temp_1, score, sep = ","))

write.table(FB_hic[, c(1:3, 10)], 
            file = "02-analysis/32-HiC_embryonic_forebrain_GSE217078_FB/forebrain_hic_longrange.bed", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

























