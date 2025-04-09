# 获得random fragments的bed文件，然后去linux用bedtools获得相关基因

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)

rf <- read.table(file = "01-data/UCR_raw/random_fragments_location.txt", 
                      sep = "\t", header = TRUE)
rf <- rf[, c(5, 2, 3, 6)]
rf <- rf[!rf$rf_name %in% c("rf.19", "rf.323"),  ]
rf_bed <- write.table(rf, file = "02-analysis/13-Network_analysis/rf.bed", 
                      col.names = FALSE, row.names = FALSE, quote = FALSE, 
                      sep = "\t")






