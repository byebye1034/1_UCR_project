# 比较UCR和其他基因组位置的GC含量

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(dbplyr)

UCR_location <- read.table(file = "01-data/UCR_raw/UCR_location.txt", 
                           sep = "\t", header = TRUE)
UCR_location <- UCR_location[!(UCR_location$UCR_name == "uc.18" | UCR_location$UCR_name == "uc.304"), ]

load("D:/R_project/UCR_project/01-data/reference/chromosomes_and_corresponding_nc_numbers.Rdata")

UCR_location_temp <- merge(UCR_location, GRCh38.p14_report, 
                           by.x = "chr_name", by.y = "Chromosome.name", all.x = TRUE)
save(UCR_location_temp, file = "01-data/UCR_raw/UCR_location.Rdata")
write.table(UCR_location_temp, file = "01-data/UCR_raw/UCR_location_refseqid.txt", 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# 获得UCR left和right bed文件
write.table(UCR_location_temp[, c(10, 6:7, 4)], 
            file = "02-analysis/14-GC_content/UCR_left.bed", sep = "\t", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(UCR_location_temp[, c(10, 8:9, 4)], 
            file = "02-analysis/14-GC_content/UCR_right.bed", sep = "\t", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

# 获得rf的bed文件
rf_location <- read.table(file = "01-data/UCR_raw/random_fragments_location.txt", 
                          sep = "\t", header = TRUE)
rf_location <- rf_location[!(rf_location$rf_name == "rf.19" | rf_location$rf_name == "rf.323"), ]
rf_location_temp <- merge(rf_location, GRCh38.p14_report, 
                          by.x = "chr_name", by.y = "Chromosome.name", all.x = TRUE)
save(rf_location_temp, file = "01-data/UCR_raw/rf_location.Rdata")
write.table(rf_location_temp[, c(8, 3:4, 6)], 
            file = "02-analysis/14-GC_content/rf.bed", sep = "\t", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(rf_location_temp, file = "01-data/UCR_raw/random_fragments_location_refseqid.txt", 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# calculate GC content ----------------------------------------------------






