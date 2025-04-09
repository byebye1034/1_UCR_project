# NG已经定义了四类devAS，统计一下每类devAS的数目
# 获得所有pattern为“d”的devAS，查看这些AS上是否富集SRSF1等的结合motif

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)

human_devAS <- read.table(file = "01-data/19-Development_alternative_splicing/human.devAS", 
                          sep = ",", header = TRUE)
print(table(human_devAS$as.type))
#    AA      AD      CE complex      RI 
# 14876   12787   31134     881   34100

# get brain dynamic AS
human_brain_devAS <- human_devAS[!(c(human_devAS$pattern.brain == "-" | human_devAS$pattern.brain == "n")), ]
human_brain_devAS <- human_brain_devAS[, c(1:9, 16)]

# 1、首先分析brain当中所有devAS的pattern：d类型占了一半以上
table(human_brain_devAS$pattern.brain)
#    d   du    u   ud 
# 4993  959 2728  944

# 2、分析AS的类型分布：主要是CE和RI
table(human_brain_devAS$as.type)
#   AA      AD      CE complex      RI 
# 1445    1451    3693      91    2944

# get brain down devAS
human_brain_down_devAS <- human_brain_devAS[(human_brain_devAS$pattern.brain == "d"), ]
table(human_brain_down_devAS$as.type)
# AA      AD      CE complex      RI 
# 818     813    1412      53    1897

# get human brain down AS bed
write.table(human_brain_down_devAS[, c(4:7)], 
            file = "02-analysis/19-Development_alternative_splicing/human_brain_down_devAS.bed", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# add chr
human_brain_down_devAS <- read.table(file = "02-analysis/19-Development_alternative_splicing/human_brain_down_devAS.bed")
human_brain_down_devAS$V1 <- paste0("chr", human_brain_down_devAS$V1)
write.table(human_brain_down_devAS, 
            file = "02-analysis/19-Development_alternative_splicing/human_brain_down_devAS_chr.bed", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# get brain up devAS
human_brain_up_devAS <- human_brain_devAS[(human_brain_devAS$pattern.brain == "u"), ]
table(human_brain_up_devAS$as.type)
#  AA      AD      CE complex      RI 
# 363     340    1531      13     481

# get human brain up devAS bed
write.table(human_brain_up_devAS[, c(4:7)], 
            file = "02-analysis/19-Development_alternative_splicing/human_brain_up_devAS.bed", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# add chr
human_brain_up_devAS$chr_id <- paste0("chr", human_brain_up_devAS$chr_id)
write.table(human_brain_up_devAS[, c(4:7)], 
            file = "02-analysis/19-Development_alternative_splicing/human_brain_up_devAS_chr.bed", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# a\ get human brain down devAS bed(with strand +/-)
#    get human brain up devAS bed(with strand +/-)

# b\ bedtools getfasta

# c\ count hexamers summary


