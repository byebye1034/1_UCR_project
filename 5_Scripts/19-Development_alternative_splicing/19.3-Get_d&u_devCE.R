# 获得所有pattern为down和up的devCE，查看exon上是否富集SRSF1的结合motif（hexamer）
# human_brain_down_devCE
# human_brain_up_devCE

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)

human_devAS <- read.table(
  file = "01-data/19-Development_alternative_splicing/human.devAS", 
  sep = ",", 
  header = TRUE
)

# get human brain devAS
human_brain_devAS <- human_devAS[!c(human_devAS$pattern.brain == "-" | human_devAS$pattern.brain == "n"), ]
human_brain_devAS <- human_brain_devAS[, c(1:9)]
print(table(human_brain_devAS$as.type))
#   AA      AD      CE complex      RI 
# 1445    1451    3693      91    2944
# 在human的brain所有的AS事件当中，CE和RI最多

print(table(human_brain_devAS$pattern.brain))
#    d   du    u   ud 
# 4993  959 2728  944
# 在human的brain所有的AS事件当中，down和up最多

# get human brain down&up devCE
human_brain_up_devCE <- human_brain_devAS[c(human_brain_devAS$as.type == "CE" & human_brain_devAS$pattern.brain == "u"), ]
human_brain_down_devCE <- human_brain_devAS[c(human_brain_devAS$as.type == "CE" & human_brain_devAS$pattern.brain == "d"), ]

# get human brain up&down devCE bed
write.table(human_brain_up_devCE[, c(4:7)], 
            file = "02-analysis/19-Development_alternative_splicing/human_brain_up_devCE.bed", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(human_brain_down_devCE[, c(4:7)], 
            file = "02-analysis/19-Development_alternative_splicing/human_brain_down_devCE.bed", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# then get fasta in linux










