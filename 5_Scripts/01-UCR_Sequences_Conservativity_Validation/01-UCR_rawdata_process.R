# UCR原始数据(hg16)下载，转换成GRCh38坐标，准备getfasta需要的bed文件

setwd("D:/R_project/UCR_project")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library(stringr)

# 导入从网站下载的UCR原始数据，转换成方便使用的格式 ----------------------------------------------

# 导入从网站下载的UCR原始数据
UCR_rawdata <- read.table(file = "D:/A_projects/UCR/UCR下载数据/UCR原始数据.csv", skip = 2, sep = ",")

# 定义列名
colnames(UCR_rawdata) <- c("name", "type", "length", "location(hg16)", 
                           "upstream gene distance", "upstream gene name", 
                           "within gene name", 
                           "downstream gene distance", "downstream gene name", 
                           "human mRNA", "human spliced EST", 
                           "any species mRNA", "any species EST", 
                           "UTR", "CDS", "intron", "10k upstream", "10k downstream", "intergenic", 
                           "RNAfold p-value", "chicken overlap", "chicken identical", 
                           "fugu overlap", "fugu identical")
# 保存为Rdata
save(UCR_rawdata, file = "D:/R_project/UCR_project/01-data/UCR_raw/UCR.Rdata")
load(file = "D:/R_project/UCR_project/01-data/UCR_raw/UCR.Rdata")

# 20230925 准备用于坐标转换的bed文件 ----------------------------------------------------------

UCR_bed <- data.frame(str_split_fixed(UCR_rawdata$`location(hg16)`, ":|-", 3))
write.table(UCR_bed, file = "UCR.bed", col.names = F, row.names = F, quote = F, sep = "\t")

# 在NCBI remap网站从hg16(NCBI34)转换成GRCh38.p14(hg38)
UCR_remapped_GRCh38 <- read.table(file = "remapped_UCR.bed")

# 导入雪寒师姐的ucr(GRCh38)列表
UCR_bed_xh <- read.table(file = "D:/1-UCR项目/UCR.bed.txt")

# 比较之后发现除了在UCR_remapped_GRCh38的数据最后多出来十个，其他的位置信息都是完全一致的
# 从这里开始使用雪寒师姐整理好的UCR_bed_xh文件进行后续处理
UCR_remapped_GRCh38 <- UCR_bed_xh

# 将uc名称放在最后一列，将染色体编号替换成NC的编号
UCR_remapped_GRCh38 <- select(UCR_remapped_GRCh38, V1, V3, V4, everything())

mapping <- read.table(file = "mapping.txt", sep = "\t")
mapping$V1 <- gsub("chromosome", "chr", mapping$V1)
mapping$V2 <- gsub(">", "", mapping$V2)
mapping$V1 <- gsub(" ", "", mapping$V1)
bed_for_getfasta <- UCR_remapped_GRCh38
bed_for_getfasta$V3 <- bed_for_getfasta$V3-1

bed_for_getfasta$V1 <- str_replace_all(bed_for_getfasta$V1, "\\bchr1\\b", "NC_000001.11")
bed_for_getfasta$V1 <- str_replace_all(bed_for_getfasta$V1, "\\bchr2\\b", "NC_000002.12")

for (i in 3:25)
{
  bed_for_getfasta$V1 <- str_replace_all(bed_for_getfasta$V1, mapping$V1[i], mapping$V2[i])
}

bed_for_getfasta <- bed_for_getfasta[order(bed_for_getfasta$V1), ]
write.table(bed_for_getfasta, file = "02-analysis/01-UCR_Conservativity_Validation/bed_for_getfasta.bed", 
            row.names = F, col.names = F, quote = F, sep = "\t")


# 20230926 blat之后的结果文件查看 --------------------------------------------------

setwd("D:/R_project/UCR_project/02-analysis/01-UCR_Conservativity_Validation/blat results/")
options(stringsAsFactors = F)

library(tidyverse)
library(stringr)

## 小鼠mouse
# 读入blat结果文件
blat_GRCm39.out <- read.table(
  file = "02-analysis/01-UCR_Conservativity_Validation/blat results/blat_GRCm39.out", 
  skip = 5, sep = "\t")

# 整理headers
header_f <- read_lines(
  file = "02-analysis/01-UCR_Conservativity_Validation/blat results/blat_GRCm39.out", 
  skip = 2, n_max = 1)

header_f <- str_split_fixed(header_f, "\t", 21)
header_f[2] <- "mis-match"
header_f[3] <- "repo.match"
header_f[5] <- "Q gap count"
header_f[6] <- "Q gap bases"
header_f[7] <- "T gap count"
header_f[8] <- "T gap bases"
header_f[10] <- "Q name"
header_f[11] <- "Q size"
header_f[12] <- "Q start"
header_f[13] <- "Q end"
header_f[14] <- "T name"
header_f[15] <- "T size"
header_f[16] <- "T start"
header_f[17] <- "T end"
header_f[18] <- "block count"
colnames(blat_GRCm39.out) <- header_f
GRCm39_complete_match <- blat_GRCm39.out
GRCm39_complete_match <- filter(GRCm39_complete_match, GRCm39_complete_match$'mis-match' == 0)

# 删除uc.188和uc.263，他们的match分别是30和44
GRCm39_part_match <- filter(GRCm39_complete_match, GRCm39_complete_match$'match' < 200)

# 481个ucr在小鼠基因组上全部都是100%保守
GRCm39_complete_match <- filter(GRCm39_complete_match, GRCm39_complete_match$'match' >= 200)

## 大鼠rat
blat_mRatBN7.2.out <- read.table(file = "blat_mRatBN7.2.out", sep = "\t", skip = 5)
header_rat <- read_lines(file = "blat_mRatBN7.2.out", skip = 2, n_max = 1)
header_rat <- str_split_fixed(header_rat, "\t", 21)
header_rat[2] <- "mis-match"
header_rat[3] <- "repo.match"
header_rat[5] <- "Q gap count"
header_rat[6] <- "Q gap bases"
header_rat[7] <- "T gap count"
header_rat[8] <- "T gap bases"
header_rat[10] <- "Q name"
header_rat[11] <- "Q size"
header_rat[12] <- "Q start"
header_rat[13] <- "Q end"
header_rat[14] <- "T name"
header_rat[15] <- "T size"
header_rat[16] <- "T start"
header_rat[17] <- "T end"
header_rat[18] <- "block count"
colnames(blat_mRatBN7.2.out) <- header_rat
mRatBN7.2_complete_match <- blat_mRatBN7.2.out
mRatBN7.2_part_match <- filter(mRatBN7.2_complete_match, mRatBN7.2_complete_match$'match' < 200 | mRatBN7.2_complete_match$'mis-match' > 0)
mRatBN7.2_complete_match <- filter(mRatBN7.2_complete_match, mRatBN7.2_complete_match$'mis-match' == 0 & mRatBN7.2_complete_match$'match' >= 200)
# 和大鼠基因组相比不完全保守的uc有两个，uc.18(231/238)和uc.304(270/272)






