# 寻找intergenic区域上的增强子位点

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(rtracklayer)

library(BiocManager)
BiocManager::install("Biostrings")

library(Biostrings)

# 读取bed文件
intergenic <- read.table(file = "01-data/10-GO_enrichment_analysis/intergenic_UCR_neighboring_gene.bed")
intergenic <- intergenic[, c(1:4)]

# 读取fasta文件
fasta <- readDNAStringSet(filepath = "01-data/UCR_raw/ucr_sequences_homemade_GRCh38p14.fa")
names(fasta) <- sub("::.*", "", names(fasta))                                   # 对fasta文件的名字进行处理

# 从FASTA文件中提取BED文件中指定的序列
intergenic_sequences <- fasta[intergenic$V4]

# 保存到一个新的FASTA文件
writeXStringSet(intergenic_sequences, "02-analysis/11-intergenic_region_enhancer_site/intergenic_sequences.fasta")












