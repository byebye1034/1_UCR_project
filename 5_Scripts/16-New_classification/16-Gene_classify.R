# 先获得两个大类：1、与蛋白编码基因重叠的coding_UCR；2、和非编码RNA重叠的ncRNA_UCR（只有一个不是lncRNA）
# 

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(dbplyr)

Gene_Non_Pr <- read.table(file = "01-data/16-New_classification/Gene_Non_Protein_Coding.bed", 
                          sep = "\t")
Gene_ncRNA <- Gene_Non_Pr[grep("RNA", Gene_Non_Pr$V6), ]

write.table(Gene_ncRNA, file = "01-data/16-New_classification/Gene_Non_Coding_RNA.bed", 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# noncoding region --------------------------------------------------------

Gene_Non_Pr <- read.table(file = "01-data/16-New_classification/Gene_Non_Protein_Coding.bed", 
                          sep = "\t")

# 查看除coding和RNA之外，基因组注释里还有什么
Gene_NonPr_NonRNA <- Gene_Non_Pr[!grepl("RNA", Gene_Non_Pr$V6), ]
table(Gene_NonPr_NonRNA$V6)
write.table(Gene_NonPr_NonRNA, file = "01-data/16-New_classification/Gene_NonPc_NonncRNA.bed", 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# 按照新思路对UCR进行分类
# 1、coding_UCR 315
# 2、ncRNA_UCR 101
# 3、intergenic_UCR 63
# note:在intergenic当中uc.329和一个假基因重叠
# 11	32176446	32176752	uc.329	11	32112049	32343409	ENSG00000227160	THEM7P	transcribed_unitary_pseudogene




