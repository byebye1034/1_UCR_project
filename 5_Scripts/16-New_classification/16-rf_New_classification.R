# rf相关的分类

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(dbplyr)
library(tidyverse)

coding_rf <- read.table(file = "02-analysis/16-New_classification/01-rf_form_protein_coding_gene.bed", 
                        sep = "\t")
ncRNA_rf <- read.table(file = "02-analysis/16-New_classification/02-rf_from_Non_Coding_RNA.bed", sep = "\t")

coding_rf <- coding_rf[, c(1:4)]
ncRNA_rf <- ncRNA_rf[, c(1:4)]

coding_ncRNA_rf <- rbind(coding_rf, ncRNA_rf)
coding_ncRNA_rf <- unique(coding_ncRNA_rf)

rf <- read.table(file = "01-data/UCR_raw/random_fragments.bed", sep = "\t")
intergenic_rf <- anti_join(rf, coding_ncRNA_rf, by = "V4")
write.table(intergenic_rf, file = "02-analysis/16-New_classification/03-rf_intergenic.bed", 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# 按照新思路进行分类的rf
# 1、coding rf 234
# 2、ncRNA rf 67
# 3、intergenic rf
# 1和2重叠的部分有22个
# [1] "rf.8"   "rf.40"  "rf.65"  "rf.70"  "rf.74"  "rf.83"  "rf.98"  "rf.105" "rf.112" "rf.115" "rf.131" "rf.141"
# [13] "rf.143" "rf.160" "rf.173" "rf.174" "rf.183" "rf.188" "rf.199" "rf.234" "rf.254" "rf.388"

# 有5个intergenic rf和假基因重叠
# 1	258083	258637	rf.1	1	257864	359681	ENSG00000228463	.	transcribed_processed_pseudogene
# 1	104107839	104108385	rf.12	1	104072983	104226678	ENSG00000215869	.	transcribed_processed_pseudogene
# 1	143491143	143491645	rf.17	1	143461221	143498767	ENSG00000274927	KMT2CP3	unprocessed_pseudogene
# 12	38202959	38203486	rf.61	12	38201566	38203792	ENSG00000254381	TUBB8P5	unprocessed_pseudogene
# 14	64383267	64384003	rf.113	14	64341673	64387986	ENSG00000234911	TEX21P	transcribed_unitary_pseudogene
# 2	85102669	85103033	rf.243	2	85102389	85102694	ENSG00000225787	LSM3P3	processed_pseudogene



