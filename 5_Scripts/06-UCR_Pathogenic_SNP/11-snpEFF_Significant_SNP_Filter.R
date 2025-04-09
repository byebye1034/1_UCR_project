# 20231125 之前使用的snpEFF对UCR上所有的SNP进行了注释，现在筛选出HIGH的SNP，用于下一步在VEP当中预测哪些与疾病相关

rm(list = ls())
setwd(dir = "D:/R_project/UCR_project")

library(tidyverse)

UCRSNPAnno <- read.table(file = "D:/A_projects/UCR/snpEFFAnno/UCRSNPAnno.vcf", 
                         sep = "\t", skip = 7)
colnames(UCRSNPAnno) <- c("CHROM",	"POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO")
HIGH <- str_detect(UCRSNPAnno$INFO, "HIGH")
UCRSNPAnnoHIGH <- UCRSNPAnno[HIGH, ]
UCRSNPAnnoHIGH <- UCRSNPAnnoHIGH %>% 
  separate(col = "INFO", into = c("Allele", "Annotation", "Annotation_Impact", "Gene_Name", 
                                  "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", 
                                  "Rank"), sep = "\\|", extra = "drop") # "\\|"是正则表达式用于匹配“|”
save(UCRSNPAnnoHIGH, file = "02-analysis/06-UCR_Pathogenic_SNP/UCRSNPAnnoHIGH.Rdata")

# 构造一个UCRSNPHIGH.vcf
load("D:/R_project/UCR_project/data/0B-SNPAnalysisData/UCRSNPAnnoHIGH.Rdata")
load("D:/R_project/UCR_project/data/filtered_UCR_SNP_passed_v2/UCRSNPVcf.Rdata")

UCRSNPHIGH <- UCRSNPVcf %>% 
  filter(ID %in% UCRSNPAnnoHIGH$ID)
write.table(UCRSNPHIGH, file = "02-analysis/06-UCR_Pathogenic_SNP/UCRSNPHIGH.vcf", 
            sep = "\t", quote = F, row.names = F, col.names = F)
# 将文件用notepad打开，把UCRSNPVcf.vcf里的前三行复制进来上传到服务器用vep预测功能




