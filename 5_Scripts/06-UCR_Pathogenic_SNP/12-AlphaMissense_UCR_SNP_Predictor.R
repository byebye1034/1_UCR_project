# snpEFF对UCR上的SNP的注释结果只有138个HIGH太低了，现在用AlphaMissense来试试，看看会不会多一些

rm(list = ls())
setwd(dir = "D:/R_project/UCR_project/")

library(tidyverse)
library(data.table)

# 准备用于VEP注释的vcf文件
UCRSNPPathogenic <- read.table(file = "02-analysis/06-UCR_Pathogenic_SNP/UCRSNPPathogenic.bed", 
                               sep = "\t")
colnames(UCRSNPPathogenic) <- c("SNPChr", "SNPStart", "SNPEnd", 
                                "REF", "ALT", "Genome", "UNIPROTID", "TRANSCRIPTID", 
                                "PROVAR", "PATHOGENICITY", "PATHO", "UCRChr", "UCRStart", "UCREnd", "UCRName")
UCRSNPPathogenicVCF <- UCRSNPPathogenic[, c(1:5, 15)]
colnames(UCRSNPPathogenicVCF) <- c("#CHROM", "START", "END", "REF", "ALT", "INFO")
UCRSNPPathogenicVCF$QUAL <- 255
UCRSNPPathogenicVCF$FILTER <- "PASS"
UCRSNPPathogenicVCF <- UCRSNPPathogenicVCF[, c(1:5, 7, 8, 6)]
write.table(UCRSNPPathogenicVCF, file = "02-analysis/06-UCR_Pathogenic_SNP/UCRSNPPatho.vcf", sep = "\t", 
            col.names = T, row.names = F, quote = F)

# 统计病理性的SNP在UCR上的比例
UCR_location <- read.table(file = "01-data/UCR_raw/UCR_location.txt", 
                          sep = "\t", header = T)

# 使用group_by和summarize进行分组和汇总
result <- UCRSNPPathogenic %>%
  group_by(UCRName) %>%
  summarize(SNPs_count = n(), Unique_SNPs_count = n_distinct(SNPStart)) # 结果将是一个新的数据框，其中包含每个基因组位置（UCRName）上发生SNP的位点数量和不重复的SNP的数量

UCR_length <- UCR_location[, c(4:5)]
UCR_length <- filter(UCR_length, UCR_length$UCR_name %in% result$UCRName)
colnames(result)[1] <- "UCR_name"
result <- merge(result, UCR_length, by = "UCR_name")
result <- result %>% 
  mutate(SNP_proportion = Unique_SNPs_count/length)
save(result, file = "02-analysis/06-UCR_Pathogenic_SNP/Patho_SNP_Percentage_in_UCR.Rdata")

# 下面这几行是干嘛的我也不清楚了
UCR_genes_type <- UCR_genes_type %>%
mutate(OCRLocation = paste(UCR_genes_type$UCRChr, ":", UCR_genes_type$UCRStart, "-", UCR_genes_type$UCREnd))
UCR_genes_type$OCRLocation <- gsub(" ", "", UCR_genes_type$OCRLocation)

UCR_pathogenic_variant_GRCh38 <- read.table(file = "02-analysis/06-UCR_Pathogenic_SNP/UCR_pathogenic_variant_GRCh38.bed", 
                                            sep = "\t")

# UCR上SNP的分布 --------------------------------------------------------------

# 现在要搞清楚UCR上的SNP是否主要分布在蛋白编码基因，是在外显子还是内含子？
# AlphaMissense预测的SNP的范围：是否包括内含子？
# 思路：构造一个包含所有SNP的bed文件，在nebula里用bedtools和exon_info做交集

library(tidyverse)

load("D:/R_project/UCR_project/02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/passed_UCR_SNPs.Rdata")
validated_UCR_SNPs <- passed_UCR_SNPs
validated_UCR_SNPs$chrom_start <- as.numeric(validated_UCR_SNPs$chrom_start)
validated_UCR_SNPs <- validated_UCR_SNPs %>% 
  mutate(chrom_end = chrom_start + 1)
validated_UCR_SNPs_bed <- validated_UCR_SNPs[, c(5, 7, 10, 1, 2, 6)]
write.table(validated_UCR_SNPs_bed, file = "02-analysis/06-UCR_Pathogenic_SNP/validated_UCR_SNP.bed", 
            sep = "\t", col.names = F, row.names = F, quote = F)

# 在nebula使用bedtools
# bedtools intersect -a validated_UCR_SNP.bed -b Exon_INFO_GeneBiotype.bed -wa -wb > SNP_within_Exon.bed

# UCR上发生的SNP大部分不是在exon上
SNP_within_exon <- read.table(file = "02-analysis/06-UCR_Pathogenic_SNP/SNP_within_Exon.bed", quote = "\t")
unique_SNP_within_exon <- SNP_within_exon[!duplicated(SNP_within_exon$V5), ]

# 结论
# 和外显子有交集的UCR：142
# 包含位于外显子上的SNP的UCR：140
# 在外显子上的SNP有5538个
# lncRNA：649 miRNA：11 protein_coding：4878
# UCR上SNP的总数目：22299
# 
# 













