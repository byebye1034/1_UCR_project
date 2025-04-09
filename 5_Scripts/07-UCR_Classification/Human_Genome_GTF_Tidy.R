# 对GRCh38基因组的基因组注释文件进行整理，用于之后的GO富集分析

rm(list = ls())
setwd(dir = "D:/R_project/UCR_project")

library(tidyverse)
library(data.table)

GRCh38.gtf <- fread(input = "02-analysis/07-UCR_Classification/Homo_sapiens.GRCh38.110.gtf", 
                    sep = "\t")
colnames(GRCh38.gtf) <- c("Chrom", "Source", "Feature", "Start", "End", 
                          "Score", "Strand", "Frame", "Attributes")

# 筛选注释文件当中的gene和exon，引物对于UCR的分类只需要这两个的位置信息
GRCh38.gtf_gene_exon_list <- filter(GRCh38.gtf, Feature == "gene" | Feature == "exon")
Gene_INFO <- filter(GRCh38.gtf_gene_exon_list, Feature == "gene")
Exon_INFO <- filter(GRCh38.gtf_gene_exon_list, Feature == "exon")

# 先处理基因信息 -----------------------------------------------------------------
With_Gene_Name <- str_detect(Gene_INFO$Attributes, "gene_name")
Gene_INFO_With_Name <- Gene_INFO[With_Gene_Name]
Gene_INFO_Without_Name <- Gene_INFO[!With_Gene_Name]


# Gene_INFO_With_Name ---------------------------------------------------
Gene_INFO_With_Name <- Gene_INFO_With_Name %>% 
  separate(col = "Attributes", into = c("GeneID", "GeneVersion", "GeneName", "GeneSource", "GeneBiotype"), 
           sep = ";", extra = "drop")
Gene_INFO_With_Name$GeneID <- gsub("gene_id", "", Gene_INFO_With_Name$GeneID)
Gene_INFO_With_Name$GeneVersion <- gsub("gene_version", "", Gene_INFO_With_Name$GeneVersion)
Gene_INFO_With_Name$GeneName <- gsub("gene_name", "", Gene_INFO_With_Name$GeneName)
Gene_INFO_With_Name$GeneSource <- gsub("gene_source", "", Gene_INFO_With_Name$GeneSource)
Gene_INFO_With_Name$GeneBiotype <- gsub("gene_biotype", "", Gene_INFO_With_Name$GeneBiotype)

# 去除双引号
for(i in 9:13){
  for(j in 1:nrow(Gene_INFO_With_Name)){
    x <- Gene_INFO_With_Name[j, i] # 赋值
    x <- as.character(x) # 化作字符串
    a <- gsub('["]', '', x) # 去双引号
    Gene_INFO_With_Name[j, i] <- a # 给数据框重新赋值
  }
}

# 去除空格
for(i in 9:13){
  for(j in 1:nrow(Gene_INFO_With_Name)){
    x <- Gene_INFO_With_Name[j, i] # 赋值
    x <- as.character(x) # 化作字符串
    a <- gsub(' ', '', x) # 去双引号
    Gene_INFO_With_Name[j, i] <- a # 给数据框重新赋值
  }
}

# Gene INFO Without Name --------------------------------------------------
Gene_INFO_Without_Name <- Gene_INFO_Without_Name %>% 
  separate(col = "Attributes", into = c("GeneID", "GeneVersion", "GeneSource", "GeneBiotype"), 
           sep = ";", extra = "drop")
Gene_INFO_Without_Name$GeneID <- gsub("gene_id", "", Gene_INFO_Without_Name$GeneID)
Gene_INFO_Without_Name$GeneVersion <- gsub("gene_version", "", Gene_INFO_Without_Name$GeneVersion)
Gene_INFO_Without_Name$GeneSource <- gsub("gene_source", "", Gene_INFO_Without_Name$GeneSource)
Gene_INFO_Without_Name$GeneBiotype <- gsub("gene_biotype", "", Gene_INFO_Without_Name$GeneBiotype)

# 去除双引号
for(i in 9:12){
  for(j in 1:nrow(Gene_INFO_Without_Name)){
    x <- Gene_INFO_Without_Name[j, i] # 赋值
    x <- as.character(x) # 化作字符串
    a <- gsub('["]', '', x) # 去双引号
    Gene_INFO_Without_Name[j, i] <- a # 给数据框重新赋值
  }
}

# 去除空格
for(i in 9:12){
  for(j in 1:nrow(Gene_INFO_Without_Name)){
    x <- Gene_INFO_Without_Name[j, i] # 赋值
    x <- as.character(x) # 化作字符串
    a <- gsub(' ', '', x) # 去双引号
    Gene_INFO_Without_Name[j, i] <- a # 给数据框重新赋值
  }
}

Gene_INFO_Without_Name <- Gene_INFO_Without_Name %>% 
  mutate(GeneName = ".")
Gene_INFO_Without_Name <- Gene_INFO_Without_Name[, c(1:10, 13, 11:12)]
Gene_INFO <- rbind(Gene_INFO_With_Name, Gene_INFO_Without_Name)
Gene_INFO <- arrange(Gene_INFO, Chrom, Start)
save(Gene_INFO, file = "data/0A-EvolutionAnalysisData/GRCh38.p14_Gene_INFO.Rdata")

# 再处理外显子信息 ----------------------------------------------------------------
# 也要区分有无基因名
With_Gene_Name <- str_detect(Exon_INFO$Attributes, "gene_name")
Exon_INFO_With_Name <- Exon_INFO[With_Gene_Name]
Exon_INFO_Without_Name <- Exon_INFO[!With_Gene_Name]

# 存在基因名的
Exon_INFO_With_Name <- Exon_INFO_With_Name %>% 
  separate(col = "Attributes", into = c("GeneID", "GeneVersion", "TranscriptID", "TranscriptVersion", 
                                        "ExonNum", "GeneName", "GeneSource", "GeneBiotype", 
                                        "TranscriptName", "TranscriptSource", "TranscriptBiotype"), 
           sep = ";", extra = "drop") # TranscriptBiotype之后的信息不一致就不要了
Exon_INFO_With_Name$GeneID <- gsub("gene_id", "", Exon_INFO_With_Name$GeneID)
Exon_INFO_With_Name$GeneVersion <- gsub("gene_version", "", Exon_INFO_With_Name$GeneVersion)
Exon_INFO_With_Name$TranscriptID <- gsub("transcript_id", "", Exon_INFO_With_Name$TranscriptID)
Exon_INFO_With_Name$TranscriptVersion <- gsub("transcript_version", "", Exon_INFO_With_Name$TranscriptVersion)
Exon_INFO_With_Name$ExonNum <- gsub("exon_number", "", Exon_INFO_With_Name$ExonNum)
Exon_INFO_With_Name$GeneName <- gsub("gene_name", "", Exon_INFO_With_Name$GeneName)
Exon_INFO_With_Name$GeneSource <- gsub("gene_source", "", Exon_INFO_With_Name$GeneSource)
Exon_INFO_With_Name$GeneBiotype <- gsub("gene_biotype", "", Exon_INFO_With_Name$GeneBiotype)
Exon_INFO_With_Name$TranscriptName <- gsub("transcript_name", "", Exon_INFO_With_Name$TranscriptName)
Exon_INFO_With_Name$TranscriptSource <- gsub("transcript_source", "", Exon_INFO_With_Name$TranscriptSource)
Exon_INFO_With_Name$TranscriptBiotype <- gsub("transcript_biotype", "", Exon_INFO_With_Name$TranscriptBiotype)

# 去除双引号
Exon_INFO_With_Name <- Exon_INFO_With_Name %>% 
  mutate_at(vars(9:19), ~str_remove_all(., '"'))

# 去除空格
Exon_INFO_With_Name <- Exon_INFO_With_Name %>% 
  mutate_at(vars(9:19), ~str_remove_all(., ' '))

# 不存在基因名的
Exon_INFO_Without_Name <- Exon_INFO_Without_Name %>% 
  separate(col = "Attributes", into = c("GeneID", "GeneVersion", "TranscriptID", "TranscriptVersion", 
                                        "ExonNum", "GeneSource", "GeneBiotype", "TranscriptSource", 
                                        "TranscriptBiotype"), 
           sep = ";", extra = "drop") # TranscriptBiotype之后的信息不一致就不要了
Exon_INFO_Without_Name$GeneID <- gsub("gene_id", "", Exon_INFO_Without_Name$GeneID)
Exon_INFO_Without_Name$GeneVersion <- gsub("gene_version", "", Exon_INFO_Without_Name$GeneVersion)
Exon_INFO_Without_Name$TranscriptID <- gsub("transcript_id", "", Exon_INFO_Without_Name$TranscriptID)
Exon_INFO_Without_Name$TranscriptVersion <- gsub("transcript_version", "", Exon_INFO_Without_Name$TranscriptVersion)
Exon_INFO_Without_Name$ExonNum <- gsub("exon_number", "", Exon_INFO_Without_Name$ExonNum)
Exon_INFO_Without_Name$GeneSource <- gsub("gene_source", "", Exon_INFO_Without_Name$GeneSource)
Exon_INFO_Without_Name$GeneBiotype <- gsub("gene_biotype", "", Exon_INFO_Without_Name$GeneBiotype)
Exon_INFO_Without_Name$TranscriptSource <- gsub("transcript_source", "", Exon_INFO_Without_Name$TranscriptSource)
Exon_INFO_Without_Name$TranscriptBiotype <- gsub("transcript_biotype", "", Exon_INFO_Without_Name$TranscriptBiotype)

Exon_INFO_Without_Name <- Exon_INFO_Without_Name %>% 
  mutate_at(vars(9:17), ~str_remove_all(., '"'))
Exon_INFO_Without_Name <- Exon_INFO_Without_Name %>% 
  mutate_at(vars(9:17), ~str_remove_all(., ' '))

Exon_INFO_Without_Name <- Exon_INFO_Without_Name %>% 
  mutate(GeneName = ".", TranscriptName = ".")
Exon_INFO_Without_Name <- Exon_INFO_Without_Name[, c(colnames(Exon_INFO_With_Name))]
Exon_INFO <- rbind(Exon_INFO_With_Name, Exon_INFO_Without_Name)
save(Exon_INFO, file = "02-analysis/07-UCR_Classification/GRCh38.p14_Exon_INFO.Rdata")

# 构造bed文件
Exon_INFO_bed <- Exon_INFO[, c(1, 4, 5, 9, 14)]
write.table(Exon_INFO_bed, file = "02-analysis/07-UCR_Classification/Exon_INFO.bed", 
            sep = "\t", quote = F, row.names = F, col.names = F)

load("D:/R_project/UCR_project/02-analysis/07-UCR_Classification/GRCh38.p14_Gene_INFO.Rdata")
Gene_INFO_bed <- Gene_INFO[, c(1, 4, 5, 9, 11)]
write.table(Gene_INFO_bed, file = "02-analysis/07-UCR_Classification/UCR classification/Gene_INFO.bed", 
            sep = "\t", quote = F, row.names = F, col.names = F)

# 构造含有GeneBiotype的bed文件
load("D:/R_project/UCR_project/02-analysis/07-UCR_Classification/GRCh38.p14_Gene_INFO.Rdata")
load("D:/R_project/UCR_project/02-analysis/07-UCR_Classification/GRCh38.p14_Exon_INFO.Rdata")

Gene_INFO_GeneBiotype_bed <- Gene_INFO[, c(1, 4, 5, 9, 11, 13)]
write.table(Gene_INFO_GeneBiotype_bed, file = "02-analysis/07-UCR_Classification/GRCh38.p14_Genome_INFO/Gene_INFO_GeneBiotype.bed", 
            sep = "\t", quote = F, row.names = F, col.names = F)

Exon_INFO_GeneBiotype_bed <- Exon_INFO[, c(1, 4, 5, 9, 11, 14, 16, 17)]
write.table(Exon_INFO_GeneBiotype_bed, file = "02-analysis/07-UCR_Classification/GRCh38.p14_Genome_INFO/Exon_INFO_GeneBiotype.bed", 
            sep = "\t", col.names = F, row.names = F, quote = F)





