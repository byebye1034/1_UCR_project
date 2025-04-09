# 对之前的物种的比对结果也进行新的identity计算

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(stringr)

UCR_location <- read.table(file = "01-data/UCR_raw/UCR_location_refseqid.txt", sep = "\t", header = TRUE)
UCR_length <- UCR_location[, c(4:5)]

# human -------------------------------------------------------------------

human <- read.table(file = "01-data/26-Evolution_newly_emerging_UCR_distribution/human_results.xls", sep = "\t")
colnames(human) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", 
                        "qend", "sstart", "send", "qseq", "sseq", "evalue", "bitscore")
human$sseqid <- gsub("::.*", "", human$sseqid)
human <- human %>% 
  group_by(sseqid) %>% 
  top_n(1, bitscore) %>% 
  top_n(1, qstart)
human <- human[, c(2, 4:6)]

human <- merge(human, UCR_length, by.x = "sseqid", by.y = "UCR_name")
human <- human %>% 
  mutate(identity = ((length.x - mismatch - gapopen - 1)/length.y)*100)

summary_human <- human %>% 
  summarise(
    ">95%" = sum(identity > 95),
    ">90%" = sum(identity > 90),
    ">80%" = sum(identity > 80)
  )

# 构建新的数据框
human_summary <- data.frame(
  species = rep("human", 3),
  identity = c(">95%", ">90%", ">80%"),
  number = unlist(summary_human)
)

save(human_summary, file = "02-analysis/26-Evolution_newly_emerging_UCR_distribution/human_summary.Rdata")

# summary -----------------------------------------------------------------

# 要替换的名称列表
species_list <- c("chimpanzee", "rat", "sheep", 
                  "horse", "mouse", "mochpri", 
                  "cow", "pig", "cat", "dog", "chicken", 
                  "turkey", "lizard", "x.tropicalis", 
                  "zebrafish", "fugu", "lungfish")

# 创建一个存储结果的列表
summary_1 <- human_summary

# 循环替换变量名并运行代码
for (species in species_list) {
  # 构建文件路径
  file_path <- paste0("01-data/26-Evolution_newly_emerging_UCR_distribution/", species, "_results.xls")
  
  # 读取数据并修改列名
  data <- read.table(file = file_path, sep = "\t")
  colnames(data) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", 
                      "qend", "sstart", "send", "qseq", "sseq", "evalue", "bitscore")
  data$sseqid <- gsub("::.*", "", data$sseqid)
  
  # 数据处理
  data <- data %>% 
    group_by(sseqid) %>% 
    top_n(1, bitscore) %>% 
    top_n(1, qstart) %>% 
    select(sseqid, length, mismatch, gapopen)
  
  data <- merge(data, UCR_length, by.x = "sseqid", by.y = "UCR_name")
  data <- data %>% 
    mutate(identity = ((length.x - mismatch - gapopen - 1)/length.y)*100)
  
  # 统计结果
  summary_data <- data %>% 
    summarise(
      ">95%" = sum(identity > 95),
      ">90%" = sum(identity > 90),
      ">80%" = sum(identity > 80)
    )
  
  # 构建新的数据框
  summary_df <- data.frame(
    species = rep(species, 3),
    identity = c(">95%", ">90%", ">80%"),
    number = unlist(summary_data)
  )
  
  # 将结果存储到列表中
  summary_1 <- rbind(summary_1, summary_df)
  
  # 保存结果到文件
  save(summary_df, file = paste0("02-analysis/26-Evolution_newly_emerging_UCR_distribution/", species, "_summary.Rdata"))
}

rownames(summary_1) <- NULL

c.elegans_summary <- data.frame(
  species = rep("c.elegans", 3),
  identity = c(">95%", ">90%", ">80%"),
  number = unlist(0, 0, 0)
)

d.melanogaster_summary <- data.frame(
  species = rep("d.melanogaster", 3),
  identity = c(">95%", ">90%", ">80%"),
  number = unlist(0, 0, 0)
)

s.cerevisiae_summary <- data.frame(
  species = rep("s.cerevisiae", 3),
  identity = c(">95%", ">90%", ">80%"),
  number = unlist(0, 0, 0)
)

summary_1 <- rbind(
  summary_1, 
  c.elegans_summary, 
  d.melanogaster_summary, 
  s.cerevisiae_summary
)

save(summary_1, file = "03-results/26-Evolution_newly_emerging_UCR_distribution/summary_1_before.Rdata")































