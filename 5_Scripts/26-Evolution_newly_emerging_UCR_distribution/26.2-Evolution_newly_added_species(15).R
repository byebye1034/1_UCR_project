# 在之前的19个物种（包括human，不包括mouse和rat）的基础上增加15个物种
# 文昌鱼 七鳃鳗 日本鲎 小型狗鱼 大白鲨 尼罗罗非鱼 大西洋鳕鱼 虹鳟
# 青蛙 平塔岛象龟 缅甸蟒
# 斑胸草稚 家鸽 鹦鹉
# 海豚

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library(stringr)

UCR_location <- read.table(file = "01-data/UCR_raw/UCR_location_refseqid.txt", sep = "\t", header = TRUE)
UCR_length <- UCR_location[, c(4:5)]

# 文昌鱼 ---------------------------------------------------------------------
# 0 没什么好说

wenchangyu <- read.table(file = "01-data/26.1-Evolution_newly_added_species/wenchangyu_results.xls")

wenchangyu_summary <- data.frame(
  species = rep("wenchangyu", 3),
  identity = c(">95%", ">90%", ">80%"),
  number = unlist(0, 0, 0)
)

# 七鳃鳗 ---------------------------------------------------------------------

qisaiman <- read.table(file = "01-data/26.1-Evolution_newly_added_species/qisaiman_results.xls", sep = "\t")
colnames(qisaiman) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", 
                        "qend", "sstart", "send", "qseq", "sseq", "evalue", "bitscore")
qisaiman$sseqid <- gsub("::.*", "", qisaiman$sseqid)
qisaiman <- qisaiman %>% 
  group_by(sseqid) %>% 
  top_n(1, bitscore) %>% 
  top_n(1, qstart)
qisaiman <- qisaiman[, c(2, 4:6)]

qisaiman <- merge(qisaiman, UCR_length, by.x = "sseqid", by.y = "UCR_name")
qisaiman <- qisaiman %>% 
  mutate(identity = ((length.x - mismatch - gapopen - 1)/length.y)*100)

summary_qisaiman <- qisaiman %>% 
  summarise(
    ">95%" = sum(identity > 95),
    ">90%" = sum(identity > 90),
    ">80%" = sum(identity > 80)
  )

# 构建新的数据框
qisaiman_summary <- data.frame(
  species = rep("qisaiman", 3),
  identity = c(">95%", ">90%", ">80%"),
  number = unlist(summary_qisaiman)
)

save(qisaiman_summary, file = "02-analysis/26.1-Evolution_newly_added_species/qisaiman_summary.Rdata")

# 0

# 日本鲎 ---------------------------------------------------------------------

# 0

ribenhou_summary <- data.frame(
  species = rep("ribenhou", 3),
  identity = c(">95%", ">90%", ">80%"),
  number = unlist(0, 0, 0)
)

# summary -----------------------------------------------------------------

# 要替换的名称列表
species_list <- c("haitun", "pingtadaoxianggui", "jiage", 
                  "banxiongcaoque", "miandianmang", "hongzun", 
                  "qingwa", "dabaisha", "xiaoxinggouyu", "niluoluofeiyu", 
                  "daxiyangxueyu")

# 创建一个存储结果的列表
summary_2 <- qisaiman_summary

# 循环替换变量名并运行代码
for (species in species_list) {
  # 构建文件路径
  file_path <- paste0("01-data/26.1-Evolution_newly_added_species/", species, "_results.xls")
  
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
  summary_2 <- rbind(summary_2, summary_df)
  
  # 保存结果到文件
  save(summary_df, file = paste0("02-analysis/26.1-Evolution_newly_added_species/", species, "_summary.Rdata"))
}

summary_2 <- rbind(
  summary_2, 
  ribenhou_summary, 
  wenchangyu_summary
)

rownames(summary_2) <- NULL
save(summary_2, file = "03-results/26-Evolution_newly_emerging_UCR_distribution/summary_2_added.Rdata")




