# 这里不止有获取SNP，还包括了第一步V1的过滤，也就是存在evidence

# snp获取 -------------------------------------------------------------------

# first methods to obtain snps --------------------------------------------

rm(list = ls())
setwd(dir = "D:/R_project/UCR_project/02-analysis/03-SNP_INFO_Retrieval")

library(tidyverse)
library(biomaRt)

# 1. 准备ucr的位置信息文件
UCR_location <- read.table(file = "UCR_remapped_renamed_for_getfasta.bed", sep = "\t")
colnames(UCR_location) <- c("chr_name", "start", "end", "UCR_name")
UCR_location$chr_name <- gsub("chr", "", UCR_location$chr_name)

# 2. 连接到Ensembl数据库，选择适当的数据库和数据集：
ensembl <- useMart(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

# 3. 使用 UCR_location 数据框中的信息作为查询条件：
ucr_snp_info <- list()  # 创建一个空的数据框来存储结果
for (i in 1:length(UCR_location)) {
  chr_name <- UCR_location$chr_name[i]
  start <- UCR_location$start[i]
  end <- UCR_location$end[i]
  
  # 获取 SNP 信息
  snp_info <- getBM(attributes = c('refsnp_id', 'refsnp_source', 'refsnp_source_description', 'chr_name','allele','chrom_start','chrom_strand', 'validated'),
                    filters = c('chr_name','start','end'),
                    values = list(chr_name, start, end),
                    mart = ensembl)
  
  UCR_name <- UCR_location$UCR_name[i]
  ucr_snp_info[[UCR_name]] <- snp_info
} # 这种方法只能获得uc.1-uc.472的数据，不知道为什么

# 472个ucr的snp数据也不是一次性获得的，一次性获得那么多数据是不允许的，所以每次获得20个，并合并为一个list

save(ucr_snp_info, file = "data/snp/ucr_snp_info.Rdata")# 471
load(file = "data/snp/ucr_snp_info.Rdata")# 471

# 将在ensembl网站直接获得的SNP信息与利用biomaRt获得的信息整合成complete_ucr_snp_info.Rdata
file_names <- list.files("D:/R_project/UCR_project/data/snp", pattern = "csv$")
file_names <- paste0("D:/R_project/UCR_project/data/snp/", file_names)
ucr_snp_info_473_482 <- list()
for (file in file_names){
  snp_info <- read.csv(file = file, header = T)
  snp_info <- snp_info %>% mutate(chr_name = "X")
  snp_info <- snp_info[, c(1, 10, 9, 17, 5, 3, 11)]
  
  uc_name <- str_sub(file, start = 35, end = 40)
  ucr_snp_info_473_482[[uc_name]] <- snp_info
}

# merge two list
ucr_snp_info <- c(ucr_snp_info, ucr_snp_info_473_482)
save(ucr_snp_info, file = "D:/R_project/UCR_project/data/snp/complete_ucr_snp_info.Rdata")

# second methods to obtain snps ----------------------------------------------

# 因为获取的数据太大，直接用getBM获得不了

library(httr)
library(jsonlite)

getVariantsInRegion <- function(chrom, start, end) {
  
  region <- sprintf("%s:%s-%s",chrom,start,end)
  server <- "https://rest.ensembl.org/overlap/region/human/"
  extension <- "?feature=variation"
  full_url <- paste0(server, region, extension)
  
  result <- GET(full_url, content_type("application/json"))
  stop_for_status(result)
  
  fromJSON(toJSON(content(result)))
}

# vars <- getVariantsInRegion(chrom = "X", start = "71248993", end = "71249202")
# full_url <- https://rest.ensembl.org/overlap/region/human/1:10537640-10537846?feature=variation

ucr_snp_info_2 <- list()
for (i in 473:481) {
  chr_name <- UCR_location$chr_name[i]
  start <- UCR_location$start[i]
  end <- UCR_location$end[i]
  
  # 获取 SNP 信息
  snp_info <- getVariantsInRegion(chrom = chr_name, start = start, end = end)
  
  UCR_name <- UCR_location$UCR_name[i]
  ucr_snp_info_2[[UCR_name]] <- snp_info
}

save(ucr_snp_info_2, file = "data/snp/ucr_snp_info(473-481).R")
load(file = "data/snp/ucr_snp_info(473-481).R")

# get UCR flanking genomic positions --------------------------------------

rm(list = ls())
setwd(dir = "D:/R_project/UCR_project/02-analysis/03-SNP_INFO_Retrieval")

library(tidyverse)

# 1. 准备ucr的位置信息文件
UCR_location <- read.table(file = "UCR_remapped_renamed_for_getfasta.bed", sep = "\t")
colnames(UCR_location) <- c("chr_name", "start", "end", "UCR_name")
UCR_location$chr_name <- gsub("chr", "", UCR_location$chr_name)

# 2. 计算每个ucr的长度
UCR_location <- UCR_location %>% mutate(length = end - start)

# 3. 获得ucr左侧 ,右侧序列位置信息
UCR_location <- UCR_location %>% 
  mutate(left_start = start - length - 1, left_end = start - 1)
UCR_location <- UCR_location %>% 
  mutate(right_start = end + 1, right_end = end + length + 1)

# 4. 获得左侧，右侧序列的snp信息
## left
ensembl <- useMart(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")


# get complete_ucr_left_snp_info.Rdata ------------------------------------

ucr_snp_info_left <- list()  # 创建一个空的数据框来存储结果

for (i in 1:471) {
  chr_name <- UCR_location$chr_name[i]
  start <- UCR_location$left_start[i]
  end <- UCR_location$left_end[i]
  
  # 获取 SNP 信息
  snp_info <- getBM(attributes = c('refsnp_id', 'refsnp_source', 'refsnp_source_description', 'chr_name','allele','chrom_start','chrom_strand', 'validated'),
                    filters = c('chr_name','start','end'),
                    values = list(chr_name, start, end),
                    mart = ensembl)
  
  UCR_name <- UCR_location$UCR_name[i]
  ucr_snp_info_left[[UCR_name]] <- snp_info
}

save(ucr_snp_info_left, file = "flanking_snps/ucr_left_snp_info.Rdata") # 471
load(file = "flanking_snps/ucr_left_snp_info.Rdata")

# merge
file_names <- list.files("D:/R_project/UCR_project/data/flanking_snps/ucr_left_snp_info_473-482", pattern = "csv$")
file_names <- paste0("D:/R_project/UCR_project/data/flanking_snps/ucr_left_snp_info_473-482/", file_names)
ucr_left_snp_info_473_482 <- list()
for (file in file_names){
  snp_info <- read.csv(file = file, header = T)
  snp_info <- snp_info %>% mutate(chr_name = "X")
  snp_info <- snp_info[, c(1, 10, 9, 17, 5, 3, 11)]
  
  uc_name <- str_sub(file, start = 71, end = 76)
  ucr_left_snp_info_473_482[[uc_name]] <- snp_info
}

ucr_left_snp_info <- c(ucr_snp_info_left, ucr_left_snp_info_473_482)
save(ucr_left_snp_info, file = "D:/R_project/UCR_project/02-analysis/03-SNP_INFO_Retrieval/flanking_snpscomplete_ucr_left_snp_info.Rdata")


# get complete_ucr_right_snp_info.Rdata -----------------------------------

library(biomaRt)
ensembl <- useMart(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

## right
ucr_snp_info_right <- list()  # 创建一个空的数据框来存储结果
for (i in 472) {
  chr_name <- UCR_location$chr_name[i]
  start <- UCR_location$right_start[i]
  end <- UCR_location$right_end[i]
  
  # 获取 SNP 信息
  snp_info <- getBM(attributes = c('refsnp_id', 'refsnp_source', 'refsnp_source_description', 'chr_name','allele','chrom_start','chrom_strand', 'validated'),
                    filters = c('chr_name','start','end'),
                    values = list(chr_name, start, end),
                    mart = ensembl)
  
  UCR_name <- UCR_location$UCR_name[i]
  ucr_snp_info_right[[UCR_name]] <- snp_info
}

save(ucr_snp_info_right, file = "flanking_snps/ucr_right_snp_info_1-471.Rdata") #471
load(file = "flanking_snps/ucr_right_snp_info_1-471.Rdata")

# merge right
file_names <- list.files("flanking_snps/ucr_right_snp_info_473-482/", pattern = "csv$")
file_names <- paste0("D:/R_project/UCR_project/02-analysis/03-SNP_INFO_Retrieval/flanking_snps/ucr_right_snp_info_473-482/", file_names)
ucr_right_snp_info_473_482 <- list()
for (file in file_names){
  snp_info <- read.csv(file = file, header = T)
  snp_info <- snp_info %>% mutate(chr_name = "X")
  snp_info <- snp_info[, c(1, 10, 9, 17, 5, 3, 11)]
  
  uc_name <- str_sub(file, start = 72, end = 77)
  ucr_right_snp_info_473_482[[uc_name]] <- snp_info
}

ucr_right_snp_info <- c(ucr_snp_info_right, ucr_right_snp_info_473_482)
save(ucr_right_snp_info, file = "D:/R_project/UCR_project/02-analysis/03-SNP_INFO_Retrieval/flanking_snps/complete_ucr_right_snp_info.Rdata")

# get random fragments and snp info ----------------------------------------------------

library("tidyverse")

seq_report <- read.table(file = "GRCh38_p14_sequence_report.tsv", sep = "\t", header = T)
genome_range <- seq_report[1:24, ]
genome_range <- genome_range[, c(13, 12)]
genome_range <- genome_range %>% 
  mutate(chr_start = c(rep(1, 24)))
colnames(genome_range)[2] <- "chr_end"
genome_range <- genome_range[, c(1, 3, 2)]

# 设置随机数生成种子以确保结果可重复
set.seed(123)

# 创建一个空的数据框来存储选择的片段
selected_segments <- data.frame(UCSC.style.name = character(0),
                                Start = numeric(0),
                                End = numeric(0),
                                Length = numeric(0))

# 需要选取的总片段数
total_segments <- 481

while (nrow(selected_segments) < total_segments) {
  # 随机选择一个染色体
  chromosome_index <- sample(1:nrow(genome_range), 1)
  chromosome <- genome_range[chromosome_index, "UCSC.style.name"]
  chr_start <- genome_range[chromosome_index, "chr_start"]
  chr_end <- genome_range[chromosome_index, "chr_end"]
  
  # 随机生成片段的长度
  segment_length <- sample(200:779, 1)
  
  # 随机生成片段的起始位置
  start_position <- sample(chr_start:(chr_end - segment_length), 1)
  
  # 计算片段的结束位置
  end_position <- start_position + segment_length
  
  # 添加选定的片段到结果数据框
  new_segment <- data.frame(UCSC.style.name = chromosome,
                            Start = start_position,
                            End = end_position,
                            Length = segment_length)
  selected_segments <- rbind(selected_segments, new_segment)
}

# 打印结果
save(selected_segments, file = "data/flanking_snps/random_fragmens.Rdata") #这个没有保存，保存了line 264的data

# get snps on random fragments
selected_segments$chr_name <- gsub("chr", "", selected_segments$UCSC.style.name)
library(tidyverse)
selected_segments <- selected_segments %>% arrange(chr_name, Start)

rf_name <- c(1:481)
rf_name <- paste("rf.", rf_name, sep = "")

selected_segments <- selected_segments %>% mutate(rf_name = rf_name)
save(selected_segments, file = "flanking_snps/random_fragmens.Rdata")

# 在后面从bravo获得所有TOPMed的snp发挥作用
# selected_segments <- selected_segments %>%
#   mutate(location = paste0(chr_name, ":", Start, "-", End))
# write.table(selected_segments, file = "D:/R_project/UCR_project/data/flanking_snps/random_fragments_location.txt", 
#             row.names = F, quote = F, sep = "\t")

ensembl <- useMart(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
random_fragments_snp_info <- list()
for (i in 449:454) {
  chr_name <- selected_segments$chr_name[i]
  start <- selected_segments$Start[i]
  end <- selected_segments$End[i]
  
  # 获取 SNP 信息
  snp_info <- getBM(attributes = c('refsnp_id', 'refsnp_source', 'refsnp_source_description', 'chr_name','allele','chrom_start','chrom_strand', 'validated'),
                    filters = c('chr_name','start','end'),
                    values = list(chr_name, start, end),
                    mart = ensembl)
  
  rf_name <- selected_segments$rf_name[i]
  random_fragments_snp_info[[rf_name]] <- snp_info
}

save(random_fragments_snp_info, file = "flanking_snps/random_fragments_snp_info.Rdata") # 1-426,430-447,449-454(共450个elements)
load(file = "flanking_snps/random_fragments_snp_info.Rdata")

# merge random fragments snp info
file_names <- list.files("D:/R_project/UCR_project/02-analysis/03-SNP_INFO_Retrieval/flanking_snps/random_fragments_snp_info_complementary/", pattern = "csv$")
file_names <- paste0("D:/R_project/UCR_project/02-analysis/03-SNP_INFO_Retrieval/flanking_snps/random_fragments_snp_info_complementary/", file_names)
random_fragments_snp_info_temp <- list()
for (file in file_names){
  snp_info <- read.csv(file = file, header = T)
  snp_info <- snp_info %>% mutate(chr_name = "X")
  snp_info <- snp_info[, c(1, 10, 9, 17, 5, 3, 11)]
  
  uc_name <- str_sub(file, start = 85, end = 90)
  random_fragments_snp_info_temp[[uc_name]] <- snp_info
}

random_fragments_snp_info <- c(random_fragments_snp_info, random_fragments_snp_info_temp)
save(random_fragments_snp_info, file = "D:/R_project/UCR_project/02-analysis/03-SNP_INFO_Retrieval/flanking_snps/complete_random_fragments_snp_info.Rdata")

# snp summary -------------------------------------------------------------

rm(list = ls())
setwd(dir = "D:/R_project/UCR_project/02-analysis/03-SNP_INFO_Retrieval")

library(tidyverse)
library(ggplot2)

load(file = "snp/complete_ucr_snp_info.Rdata")
load("flanking_snps/complete_ucr_right_snp_info.Rdata")
load("/flanking_snps/complete_ucr_left_snp_info.Rdata")
load("/flanking_snps/complete_random_fragments_snp_info.Rdata")

# filter_evidence_SNP -----------------------------------------------------

library(tidyverse)

# filter_ucr_snp_info ------------------------------------------------------------

# 将chr_name列的数据类型设置为character
ucr_snp_info <- lapply(ucr_snp_info, function(df) {
  if ("validated" %in% colnames(df)) {
    df$chr_name <- as.character(df$chr_name)
  } else if ("Evidence" %in% colnames(df)) {
    df$chr_name <- as.character(df$chr_name)
  }
  return(df)
})

# 合并数据框
combined_ucr_snp_info <- bind_rows(ucr_snp_info, .id = "list_type")

# 将所有空白值修改为缺失值
combined_ucr_snp_info <- combined_ucr_snp_info %>% mutate_all(~ ifelse(. == "", NA, .))

# 筛选validated和Evidence不同时为NA的行
filtered_combined_ucr_snp_info <- combined_ucr_snp_info %>%
  filter(!is.na(validated) | !is.na(Evidence))

save(filtered_combined_ucr_snp_info, file = "filtered_UCR_SNP_with_evidence_v1/filtered_ucr_snp_info_v1.Rdata")

# filter_ucr_left_snp_info -------------------------------------------------------

# 将chr_name列的数据类型设置为character
ucr_left_snp_info <- lapply(ucr_left_snp_info, function(df) {
  if ("validated" %in% colnames(df)) {
    df$chr_name <- as.character(df$chr_name)
  } else if ("Evidence" %in% colnames(df)) {
    df$chr_name <- as.character(df$chr_name)
  }
  return(df)
})

# 合并数据框
combined_ucr_left_snp_info <- bind_rows(ucr_left_snp_info, .id = "list_type")

# 将所有空白值修改为缺失值
combined_ucr_left_snp_info <- combined_ucr_left_snp_info %>% mutate_all(~ ifelse(. == "", NA, .))

# 筛选validated和Evidence不同时为NA的行
filtered_combined_ucr_left_snp_info <- combined_ucr_left_snp_info %>%
  filter(!is.na(validated) | !is.na(Evidence))

save(filtered_combined_ucr_left_snp_info, file = "filtered_UCR_SNP_with_evidence_v1/filtered_ucr_left_snp_info_v1.Rdata")

# filter_ucr_right_snp_info ------------------------------------------------------

# 将chr_name列的数据类型设置为character
ucr_right_snp_info <- lapply(ucr_right_snp_info, function(df) {
  if ("validated" %in% colnames(df)) {
    df$chr_name <- as.character(df$chr_name)
  } else if ("Evidence" %in% colnames(df)) {
    df$chr_name <- as.character(df$chr_name)
  }
  return(df)
})

# 合并数据框
combined_ucr_right_snp_info <- bind_rows(ucr_right_snp_info, .id = "list_type")

# 将所有空白值修改为缺失值
combined_ucr_right_snp_info <- combined_ucr_right_snp_info %>% mutate_all(~ ifelse(. == "", NA, .))

# 筛选validated和Evidence不同时为NA的行
filtered_combined_ucr_right_snp_info <- combined_ucr_right_snp_info %>%
  filter(!is.na(validated) | !is.na(Evidence))

save(filtered_combined_ucr_right_snp_info, file = "filtered_UCR_SNP_with_evidence_v1/filtered_ucr_right_snp_info_v1.Rdata")

# filter_random_fragments_snp_info ----------------------------------------

# 使用map_lgl函数检查每个数据框是否为空
# 如果为空则返回FALSE，不为空则返回TRUE
keep_data_frames <- map_lgl(random_fragments_snp_info, ~ !all(unlist(.x) %in% c(NA, "")))

# 使用keep函数从列表中保留具有数据的元素
random_fragments_snp_info <- random_fragments_snp_info[keep_data_frames]

# 将chr_name列的数据类型设置为character
random_fragments_snp_info <- lapply(random_fragments_snp_info, function(df) {
  if ("validated" %in% colnames(df)) {
    df$chr_name <- as.character(df$chr_name)
  } else if ("Evidence" %in% colnames(df)) {
    df$chr_name <- as.character(df$chr_name)
  }
  return(df)
})

# 合并数据框
combined_random_fragments_snp_info <- bind_rows(random_fragments_snp_info, .id = "list_type")

# 将所有空白值修改为缺失值
combined_random_fragments_snp_info <- combined_random_fragments_snp_info %>% mutate_all(~ ifelse(. == "", NA, .))

# 筛选validated和Evidence不同时为NA的行
filtered_combined_random_fragments_snp_info <- combined_random_fragments_snp_info %>%
  filter(!is.na(validated) | !is.na(Evidence))

save(filtered_combined_random_fragments_snp_info, file = "filtered_UCR_SNP_with_evidence_v1/filtered_random_fragments_snp_info_v1.Rdata")

# standardize_step_format_snp -------------------------------------------------------

setwd(dir = "D:/R_project/UCR_project/02-analysis/03-SNP_INFO_Retrieval")
rm(list = ls())

library(tidyverse)

load("filtered_UCR_SNP_with_evidence_v1/filtered_ucr_snp_info_v1.Rdata")
load("filtered_UCR_SNP_with_evidence_v1/filtered_ucr_right_snp_info_v1.Rdata")
load("filtered_UCR_SNP_with_evidence_v1/filtered_ucr_left_snp_info_v1.Rdata")
load("filtered_UCR_SNP_with_evidence_v1/filtered_random_fragments_snp_info_v1.Rdata")

## 把getBM和直接从网页上获得的snp信息整理成统一格式
## ucr_snp
temp <- filtered_combined_ucr_snp_info[c(23441:23788), ]
temp$refsnp_id <- temp$Variant.ID
temp$refsnp_source <- temp$Source
temp$refsnp_source_description <- temp$Class
temp$allele <- temp$Alleles
temp$chrom_start <- str_sub(temp$Location, 3, str_length(temp$Location))
temp$validated <- temp$Evidence
temp <- temp[, c(-10:-15)]
filtered_combined_ucr_snp_info <- filtered_combined_ucr_snp_info[, c(-10:-15)]
filtered_combined_ucr_snp_info[c(23441:23788), ] <- temp
save(filtered_combined_ucr_snp_info, file = "filtered_UCR_SNP_with_evidence_v1/filtered_ucr_snp_info_clean_v1.Rdata")

## ucr_left_snp
colSums(is.na(filtered_combined_ucr_left_snp_info))
colSums(is.na(filtered_combined_ucr_left_snp_info[c(27227:27567), ])) # 确定从27227行开始下面全是从网页获得的
temp <- filtered_combined_ucr_left_snp_info[c(27227:27567), ]
temp$refsnp_id <- temp$Variant.ID
temp$refsnp_source <- temp$Source
temp$refsnp_source_description <- temp$Class
temp$allele <- temp$Alleles
temp$chrom_start <- str_sub(temp$Location, 3, str_length(temp$Location))
temp$validated <- temp$Evidence
temp <- temp[, c(-10:-15)]
filtered_combined_ucr_left_snp_info <- filtered_combined_ucr_left_snp_info[, c(-10:-15)]
filtered_combined_ucr_left_snp_info[c(27227:27567), ] <- temp
save(filtered_combined_ucr_left_snp_info, file = "filtered_UCR_SNP_with_evidence_v1/filtered_ucr_left_snp_info_clean_v1.Rdata")

## ucr_right_snp
colSums(is.na(filtered_combined_ucr_right_snp_info))
colSums(is.na(filtered_combined_ucr_left_snp_info[c(27481:27813), ])) # 确定从27481行开始下面全是从网页获得的
temp <- filtered_combined_ucr_right_snp_info %>% 
  filter(is.na(refsnp_id))
temp$refsnp_id <- temp$Variant.ID
temp$refsnp_source <- temp$Source
temp$refsnp_source_description <- temp$Class
temp$allele <- temp$Alleles
temp$chrom_start <- str_sub(temp$Location, 3, str_length(temp$Location))
temp$validated <- temp$Evidence
temp <- temp[, c(-10:-15)]
filtered_combined_ucr_right_snp_info <- filtered_combined_ucr_right_snp_info %>% 
  filter(!is.na(refsnp_id))
filtered_combined_ucr_right_snp_info <- filtered_combined_ucr_right_snp_info[, c(-10:-15)]
filtered_combined_ucr_right_snp_info <- rbind(filtered_combined_ucr_right_snp_info, temp)
save(filtered_combined_ucr_right_snp_info, file = "filtered_UCR_SNP_with_evidence_v1/filtered_ucr_right_snp_info_clean_v1.Rdata")

## random_fragments_snp
colSums(is.na(filtered_combined_random_fragments_snp_info))
colSums(is.na(filtered_combined_random_fragments_snp_info[c(46244:47663), ]))
temp <- filtered_combined_random_fragments_snp_info %>% 
  filter(is.na(refsnp_id))
temp$refsnp_id <- temp$Variant.ID
temp$refsnp_source <- temp$Source
temp$refsnp_source_description <- temp$Class
temp$allele <- temp$Alleles
temp$chrom_start <- str_sub(temp$Location, 3, str_length(temp$Location))
temp$validated <- temp$Evidence
temp$chr_name <- str_sub(temp$Location, 1, 1) # rf当中chr_name和Location当中的染色体数目不一样，Location应该是对的，chr_name当中的X好像是我自己添加的
temp <- temp[, c(-10:-15)]
filtered_combined_random_fragments_snp_info <- filtered_combined_random_fragments_snp_info %>% 
  filter(!is.na(refsnp_id))
filtered_combined_random_fragments_snp_info <- filtered_combined_random_fragments_snp_info[, c(-10:-15)]
filtered_combined_random_fragments_snp_info <- rbind(filtered_combined_random_fragments_snp_info, temp)
save(filtered_combined_random_fragments_snp_info, file = "filtered_UCR_SNP_with_evidence_v1/filtered_rf_snp_info_clean_v1.Rdata")















