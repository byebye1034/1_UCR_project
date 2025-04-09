# 这里是SNP的第二步QC：筛选所有存在证据支持的SNP放在文件passed_XXX_SNPs.Rdata文件中

rm(list = ls())
setwd(dir = "D:/R_project/UCR_project/02-analysis/04-SNP_Quality_Control/")

library(tidyverse)

load("D:/R_project/UCR_project/02-analysis/03-SNP_INFO_Retrieval/filtered_UCR_SNP_with_evidence_v1/filtered_ucr_snp_info_clean_v1.Rdata")
load("D:/R_project/UCR_project/02-analysis/03-SNP_INFO_Retrieval/filtered_UCR_SNP_with_evidence_v1/filtered_ucr_right_snp_info_clean_v1.Rdata")
load("D:/R_project/UCR_project/02-analysis/03-SNP_INFO_Retrieval/filtered_UCR_SNP_with_evidence_v1/filtered_ucr_left_snp_info_clean_v1.Rdata")
load("D:/R_project/UCR_project/02-analysis/03-SNP_INFO_Retrieval/filtered_UCR_SNP_with_evidence_v1/filtered_rf_snp_info_clean_v1.Rdata")

# UCR SNPs ----------------------------------------------------------------

# 准备所有passed的TOPMed的SNP，用于第三步筛选
snp_file_names <- list.files(path = "D:/R_project/UCR_project/02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/ucr_snp_TOPMed_info")
snp_file_names <- paste0("D:/R_project/UCR_project/02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/ucr_snp_TOPMed_info/", snp_file_names)

UCR_TOPMed_snp <- data.frame()
for (file in snp_file_names){
  temper <- read.table(file = file, sep = "\t", skip = 16)
  UCR_TOPMed_snp <- rbind(UCR_TOPMed_snp, temper)
}
UCR_TOPMed_snp <- UCR_TOPMed_snp[, -8]
colnames(UCR_TOPMed_snp) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")
passed_UCR_TOPMed_snp <- filter(UCR_TOPMed_snp, FILTER == "PASS")
passed_UCR_TOPMed_snp$POS <- as.character(passed_UCR_TOPMed_snp$POS)
passed_UCR_TOPMed_snp <- passed_UCR_TOPMed_snp %>% 
  mutate(LOCATION = paste0(CHROM, sep = ":", POS))

# 第一步：筛选validated中包含gnomAD的行
passed_step1 <- filtered_combined_ucr_snp_info %>%
  filter(grepl("gnomAD", validated, ignore.case = TRUE))

# 第二步：筛选temp当中validated列包含1000Genomes的行
temp <- filtered_combined_ucr_snp_info %>%
  anti_join(passed_step1, by = "refsnp_id")  # 请替换 "your_unique_identifier_column" 为你数据框中的唯一标识列

passed_step2 <- temp %>%
  filter(grepl("1000Genomes", validated, ignore.case = TRUE))

# 第三步：筛选temp当中validated列包含TOPMed的行step3，然后筛选其中pass的部分passed_step3
temp <- temp %>%
  anti_join(passed_step2, by = "refsnp_id")

step3 <- temp %>%
  filter(grepl("TOPMed", validated, ignore.case = TRUE)) %>% 
  mutate(location = paste0(chr_name, sep = ":", chrom_start))

passed_step3 <- step3 %>% 
  filter(step3$location %in% passed_UCR_TOPMed_snp$LOCATION)
passed_step3 <- passed_step3[, -10]

# 第四步：具有Frequency、Phenotype_or_Disease和cited证据的也都保留
temp <- temp %>% 
  anti_join(step3, by = "refsnp_id")

passed_step4 <- temp %>% 
  filter(grepl("Frequency|Phenotype_or_Disease|cited", validated, ignore.case = TRUE))

passed_UCR_SNPs <- rbind(passed_step1, passed_step2, passed_step3, passed_step4)
save(passed_UCR_SNPs, file = "filtered_TOPMed_SNPs/passed_UCR_SNPs.Rdata")

# UCR left SNPs -----------------------------------------------------------

# 准备所有passed的TOPMed的SNP，用于第三步筛选
snp_file_names <- list.files(path = "D:/R_project/UCR_project/02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/ucr_left_snp_info")
snp_file_names <- paste0("D:/R_project/UCR_project/02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/ucr_left_snp_info/", snp_file_names)

UCR_left_TOPMed_snp <- data.frame()
for (file in snp_file_names){
  temper <- read.table(file = file, sep = "\t", skip = 16)
  UCR_left_TOPMed_snp <- rbind(UCR_left_TOPMed_snp, temper)
}
UCR_left_TOPMed_snp <- UCR_left_TOPMed_snp[, -8]
colnames(UCR_left_TOPMed_snp) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")
passed_UCR_left_TOPMed_snp <- filter(UCR_left_TOPMed_snp, FILTER == "PASS")
passed_UCR_left_TOPMed_snp$POS <- as.character(passed_UCR_left_TOPMed_snp$POS)
passed_UCR_left_TOPMed_snp <- passed_UCR_left_TOPMed_snp %>% 
  mutate(LOCATION = paste0(CHROM, sep = ":", POS))

# 第一步：筛选validated中包含gnomAD的行
passed_step1 <- filtered_combined_ucr_left_snp_info %>%
  filter(grepl("gnomAD", validated, ignore.case = TRUE))

# 第二步：筛选temp当中validated列包含1000Genomes的行
temp <- filtered_combined_ucr_left_snp_info %>%
  anti_join(passed_step1, by = "refsnp_id")  # 请替换 "your_unique_identifier_column" 为你数据框中的唯一标识列

passed_step2 <- temp %>%
  filter(grepl("1000Genomes", validated, ignore.case = TRUE))

# 第三步：筛选temp当中validated列包含TOPMed的行step3，然后筛选其中pass的部分passed_step3
temp <- temp %>%
  anti_join(passed_step2, by = "refsnp_id")

step3 <- temp %>%
  filter(grepl("TOPMed", validated, ignore.case = TRUE)) %>% 
  mutate(location = paste0(chr_name, sep = ":", chrom_start))

passed_step3 <- step3 %>% 
  filter(step3$location %in% passed_UCR_left_TOPMed_snp$LOCATION)
passed_step3 <- passed_step3[, -10]

# 第四步：具有Frequency、Phenotype_or_Disease和cited证据的也都保留
temp <- temp %>% 
  anti_join(step3, by = "refsnp_id")

passed_step4 <- temp %>% 
  filter(grepl("Frequency|Phenotype_or_Disease|cited", validated, ignore.case = TRUE))

passed_UCR_left_SNPs <- rbind(passed_step1, passed_step2, passed_step3, passed_step4)
save(passed_UCR_left_SNPs, file = "filtered_TOPMed_SNPs/passed_UCR_left_SNPs.Rdata")

# UCR right SNPs ----------------------------------------------------------

rm(list = ls())
library(tidyverse)

# 准备所有passed的TOPMed的SNP，用于第三步筛选
snp_file_names <- list.files(path = "filtered_TOPMed_SNPs/ucr_right_snp_info")
snp_file_names <- paste0("D:/R_project/UCR_project/02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/ucr_right_snp_info/", snp_file_names)

UCR_right_TOPMed_snp <- data.frame()
for (file in snp_file_names){
  temper <- read.table(file = file, sep = "\t", skip = 16)
  UCR_right_TOPMed_snp <- rbind(UCR_right_TOPMed_snp, temper)
}
UCR_right_TOPMed_snp <- UCR_right_TOPMed_snp[, -8]
colnames(UCR_right_TOPMed_snp) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")
passed_UCR_right_TOPMed_snp <- filter(UCR_right_TOPMed_snp, FILTER == "PASS")
passed_UCR_right_TOPMed_snp$POS <- as.character(passed_UCR_right_TOPMed_snp$POS)
passed_UCR_right_TOPMed_snp <- passed_UCR_right_TOPMed_snp %>% 
  mutate(LOCATION = paste0(CHROM, sep = ":", POS))

# 第一步：筛选validated中包含gnomAD的行
passed_step1 <- filtered_combined_ucr_right_snp_info %>%
  filter(grepl("gnomAD", validated, ignore.case = TRUE))

# 第二步：筛选temp当中validated列包含1000Genomes的行
temp <- filtered_combined_ucr_right_snp_info %>%
  anti_join(passed_step1, by = "refsnp_id")  # 请替换 "your_unique_identifier_column" 为你数据框中的唯一标识列

passed_step2 <- temp %>%
  filter(grepl("1000Genomes", validated, ignore.case = TRUE))

# 第三步：筛选temp当中validated列包含TOPMed的行step3，然后筛选其中pass的部分passed_step3
temp <- temp %>%
  anti_join(passed_step2, by = "refsnp_id")

step3 <- temp %>%
  filter(grepl("TOPMed", validated, ignore.case = TRUE)) %>% 
  mutate(location = paste0(chr_name, sep = ":", chrom_start))

passed_step3 <- step3 %>% 
  filter(step3$location %in% passed_UCR_right_TOPMed_snp$LOCATION)
passed_step3 <- passed_step3[, -10]

# 第四步：具有Frequency、Phenotype_or_Disease和cited证据的也都保留
temp <- temp %>% 
  anti_join(step3, by = "refsnp_id")

passed_step4 <- temp %>% 
  filter(grepl("Frequency|Phenotype_or_Disease|cited", validated, ignore.case = TRUE))

passed_UCR_right_SNPs <- rbind(passed_step1, passed_step2, passed_step3, passed_step4)
save(passed_UCR_right_SNPs, file = "filtered_TOPMed_SNPs/passed_UCR_right_SNPs.Rdata")

# Random Fragments SNPs ---------------------------------------------------

rm(list = ls())

library(tidyverse)

# 准备所有passed的TOPMed的SNP，用于第三步筛选
snp_file_names <- list.files(path = "filtered_TOPMed_SNPs/rf_snp_info/")
snp_file_names <- paste0("D:/R_project/UCR_project/02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/rf_snp_info/", snp_file_names)

large_files <- character(0)
# 遍历文件列表
for (file in snp_file_names) {
  # 获取文件的大小信息
  file_info <- file.info(file)
  
  # 如果文件大小小于1.5k（1500字节），将文件名添加到small_files向量中
  if (file_info$size > 1500) {
    large_files <- c(large_files, file)
  }
}

rf_TOPMed_snp <- data.frame()
for (file in large_files){
  temper <- read.table(file = file, sep = "\t", skip = 16)
  rf_TOPMed_snp <- rbind(rf_TOPMed_snp, temper)
}
rf_TOPMed_snp <- rf_TOPMed_snp[, -8]
colnames(rf_TOPMed_snp) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")
passed_rf_TOPMed_snp <- filter(rf_TOPMed_snp, FILTER == "PASS")
passed_rf_TOPMed_snp$POS <- as.character(passed_rf_TOPMed_snp$POS)
passed_rf_TOPMed_snp <- passed_rf_TOPMed_snp %>% 
  mutate(LOCATION = paste0(CHROM, sep = ":", POS))

# 第一步：筛选validated中包含gnomAD的行
passed_step1 <- filtered_combined_random_fragments_snp_info %>%
  filter(grepl("gnomAD", validated, ignore.case = TRUE))

# 第二步：筛选temp当中validated列包含1000Genomes的行
temp <- filtered_combined_random_fragments_snp_info %>%
  anti_join(passed_step1, by = "refsnp_id")  # 请替换 "your_unique_identifier_column" 为你数据框中的唯一标识列

passed_step2 <- temp %>%
  filter(grepl("1000Genomes", validated, ignore.case = TRUE))

# 第三步：筛选temp当中validated列包含TOPMed的行step3，然后筛选其中pass的部分passed_step3
temp <- temp %>%
  anti_join(passed_step2, by = "refsnp_id")

step3 <- temp %>%
  filter(grepl("TOPMed", validated, ignore.case = TRUE)) %>% 
  mutate(location = paste0(chr_name, sep = ":", chrom_start))

passed_step3 <- step3 %>% 
  filter(step3$location %in% passed_rf_TOPMed_snp$LOCATION)
passed_step3 <- passed_step3[, -10]

# 第四步：具有Frequency、Phenotype_or_Disease和cited证据的也都保留
temp <- temp %>% 
  anti_join(step3, by = "refsnp_id")

passed_step4 <- temp %>% 
  filter(grepl("Frequency|Phenotype_or_Disease|cited", validated, ignore.case = TRUE))

passed_rf_SNPs <- rbind(passed_step1, passed_step2, passed_step3, passed_step4)
save(passed_rf_SNPs, file = "filtered_TOPMed_SNPs/passed_rf_SNPs.Rdata")
















































































