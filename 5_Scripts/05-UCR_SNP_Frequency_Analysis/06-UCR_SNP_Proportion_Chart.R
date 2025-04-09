rm(list = ls())
setwd(dir = "D:/R_project/UCR_project")

library(tidyverse)
library(ggplot2)

load("02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/passed_UCR_SNPs.Rdata")
load("02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/passed_UCR_left_SNPs.Rdata")
load("02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/passed_UCR_right_SNPs.Rdata")
load("02-analysis/04-SNP_Quality_Control/filtered_TOPMed_SNPs/passed_rf_SNPs.Rdata")

# UCR_SNP_proportion
UCR_SNP_proportion <- passed_UCR_SNPs %>% 
  count(list_type) %>% 
  rename(value = list_type, count = n) %>% 
  arrange(desc(count))
colnames(UCR_SNP_proportion) <- c("UCR_name", "SNP_num")

UCR_location <- read.table(file = "D:/R_project/UCR_project/data/flanking_snps/UCR_location.txt", 
                           sep = "\t", header = T)
UCR_SNP_proportion <- merge(UCR_SNP_proportion, UCR_location, by = "UCR_name")
UCR_SNP_proportion <- UCR_SNP_proportion[, c(1, 2, 6)]
colnames(UCR_SNP_proportion)[3] <- "UCR_length"
UCR_SNP_proportion <- UCR_SNP_proportion %>% 
  mutate(SNP_proportion = SNP_num/UCR_length)
save(UCR_SNP_proportion, file = "02-analysis/05-UCR_SNP_Frequency_Analysis/UCR_SNP_proportion/UCR_SNP_proportion.Rdata")

# UCR_left_SNP_proportion
UCR_left_SNP_proportion <- passed_UCR_left_SNPs %>% 
  count(list_type) %>% 
  rename(value = list_type, count = n) %>% 
  arrange(desc(count))
colnames(UCR_left_SNP_proportion) <- c("UCR_name", "SNP_num")

UCR_left_SNP_proportion <- merge(UCR_left_SNP_proportion, UCR_location, by = "UCR_name")
UCR_left_SNP_proportion <- UCR_left_SNP_proportion[, c(1, 2, 6)]
colnames(UCR_left_SNP_proportion)[3] <- "UCR_length"
UCR_left_SNP_proportion <- UCR_left_SNP_proportion %>% 
  mutate(SNP_proportion = SNP_num/UCR_length)
save(UCR_left_SNP_proportion, file = "02-analysis/05-UCR_SNP_Frequency_Analysis/UCR_SNP_proportion/UCR_left_SNP_proportion.Rdata")

# UCR_right_SNP_proportion
UCR_right_SNP_proportion <- passed_UCR_right_SNPs %>% 
  count(list_type) %>% 
  rename(value = list_type, count = n) %>% 
  arrange(desc(count))
colnames(UCR_right_SNP_proportion) <- c("UCR_name", "SNP_num")

UCR_right_SNP_proportion <- merge(UCR_right_SNP_proportion, UCR_location, by = "UCR_name")
UCR_right_SNP_proportion <- UCR_right_SNP_proportion[, c(1, 2, 6)]
colnames(UCR_right_SNP_proportion)[3] <- "UCR_length"
UCR_right_SNP_proportion <- UCR_right_SNP_proportion %>% 
  mutate(SNP_proportion = SNP_num/UCR_length)
save(UCR_right_SNP_proportion, file = "02-analysis/05-UCR_SNP_Frequency_Analysis/UCR_SNP_proportion/UCR_right_SNP_proportion.Rdata")

# Random_Fragments_SNP_proportion
RF_SNP_proportion <- passed_rf_SNPs %>% 
  count(list_type) %>% 
  rename(value = list_type, count = n) %>% 
  arrange(desc(count))
colnames(RF_SNP_proportion) <- c("RF_name", "SNP_num")

load("D:/R_project/UCR_project/02-analysis/03-SNP_INFO_Retrieval/flanking_snps/random_fragmens.Rdata")
colnames(selected_segments)[6] <- "RF_name"
RF_SNP_proportion <- merge(RF_SNP_proportion, selected_segments, by = "RF_name")

RF_SNP_proportion <- RF_SNP_proportion[, c(1, 2, 6)]
colnames(RF_SNP_proportion)[3] <- "RF_length"
RF_SNP_proportion <- RF_SNP_proportion %>% 
  mutate(SNP_proportion = SNP_num/RF_length)
save(RF_SNP_proportion, file = "02-analysis/05-UCR_SNP_Frequency_Analysis/UCR_SNP_proportion/RF_SNP_proportion.Rdata")

# 比较各组之间是否存在显著差异
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(car)

load("02-analysis/05-UCR_SNP_Frequency_Analysis/UCR_SNP_proportion/UCR_left_SNP_proportion.Rdata")
load("02-analysis/05-UCR_SNP_Frequency_Analysis/UCR_SNP_proportion/UCR_right_SNP_proportion.Rdata")
load("02-analysis/05-UCR_SNP_Frequency_Analysis/UCR_SNP_proportion/UCR_SNP_proportion.Rdata")
load("02-analysis/05-UCR_SNP_Frequency_Analysis/UCR_SNP_proportion/RF_SNP_proportion.Rdata")

temp <- c(RF_SNP_proportion$SNP_proportion, c(rep("0", 43)))
df <- cbind(Group1 = UCR_SNP_proportion$SNP_proportion, 
            Group2 = UCR_left_SNP_proportion$SNP_proportion, 
            Group3 = UCR_right_SNP_proportion$SNP_proportion, 
            Group4 = temp)
colnames(df) <- c("UCR", "UCR_left", "UCR_right", "Random_fragments")

# 进行一元方差分析
# UCRvsUCR_left:4.68e-05 ***
UCRvsUCR_left <- aov(UCR ~ UCR_left, data = df)

# 显示方差分析结果
summary(anova_result)
#              Df Sum Sq   Mean Sq F value   Pr(>F)    
# UCR_left    434 1.1138 0.0025664   2.697 4.68e-05 ***
# Residuals    46 0.0438 0.0009516                     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# UCRvsUCR_right:*
UCRvsUCR_right <- aov(UCR ~ UCR_right, data = df)
summary(UCRvsUCR_right)
#              Df Sum Sq  Mean Sq F value Pr(>F)  
# UCR_right   447 1.1107 0.002485   1.748 0.0258 *
# Residuals    33 0.0469 0.001421                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# UCRvsRF:*
UCRvsRF <- aov(UCR ~ Random_fragments, data = df)
summary(UCRvsRF)
#                   Df Sum Sq  Mean Sq F value Pr(>F)  
# Random_fragments 428 1.0738 0.002509   1.558 0.0252 *
# Residuals         52 0.0838 0.001611                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# UCR_leftvsUCR_right:0.0668
UCR_leftvsUCR_right <- aov(UCR_left ~ UCR_right, data = df)
summary(UCR_leftvsUCR_right)

write.table(df, file = "D:/R_project/UCR_project/data/UCR_SNP_proportion/UCR_SNP_proportion.txt", 
            row.names = F, quote = F)





