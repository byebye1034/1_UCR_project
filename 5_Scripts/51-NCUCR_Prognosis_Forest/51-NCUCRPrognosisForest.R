# ncUCR overlapping lncRNA prognosis in glioma
# check the association of differential exp ncUCR genes and survival

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library(forestplot)
library(ezcox)
library(grid)

# df <- read.delim(
#   file = "01-data/51-NCUCR_Prognosis_Forest/LINC01019_SATB1-AS1_LGG.tsv"
# )[-2]
# df2 <- df[c(1, 2)]
# surveData <- read.delim(file = "01-data/51-NCUCR_Prognosis_Forest/TCGA-LGG.survival.tsv")
# df2 <- df2 %>% 
#   dplyr::inner_join(surveData, by = "sample")
# colnames(df2)[2] <- "values"
# 
# unicox_res_genes <- ezcox::ezcox(
#   df2 %>% 
#     select(all_of(c("values", paste0(measure, ".time"), measure))) %>% 
#     na.omit(), 
#   covariates = "values", 
#   time = paste0(measure, ".time"), 
#   status = measure, 
#   verbose = FALSE
# )

# TCGA LGG ----------------------------------------------------------------

# import the DEG in LGG
LGGncUCRDEG <- read.table(file = "03-results/45-GBM/GSE147352_LGG_DEG.txt")
LGGncUCRDEG <- LGGncUCRDEG[LGGncUCRDEG$type != "nosig", ]
LGGncUCRDEG$gene <- rownames(LGGncUCRDEG)

# import TCGA exp and survival data
TCGALGGTPM <- read.table(file = "01-data/51-NCUCR_Prognosis_Forest/TCGA-LGG.star_tpm.tsv", header = TRUE)
TCGALGGTPM$Ensembl_ID <- gsub("\\.\\d+", "", TCGALGGTPM$Ensembl_ID)
TCGALGGSurve <- read.table(file = "01-data/51-NCUCR_Prognosis_Forest/TCGA-LGG.survival.tsv", header = TRUE)

TCGALGGTPM <- TCGALGGTPM %>% 
  filter(Ensembl_ID %in% rownames(LGGncUCRDEG))

TCGALGGTPM <- column_to_rownames(TCGALGGTPM, var = "Ensembl_ID")
TCGALGGTPM <- as.data.frame(t(TCGALGGTPM))
TCGALGGTPM <- 2^TCGALGGTPM - 1
TCGALGGTPM$sample <- rownames(TCGALGGTPM)
TCGALGGTPM$sample <- gsub("\\.", "-", TCGALGGTPM$sample)

LGGdata <- TCGALGGSurve %>% 
  dplyr::inner_join(TCGALGGTPM, by = "sample")
save(LGGdata, file = "02-analysis/51-NCUCR_Prognosis_Forest/LGG_temp_data.rdata")

# cox
allRssults <- list()

for(gene in LGGncUCRDEG$gene){
  # 使用 ezcox 进行单因素分析
  unicox_res_genes <- ezcox::ezcox(
    LGGdata %>% 
      select(all_of(c(gene, "OS.time", "OS"))) %>% 
      na.omit(), 
    covariates = gene, 
    time = "OS.time", 
    status = "OS", 
    verbose = FALSE
  )
  
  # 将当前基因的分析结果添加到列表中
  allRssults[[gene]] <- unicox_res_genes
}

# 将列表中的所有结果按行合并
final_result <- do.call(rbind, allRssults)

# 查看最终结果
print(final_result)

final_result <- final_result %>% 
  mutate(HR_log = log(.data$HR)) %>% 
  mutate(lower_95_log = log(.data$lower_95)) %>% 
  mutate(upper_95_log = log(.data$upper_95)) %>% 
  mutate(type = if_else(.data$p.value < 0.05 & .data$HR_log > 0, "Risky", if_else(.data$p.value < 0.05 & .data$HR_log < 0, "Protective", "NS"))) %>% 
  mutate(type = factor(type, levels = c("NS", "Risky", "Protective")))

## visualization
final_result <- subset(final_result, final_result$type != "NS")
final_result$p.value <- if_else(
  final_result$p.value < 0.01, "<0.01", 
  format(round(as.numeric(final_result$p.value), 3), nsmall = 2)
)
final_result$HR_log <- round(final_result$HR_log, 2)
final_result$lower_95_log <- round(final_result$lower_95_log, 2)
final_result$upper_95_log <- round(final_result$upper_95_log, 2)
final_result$HazardRadio <- paste0(final_result$HR_log, "(", final_result$lower_95_log, "-", final_result$upper_95_log, ")")
final_result_LGG <- final_result

# TCGA GBM ----------------------------------------------------------------

# import the DEG in GBM
GBMncUCRDEG <- read.table(file = "03-results/45-GBM/GSE147352_GBM_DEG.txt")
GBMncUCRDEG <- GBMncUCRDEG[GBMncUCRDEG$type != "nosig", ]
GBMncUCRDEG$gene <- rownames(GBMncUCRDEG)

# import TCGA exp and survival data
TCGAGBMTPM <- read.table(file = "01-data/51-NCUCR_Prognosis_Forest/TCGA-GBM.star_tpm.tsv", header = TRUE)
TCGAGBMTPM$Ensembl_ID <- gsub("\\.\\d+", "", TCGAGBMTPM$Ensembl_ID)
TCGAGBMSurve <- read.table(file = "01-data/51-NCUCR_Prognosis_Forest/TCGA-GBM.survival.tsv", header = TRUE)

TCGAGBMTPM <- TCGAGBMTPM %>% 
  filter(Ensembl_ID %in% rownames(GBMncUCRDEG))

TCGAGBMTPM <- column_to_rownames(TCGAGBMTPM, var = "Ensembl_ID")
TCGAGBMTPM <- as.data.frame(t(TCGAGBMTPM))
TCGAGBMTPM <- 2^TCGAGBMTPM - 1
TCGAGBMTPM$sample <- rownames(TCGAGBMTPM)
TCGAGBMTPM$sample <- gsub("\\.", "-", TCGAGBMTPM$sample)

GBMdata <- TCGAGBMSurve %>% 
  dplyr::inner_join(TCGAGBMTPM, by = "sample")
save(GBMdata, file = "02-analysis/51-NCUCR_Prognosis_Forest/GBM_temp_data.rdata")

# cox
allRssults <- list()

for(gene in GBMncUCRDEG$gene){
  # 使用 ezcox 进行单因素分析
  unicox_res_genes <- ezcox::ezcox(
    GBMdata %>% 
      select(all_of(c(gene, "OS.time", "OS"))) %>% 
      na.omit(), 
    covariates = gene, 
    time = "OS.time", 
    status = "OS", 
    verbose = FALSE
  )
  
  # 将当前基因的分析结果添加到列表中
  allRssults[[gene]] <- unicox_res_genes
}

# 将列表中的所有结果按行合并
final_result <- do.call(rbind, allRssults)

# 查看最终结果
print(final_result)

final_result <- final_result %>% 
  mutate(HR_log = log(.data$HR)) %>% 
  mutate(lower_95_log = log(.data$lower_95)) %>% 
  mutate(upper_95_log = log(.data$upper_95)) %>% 
  mutate(type = if_else(.data$p.value < 0.05 & .data$HR_log > 0, "Risky", if_else(.data$p.value < 0.05 & .data$HR_log < 0, "Protective", "NS"))) %>% 
  mutate(type = factor(type, levels = c("NS", "Risky", "Protective")))

## visualization
final_result <- subset(final_result, final_result$type != "NS")
final_result$p.value <- if_else(
  final_result$p.value < 0.01, "<0.01", 
  format(round(as.numeric(final_result$p.value), 3), nsmall = 2)
)
final_result$HR_log <- round(final_result$HR_log, 2)
final_result$lower_95_log <- round(final_result$lower_95_log, 2)
final_result$upper_95_log <- round(final_result$upper_95_log, 2)
final_result$HazardRadio <- paste0(final_result$HR_log, "(", final_result$lower_95_log, "-", final_result$upper_95_log, ")")
final_result_GBM <- final_result

## merge to 1 plot
final_result_glioma <- rbind(final_result_LGG, final_result_GBM)

label <- as.matrix(final_result_glioma[, c(1, 11, 17)])

txt <- fpTxtGp(
  label = gpar(fontsize = 7),  # 标签文字大小
  ticks = gpar(fontsize = 7),  # 刻度文字大小
  xlab = gpar(fontsize = 7),   # x轴标签文字大小
  title = gpar(fontsize = 7)   # 标题文字大小
)

forestplot(
  labeltext = label, 
  mean = final_result_glioma$HR_log, 
  lower = final_result_glioma$lower_95_log, 
  upper = final_result_glioma$upper_95_log, 
  hrzl_lines = TRUE, 
  size = 0.25, 
  ci_column = 2, 
  txt_gp = txt
)

pdf(file = "03-results/51-NCUCR_Prognosis_Forest/forestplot.pdf", 
    width = 960/254, height = 600/254)

dev.off()



