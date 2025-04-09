# construct ANN model

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library(neuralnet)
library(pROC)
library(caret)
library(readxl)
library(ggview)
library(data.table)

lwd_pt <- .pt*72.27/96

set.seed(123)
n <- 200
data <- data.frame(
  gene1 = rnorm(n),
  gene2 = rnorm(n),
  gene3 = rnorm(n),
  diagnosis = factor(sample(0:1, n, replace = TRUE))
)

# 数据预处理：将诊断标签转换为数值型
data$diagnosis <- as.numeric(data$diagnosis) - 1

# 划分训练集和验证集
trainIndex <- createDataPartition(data$diagnosis, p = 0.7, list = FALSE)
trainData <- data[trainIndex, ]
validData <- data[-trainIndex, ]

# 构建公式
formula <- as.formula(paste("diagnosis ~", paste(names(trainData)[!names(trainData) %in% "diagnosis"], collapse = " + ")))

# 构建人工神经网络模型
set.seed(123)
nnModel <- neuralnet(formula, 
                     data = trainData, 
                     hidden = c(3),  # 隐藏层节点数
                     linear.output = FALSE)

# 在验证集上进行预测
predictions <- compute(nnModel, validData[, !names(validData) %in% "diagnosis"])
predictedProbs <- predictions$net.result

# 绘制AUC曲线
rocCurve <- roc(validData$diagnosis, predictedProbs)
plot(rocCurve, main = "AUC Curve for LGG Diagnosis Model", col = "blue")
aucValue <- auc(rocCurve)
legend("bottomright", legend = paste("AUC =", round(aucValue, 3)), lty = 1, col = "blue")

# GBM ---------------------------------------------------------------------

## import GSE147352 reads count
load("D:/R_project/UCR_project/01-data/45-GBM/GSE147352_ReadsCount.rdata")

# get sample information
SampleInfo <- read_xlsx(path = "01-data/45-GBM/GSE147352_clininfo.xlsx", sheet = 1)
table(SampleInfo$Type)
# GBM 85 LGG 18 Normal 15
GBMSampleInfo <- filter(SampleInfo, Type == "Normal" | Type == "GBM")
LGGSampleInfo <- filter(SampleInfo, Type == "Normal" | Type == "LGG")

## import genes
GBMGene <- c("ENSG00000224243", "ENSG00000257545")

GBMdata <- ReadsCount[GBMGene, ]
GBMdata <- select(GBMdata, GBMSampleInfo$`sample ID`)
GBMdata <- as.data.frame(t(GBMdata))
GBMdata$sample_ID <- rownames(GBMdata)
colnames(GBMSampleInfo)[1] <- "sample_ID"

GBMdata <- merge(GBMdata, GBMSampleInfo[, c(1:2)], by = "sample_ID")
GBMdata$Type <- gsub("GBM", "1", GBMdata$Type)
GBMdata$Type <- gsub("Normal", "0", GBMdata$Type)

# 数据预处理：将诊断标签转换为数值型
GBMdata$Type <- as.numeric(GBMdata$Type)
GBMdata <- GBMdata[, c(2:4)]

GBMTrainData <- GBMdata
trainData <- GBMTrainData

# 
TCGA_GTEx_GBM_validation_set$Type <- gsub("GBM", "1", TCGA_GTEx_GBM_validation_set$Type)
TCGA_GTEx_GBM_validation_set$Type <- gsub("Normal", "0", TCGA_GTEx_GBM_validation_set$Type)
TCGA_GTEx_GBM_validation_set$Type <- as.numeric(TCGA_GTEx_GBM_validation_set$Type)

GBMValidData <- TCGA_GTEx_GBM_validation_set
for (i in 1:2) {
  GBMValidData[[i]] <- as.numeric(GBMValidData[[i]])
}
GBMValidData[, 1:2] <- 2^GBMValidData[, 1:2] - 1
validData <- GBMValidData
validData <- as.numeric(validData)

save(trainData, validData, 
     file = "03-results/52-ANN_Model/GBM_ANN.rdata")

# # 划分训练集和验证集
# trainIndex <- createDataPartition(GBMdata$Type, p = 0.7, list = FALSE)
# trainData <- GBMdata[trainIndex, ]
# validData <- GBMdata[-trainIndex, ]

# 构建公式
formula <- as.formula(paste("Type ~", paste(names(trainData)[!names(trainData) %in% "Type"], collapse = " + ")))

# 构建人工神经网络模型
set.seed(123)
nnModel <- neuralnet(formula, 
                     data = trainData, 
                     hidden = c(4),  # 隐藏层节点数
                     linear.output = FALSE)

# 在验证集上进行预测
predictions <- compute(nnModel, validData[, !names(validData) %in% "Type"])
predictedProbs <- predictions$net.result

# 绘制AUC曲线
rocCurve <- roc(validData$Type, predictedProbs)
# plot(rocCurve, main = "AUC Curve for GBM Diagnosis Model", col = "blue")
# aucValue <- auc(rocCurve)
# legend("bottomright", legend = paste("AUC =", round(aucValue, 3)), lty = 1, col = "blue")

# 将 ROC 曲线数据转换为数据框
rocData <- data.frame(
  FPR = rocCurve$specificities,
  TPR = rocCurve$sensitivities
)

# 使用 ggplot2 绘制 AUC 曲线
plotAUCCurve <- function(data, tumor){
  ggplot(data = data, mapping = aes(x = 1-FPR, y = TPR)) +
    geom_line(color = "red") +
    geom_abline(intercept = 0, slope = 1, 
                linetype = "dashed", color = "gray") +
    
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    
    labs(title = paste0("AUC Curve for ", tumor, " Diagnosis Model"),
         x = "False Positive Rate",
         y = "True Positive Rate") +
    
    annotate("text", x = 0.7, y = 0.2, 
             label = paste("AUC =", round(auc(rocCurve), 3)), 
             color = "red", size = 7/lwd_pt) +
    
    theme(
      panel.grid = element_blank(), 
      panel.background = element_blank(), 
      text = element_text(size = 7, color = "#000000"), 
      line = element_line(linewidth = 0.5/lwd_pt, color = "#000000"), 
      
      axis.title = element_text(size = 7, color = "#000000"), 
      axis.text = element_text(size = 7, color = "#000000"), 
      axis.ticks = element_line(linewidth = 0.5/lwd_pt, color = "#000000"), 
      axis.line = element_line(linewidth = 0.5/lwd_pt, color = "#000000"),
      
      legend.position = "top", 
      legend.key.size = unit(0.25, "cm"), 
      legend.title = element_text(size = 7, color = "#000000"), 
      legend.text = element_text(size = 7, color = "#000000"), 
      
      plot.title = element_text(size = 7, color = "#000000"),
      
      aspect.ratio = 1:1
    )
}

plotAUCCurve(data = rocData, tumor = "GBM")

pdf(file = "03-results/52-ANN_Model/AUC_GBM.pdf", width = 480/254, height = 480/254)
plotAUCCurve(data = rocData, tumor = "GBM")
dev.off()

# LGG ---------------------------------------------------------------------

# get sample information
SampleInfo <- read_xlsx(path = "01-data/45-GBM/GSE147352_clininfo.xlsx", sheet = 1)
table(SampleInfo$Type)
# GBM 85 LGG 18 Normal 15
LGGSampleInfo <- filter(SampleInfo, Type == "Normal" | Type == "LGG")

## import genes
LGGGene <- c("ENSG00000204929", "ENSG00000224592", "ENSG00000231764", 
             "ENSG00000234323", "ENSG00000257585", "ENSG00000257986")

LGGdata <- ReadsCount[LGGGene, ]
LGGdata <- select(LGGdata, LGGSampleInfo$`sample ID`)
LGGdata <- as.data.frame(t(LGGdata))
LGGdata$sample_ID <- rownames(LGGdata)
colnames(LGGSampleInfo)[1] <- "sample_ID"

LGGdata <- merge(LGGdata, LGGSampleInfo[, c(1:2)], by = "sample_ID")
LGGdata$Type <- gsub("LGG", "1", LGGdata$Type)
LGGdata$Type <- gsub("Normal", "0", LGGdata$Type)

# 数据预处理：将诊断标签转换为数值型
LGGdata$Type <- as.numeric(LGGdata$Type)
LGGdata <- LGGdata[, c(2:8)]

LGGTrainData <- LGGdata
trainData <- LGGTrainData

TCGA_GTEx_LGG_validation_set$Type <- gsub("LGG", "1", TCGA_GTEx_LGG_validation_set$Type)
TCGA_GTEx_LGG_validation_set$Type <- gsub("Normal", "0", TCGA_GTEx_LGG_validation_set$Type)
TCGA_GTEx_LGG_validation_set$Type <- as.numeric(TCGA_GTEx_LGG_validation_set$Type)

LGGValidData <- TCGA_GTEx_LGG_validation_set
# 确保前六列是数值型
for (i in 1:6) {
  LGGValidData[[i]] <- as.numeric(LGGValidData[[i]])
}
LGGValidData[, 1:6] <- 2^LGGValidData[, 1:6] - 1
validData <- LGGValidData
validData <- as.numeric(validData)

# # 划分训练集和验证集
# trainIndex <- createDataPartition(LGGdata$Type, p = 0.6, list = FALSE)
# trainData <- LGGdata[trainIndex, ]
# validData <- LGGdata[-trainIndex, ]

# 构建公式
formula <- as.formula(paste("Type ~", paste(names(trainData)[!names(trainData) %in% "Type"], collapse = " + ")))

# 构建人工神经网络模型
set.seed(123)


pdf(file = "03-results/52-ANN_Model/AUC_LGG.pdf", width = 480/254, height = 480/254)

nnModel <- neuralnet(formula, 
                     data = trainData, 
                     hidden = c(4),  # 隐藏层节点数
                     linear.output = FALSE)

# 在验证集上进行预测
predictions <- compute(nnModel, validData[, !names(validData) %in% "Type"])
predictedProbs <- predictions$net.result

# plot AUC curve
rocCurve <- roc(validData$Type, predictedProbs)
# plot(rocCurve, main = "AUC Curve for GBM Diagnosis Model", col = "blue")
# aucValue <- auc(rocCurve)
# legend("bottomright", legend = paste("AUC =", round(aucValue, 3)), lty = 1, col = "blue")

# ROC curve data to data.frame
rocData <- data.frame(
  FPR = rocCurve$specificities,
  TPR = rocCurve$sensitivities
)

plotAUCCurve(data = rocData, tumor = "LGG")

dev.off()

plotAUCCurve(data = rocData, tumor = "LGG") + canvas(width = 4.8, height = 4.8, units = "cm", dpi = 1200)

# TCGA-GTEx data ----------------------------------------------------------

## GTEx sample info
GTEx_phenotype <- read.table(
  file = "01-data/52-ANN_Model/GTEx/GTEX_phenotype", sep = "\t", header = TRUE
)
Brain_tissue <- GTEx_phenotype %>% 
  filter(X_primary_site == "Brain")

## GTEx exp
GTEx_exp <- fread(
  file = "01-data/52-ANN_Model/GTEx/gtex_gene_expected_count.gz", 
  sep = "\t", header = TRUE
)

LGGGene <- c("ENSG00000204929", "ENSG00000224592", "ENSG00000231764", 
             "ENSG00000234323", "ENSG00000257585", "ENSG00000257986")
GBMGene <- c("ENSG00000224243", "ENSG00000257545")
Genes <- c(LGGGene, GBMGene)

GTEx_exp$sample <- str_sub(GTEx_exp$sample, 1, 15)
GTEx_exp[1, 1]
GTEx_exp <- GTEx_exp %>% 
  filter(sample %in% Genes)
save(GTEx_exp, file = "02-analysis/52-ANN_Model/GTEx_exp.rdata")
load(file = "02-analysis/52-ANN_Model/GTEx_exp.rdata")

GTEx_exp <- as.data.frame(GTEx_exp)

Brain_tissue_id <- Brain_tissue$Sample
Brain_tissue_id <- intersect(Brain_tissue_id, colnames(GTEx_exp))
Brain_tissue_id <- c("sample", Brain_tissue_id)

GTEx_brain_exp <- GTEx_exp[, Brain_tissue_id]

rownames(GTEx_brain_exp) <- GTEx_brain_exp$sample
GTEx_brain_exp <- t(GTEx_brain_exp)
GTEx_brain_exp <- GTEx_brain_exp[-1, ]
GTEx_brain_exp <- as.data.frame(GTEx_brain_exp)
GTEx_brain_exp <- mutate(GTEx_brain_exp, Type = "Normal")
save(GTEx_brain_exp, file = "02-analysis/52-ANN_Model/Validation_set_GTEx.rdata")

## TCGA LGG exp
TCGA_LGG_exp <- read.table(
  file = "01-data/52-ANN_Model/TCGA-LGG/TCGA-LGG.star_counts.tsv.gz", 
  sep = "\t", header = TRUE
)

TCGA_LGG_sampleinfo <- read.table(
  file = "01-data/52-ANN_Model/TCGA-LGG/TCGA-LGG.clinical.tsv.gz", 
  sep = "\t", header = TRUE
)

TCGA_LGG_exp$Ensembl_ID <- str_sub(TCGA_LGG_exp$Ensembl_ID, 1, 15)
TCGA_LGG_exp <- TCGA_LGG_exp[TCGA_LGG_exp$Ensembl_ID %in% Genes, ]
TCGA_LGG_exp <- as.data.frame(TCGA_LGG_exp)
TCGA_LGG_exp <- t(TCGA_LGG_exp)
colnames(TCGA_LGG_exp) <- TCGA_LGG_exp[1, ]
TCGA_LGG_exp <- TCGA_LGG_exp[-1, ]
TCGA_LGG_exp <- as.data.frame(TCGA_LGG_exp)
TCGA_LGG_exp <- mutate(TCGA_LGG_exp, Type = "LGG")
save(TCGA_LGG_exp, file = "02-analysis/52-ANN_Model/Validation_set_TCGA_LGG.rdata")

GTEx_brain_exp <- GTEx_brain_exp %>% 
  select("ENSG00000204929", "ENSG00000224592", "ENSG00000231764", 
         "ENSG00000234323", "ENSG00000257585", "ENSG00000257986", "Type")
TCGA_LGG_exp <- TCGA_LGG_exp %>% 
  select("ENSG00000204929", "ENSG00000224592", "ENSG00000231764", 
         "ENSG00000234323", "ENSG00000257585", "ENSG00000257986", "Type")
TCGA_GTEx_LGG_validation_set <- rbind(GTEx_brain_exp, TCGA_LGG_exp)
save(TCGA_GTEx_LGG_validation_set, 
     file = "02-analysis/52-ANN_Model/TCGA_GTEx_LGG_validation_set.rdata")

## TCGA GBM exp
TCGA_GBM_exp <- read.table(
  file = "01-data/52-ANN_Model/TCGA-GBM/TCGA-GBM.star_counts.tsv.gz", 
  sep = "\t", header = TRUE
)

TCGA_GBM_sampleinfo <- read.table(
  file = "01-data/52-ANN_Model/TCGA-GBM/TCGA-GBM.clinical.tsv.gz", 
  sep = "\t", header = TRUE
)
TCGA_GBM_sampleinfo <- TCGA_GBM_sampleinfo %>% 
  mutate(Type = str_sub(sample, nchar(sample)-3, nchar(sample)))
table(TCGA_GBM_sampleinfo$Type)

TCGA_GBM_exp$Ensembl_ID <- str_sub(TCGA_GBM_exp$Ensembl_ID, 1, 15)
TCGA_GBM_exp <- TCGA_GBM_exp[TCGA_GBM_exp$Ensembl_ID %in% Genes, ]
TCGA_GBM_exp <- as.data.frame(TCGA_GBM_exp)
TCGA_GBM_exp <- t(TCGA_GBM_exp)
colnames(TCGA_GBM_exp) <- TCGA_GBM_exp[1, ]
TCGA_GBM_exp <- TCGA_GBM_exp[-1, ]
TCGA_GBM_exp <- as.data.frame(TCGA_GBM_exp)
TCGA_GBM_exp$sample <- rownames(TCGA_GBM_exp)
TCGA_GBM_exp <- TCGA_GBM_exp %>% 
  mutate(Type = str_sub(sample, nchar(sample)-2, nchar(sample)))
table(TCGA_GBM_exp$Type)
TCGA_GBM_exp <- TCGA_GBM_exp %>% 
  mutate(Type1 = if_else(Type == "11A", "Normal", "GBM"))
table(TCGA_GBM_exp$Type1)
TCGA_GBM_exp <- TCGA_GBM_exp[, c(1:8, 11)]
colnames(TCGA_GBM_exp)[9] <- "Type"
save(TCGA_GBM_exp, file = "02-analysis/52-ANN_Model/Validation_set_TCGA_GBM.rdata")

GTEx_brain_exp <- GTEx_brain_exp %>% 
  select("ENSG00000224243", "ENSG00000257545", "Type")
TCGA_GBM_exp <- TCGA_GBM_exp %>% 
  select("ENSG00000224243", "ENSG00000257545", "Type")
TCGA_GTEx_GBM_validation_set <- rbind(GTEx_brain_exp, TCGA_GBM_exp)
save(TCGA_GTEx_GBM_validation_set, 
     file = "02-analysis/52-ANN_Model/TCGA_GTEx_GBM_validation_set.rdata")





