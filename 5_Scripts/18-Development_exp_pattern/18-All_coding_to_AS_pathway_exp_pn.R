# 在富集分析之后，查看all coding genes在所有器官当中的pattern
# 发现在大脑当中AS pathway的基因有specific pattern（这里我想绘制一个热图，不展示基因名，因为基因太多）
# 在mouse、rat的brain当中查看AS pathway的基因的pattern

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(pheatmap)
library(stringr)
library(readxl)

# get human rpkm
human_rpkm <- read.table(file = "01-data/18-Development_exp_pattern/Human_rpkm.txt", 
                         sep = " ", header = TRUE)
colnames(human_rpkm)

cerebellum_human_rpkm <- human_rpkm[, c(1, 55:112)]
heart_human_rpkm <- human_rpkm[, c(1, 113:156)]
kidney_human_rpkm <- human_rpkm[, c(1, 157:192)]
liver_human_rpkm <- human_rpkm[, c(1, 193:241)]
overy_human_rpkm <- human_rpkm[, c(1, 242:259)]
testis_human_rpkm <- human_rpkm[, c(1, 260:298)]





















