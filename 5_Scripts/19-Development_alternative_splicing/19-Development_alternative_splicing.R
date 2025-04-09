#

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)

human_psi <- read.table(file = "01-data/19-Development_alternative_splicing/human.psi", 
                        sep = ",", header = TRUE)
print(colnames(human_psi))

## brain
brain_human_psi <- human_psi[, c(grep(pattern = "Brain", colnames(human_psi)))]
# 提取列名中的数字
nums <- as.numeric(sub(".*\\.(\\d+)$", "\\1", names(brain_human_psi)))
# 根据提取的数字对列进行排序
brain_human_psi <- brain_human_psi[, order(nums)]
brain_human_psi <- cbind(human_psi$X, brain_human_psi)

## cerebellum
cerebellum_human_psi <- human_psi[, c(grep(pattern = "Cerebellum", colnames(human_psi)))]
nums <- as.numeric(sub(".*\\.(\\d+)$", "\\1", names(cerebellum_human_psi)))
cerebellum_human_psi <- cerebellum_human_psi[, order(nums)]
cerebellum_human_psi <- cbind(human_psi$X, cerebellum_human_psi)

## heart
heart_human_psi <- human_psi[, c(grep(pattern = "Heart", colnames(human_psi)))]
nums <- as.numeric(sub(".*\\.(\\d+)$", "\\1", names(heart_human_psi)))
heart_human_psi <- heart_human_psi[, order(nums)]
heart_human_psi <- cbind(human_psi$X, heart_human_psi)

## kidney
kidney_human_psi <- human_psi[, c(grep(pattern = "Kidney", colnames(human_psi)))]
nums <- as.numeric(sub(".*\\.(\\d+)$", "\\1", names(kidney_human_psi)))
kidney_human_psi <- kidney_human_psi[, order(nums)]
kidney_human_psi <- cbind(human_psi$X, kidney_human_psi)

## liver
liver_human_psi <- human_psi[, c(grep(pattern = "Liver", colnames(human_psi)))]
nums <- as.numeric(sub(".*\\.(\\d+)$", "\\1", names(liver_human_psi)))
liver_human_psi <- liver_human_psi[, order(nums)]
liver_human_psi <- cbind(human_psi$X, liver_human_psi)

## ovary
ovary_human_psi <- human_psi[, c(grep(pattern = "Ovary", colnames(human_psi)))]
nums <- as.numeric(sub(".*\\.(\\d+)$", "\\1", names(ovary_human_psi)))
ovary_human_psi <- ovary_human_psi[, order(nums)]
ovary_human_psi <- cbind(human_psi$X, ovary_human_psi)

## testis
testis_human_psi <- human_psi[, c(grep(pattern = "Testis", colnames(human_psi)))]
nums <- as.numeric(sub(".*\\.(\\d+)$", "\\1", names(testis_human_psi)))
testis_human_psi <- testis_human_psi[, order(nums)]
testis_human_psi <- cbind(human_psi$X, testis_human_psi)

save(brain_human_psi, 
     cerebellum_human_psi, 
     heart_human_psi, 
     kidney_human_psi, 
     liver_human_psi, 
     ovary_human_psi, 
     testis_human_psi, 
     file = "01-data/19-Development_alternative_splicing/human_organ_psi.Rdata")













