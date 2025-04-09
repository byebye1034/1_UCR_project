#
# https://mp.weixin.qq.com/s/LB5qqUShvm87xQ_hz7IM1g

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(WGCNA)
library(tidyverse)
library(stringr)
library(forcats)
library(openxlsx)
library(edgeR)
library(limma)

lwd_pt <- .pt*72.27/96

# data preparation --------------------------------------------------------

# exprDataMatrix
load("D:/R_project/UCR_project/01-data/45-GBM/GSE147352_ReadsCount.rdata")
readsCount <- ReadsCount
rm(ReadsCount)

# clinicalInfo
sampleInfo <- read.xlsx(xlsxFile = "01-data/45-GBM/GSE147352_clininfo.xlsx", sheet = 1)
table(sampleInfo$Type)
# GBM 85 LGG 18 Normal 15
sampleInfoGBM <- filter(sampleInfo, Type == "Normal" | Type == "GBM")
sampleInfoLGG <- filter(sampleInfo, Type == "Normal" | Type == "LGG")

save(sampleInfoGBM, sampleInfoLGG, file = "02-analysis/46-NCUCR_Genes_WGCNA/sampleInfoGBM_sampleInfoLGG.rdata")

readsCountGBM <- select(readsCount, sampleInfoGBM$sample.ID)
readsCountLGG <- select(readsCount, sampleInfoLGG$sample.ID)

# normalization and 2cpm
normFactorsGBM <- calcNormFactors(readsCountGBM, method = "TMM")
normFactorsLGG <- calcNormFactors(readsCountLGG, method = "TMM")

lcpmGBM <- edgeR::cpm(readsCountGBM, log = TRUE) # directly shift to log2(cpm + 1)
lcpmLGG <- edgeR::cpm(readsCountLGG, log = TRUE)

# row:sample column:gene
lcpmGBM <- data.frame(t(lcpmGBM))
lcpmLGG <- data.frame(t(lcpmLGG))

# get var of all genes
geneVarLGG <- apply(lcpmLGG, 2, var)
geneVarGBM <- apply(lcpmGBM, 2, var)

# var 2 data.frame
geneVarLGG <- data.frame(gene = names(geneVarLGG), variance = geneVarLGG)
geneVarGBM <- data.frame(gene = names(geneVarGBM), variance = geneVarGBM)

# get 75%
varianceThreshold <- quantile(geneVarLGG$variance, 0.75)
varianceThreshold <- quantile(geneVarGBM$variance, 0.75)

# 筛选出方差排在前 75% 的基因
top_genes <- geneVarGBM[geneVarGBM$variance >= varianceThreshold, ]
top_genes <- geneVarLGG[geneVarLGG$variance >= varianceThreshold, ]

# 输出筛选结果
lcpmLGG <- lcpmLGG[, names(lcpmLGG) %in% top_genes$gene]
lcpmGBM <- lcpmGBM[, names(lcpmGBM) %in% top_genes$gene]

save(lcpmGBM, lcpmLGG, file = "02-analysis/46-NCUCR_Genes_WGCNA/lcpmGBM_lcpmLGG.rdata")

# GBM data preparation for net building on server -------------------------

# judge data quality
gsg <- goodSamplesGenes(lcpmGBM, verbose = 3)
gsg$allOK

# 哎呀！质量不好，咱们执行以下代码块
if (!gsg$allOK)
{ 
  # 如果goodGenes属性中有不好的基因，打印并移除它们
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(lcpmGBM)[!gsg$goodGenes], collapse = ", "))) 
  
  # 如果goodSamples属性中有不好的样本，打印并移除它们
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(lcpmGBM)[!gsg$goodSamples], collapse = ", "))) 
  
  # 根据质量检查结果对象中的goodSamples和goodGenes属性，过滤数据集exp_tpm
  lcpmGBM <- lcpmGBM[gsg$goodSamples, gsg$goodGenes] 
}

# rejudge data quality
gsg <- goodSamplesGenes(lcpmGBM, verbose = 3)
gsg$allOK

# sample cluster tree
sampleTree <- hclust(dist(lcpmGBM), method = "average")

# Assign different colours to the phenotypic information for each sample and generate colour vectors
traitColors <- numbers2colors(as.numeric(factor(sampleInfoGBM$Type)), 
                              colors = rainbow(length(table(sampleInfoGBM$Type))), 
                              signed = FALSE)

pdf(file = "03-results/46-NCUCR_Genes_WGCNA/01-sampleTreeGBM.pdf", width = 7200/254, height = 1800/254)

plotDendroAndColors(sampleTree, 
                    traitColors,
                    groupLabels = names(sampleInfoGBM),
                    main = "Sample dendrogram and trait heatmap")

dev.off()

# remove SYC118 sample
lcpmGBM <- filter(lcpmGBM, !(rownames(lcpmGBM) %in% c("SYC118")))
sampleInfoGBM <- filter(sampleInfoGBM, !(sampleInfoGBM$sample.ID %in% c("SYC118")))
# recluster and plot

# sample cluster tree again
sampleTreeFinal <- hclust(dist(lcpmGBM), method = "average")

# Assign different colours to the phenotypic information for each sample and generate colour vectors
traitColors <- numbers2colors(as.numeric(factor(sampleInfoGBM$Type)), 
                              colors = rainbow(length(table(sampleInfoGBM$Type))), 
                              signed = FALSE)

pdf(file = "03-results/46-NCUCR_Genes_WGCNA/02-sampleTreeFianlGBM.pdf", width = 7200/254, height = 1800/254)

plotDendroAndColors(sampleTreeFinal, 
                    traitColors,
                    groupLabels = names(sampleInfoGBM),
                    main = "Sample dendrogram and trait heatmap")

dev.off()

# 将层次聚类树中的样本分成不同的簇
clustGBM <- cutreeStatic(sampleTreeFinal, cutHeight = 320, minSize = 5)
table(clustGBM)

keepSamples <- clustGBM == 1 | clustGBM == 2
lcpmGBM <- lcpmGBM[keepSamples, ]

nGenes <- ncol(lcpmGBM)
nGenes

nSamples <- nrow(lcpmGBM)
nSamples

dim(lcpmGBM)

## 挑选最佳软阈值
# 设置 power 参数选择范围，可以自行修改设定
powersGBM <- c(c(1:10), seq(from = 12, to = 20, by = 2))

# 选择最佳软阈值，获取各个阈值下的 R^2 和平均连接度
sftGBM <- pickSoftThreshold(lcpmGBM, powerVector = powersGBM, verbose = 5)
save(sftGBM, file = "02-analysis/46-NCUCR_Genes_WGCNA/sftGBM.rdata")

# 我们看一下软阈值选择结果
sftGBM
sftGBM$powerEstimate
sftGBM$fitIndices

# 绘制软阈值和拟合指标的关系图
pdf(file = "03-results/46-NCUCR_Genes_WGCNA/scaleIndependence_meanConnectivity_GBM.pdf", width = 900/254, height = 450/254)

par(mfrow = c(1,2), cex = 0.4)  # 将绘图区域分为1行2列

plot(sftGBM$fitIndices[, 1], 
     -sign(sftGBM$fitIndices[, 3]) * sftGBM$fitIndices[, 2],
     xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = paste("Scale independence"))

text(sftGBM$fitIndices[, 1], 
     -sign(sftGBM$fitIndices[, 3]) * sftGBM$fitIndices[, 2],
     labels = powersGBM, col = "red")

# 添加 R^2 水平线，使用 R^2 阈值为 0.90，官网建议最好是0.85或以上
abline(h = 0.90, col = "red")

# 绘制软阈值对平均连接度的影响图
plot(sftGBM$fitIndices[, 1], 
     sftGBM$fitIndices[, 5],
     xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))

text(sftGBM$fitIndices[, 1], 
     sftGBM$fitIndices[, 5], 
     labels = powersGBM, col = "red")

dev.off()

# build co-expression network
# net <- blockwiseModules(lcpmGBM, power = sftGBM$powerEstimate,
#                         maxBlockSize = 20000, TOMType = "unsigned",
#                         minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
#                         numericLabels = TRUE, pamRespectsDendro = FALSE,
#                         saveTOMs = F, verbose = 3)
# 
# save(net, file = "netGBM.rdata")

# GBM network analysis ----------------------------------------------------

load("D:/R_project/UCR_project/03-results/46-NCUCR_Genes_WGCNA/GSE146352/minModuleSize100/netGBM100.rdata")
load("D:/R_project/UCR_project/03-results/46-NCUCR_Genes_WGCNA/GSE146352/minModuleSize100/sftGBM.rdata")

# look module status
table(net$colors)

# 0 Indicates genes that are not assigned into any module. 

## module visualization
moduleColors <- labels2colors(net$colors)
table(moduleColors)

# hierarchical clustering tree
pdf(file = "03-results/46-NCUCR_Genes_WGCNA/GSE146352/minModuleSize100/04-clustersGBM.pdf", 
    width = 1350/254, height = 900/254)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, 
                    cex.colorLabels = 0.4, cex.dendroLabels = 0.4, cex.rowText = 0.4)
dev.off()

## corelation and significance

sampleInfoGBM$Type <- factor(sampleInfoGBM$Type)  # 将sampleInfoGBM$Type转换为因子变量
levels(sampleInfoGBM$Type)  # 显示因子变量的水平

# 绘制热图
pdf(file = "03-results/46-NCUCR_Genes_WGCNA/GSE146352/minModuleSize100/ModuleTraitRelationshipsGSE147352GBM.pdf", 
    width = 900/254, height = 900/254)

if(T){ 
  nGenes = ncol(lcpmGBM)  # 基因数目
  nSamples = nrow(lcpmGBM)  # 样本数目
  design <- model.matrix(~0+sampleInfoGBM$Type)  # 构建模型矩阵
  colnames(design)= levels(sampleInfoGBM$Type)  # 修改列名，获取分组信息
  MES0 <- moduleEigengenes(lcpmGBM, moduleColors)$eigengenes  # 计算模块特征向量
  MEs = orderMEs(MES0)  # 对模块特征向量排序
  moduleTraitCor <- cor(MEs, design, use = "p")  # 计算模块特征向量与表型的相关系数矩阵
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)  # 计算相关系数矩阵的p值
  textMatrix = paste(signif(moduleTraitCor,2), "(", 
                     signif(moduleTraitPvalue,1), ")", sep = "")  # 构建绘图时用的文本矩阵
  dim(textMatrix) = dim(moduleTraitCor)  # 修改文本矩阵的维度，与相关系数矩阵相同
  par(mar=c(3, 9, 3, 3), cex = 0.4)  # 设置绘图边距
  labeledHeatmap(Matrix = moduleTraitCor,  # 绘制带标签的热图
                 xLabels = colnames(design),  # x轴标签
                 xLabelsAngle = 0, xLabelsAdj = 0, 
                 yLabels = names(MEs),  # y轴标签
                 ySymbols = names(MEs),  # y轴符号
                 colorLabels = FALSE,  # 不显示颜色标签
                 colors = blueWhiteRed(50),  # 颜色范围
                 textMatrix = textMatrix,  # 显示文本矩阵
                 setStdMargins = FALSE,  # 不设置标准边距
                 cex.text = 0.6,  # 文本大小
                 zlim = c(-1,1),  # 颜色映射范围
                 main = paste("Module-trait relationships"))  # 绘图标题
}

dev.off()

# check Green and turquoise module genes
# ncUCR genes to check
ncUCRGenes <- read.table(file = "02-analysis/45-GBM/NCUCR_genes_id.txt", fill = NA)

# green
module="green"  # 选择要导出的模块
probes = colnames(lcpmGBM)  # 获取基因名称
inModule = (moduleColors == module)  # 找到属于当前模块的基因
modProbes=probes[inModule]  # 提取属于当前模块的基因名称
head(modProbes)  # 显示基因名称前几行
genesWithinGreenMod <- modProbes
head(genesWithinGreenMod)

print(intersect(ncUCRGenes$V1, genesWithinGreenMod))
# [1] "ENSG00000204929" "ENSG00000245526"

# turquoise
module="turquoise"  # 选择要导出的模块
probes = colnames(lcpmGBM)  # 获取基因名称
inModule = (moduleColors == module)  # 找到属于当前模块的基因
modProbes=probes[inModule]  # 提取属于当前模块的基因名称
head(modProbes)  # 显示基因名称前几行
genesWithinTurquoiseMod <- modProbes
head(genesWithinTurquoiseMod)

intersect(ncUCRGenes$V1, genesWithinTurquoiseMod)
# [1] "ENSG00000224592" "ENSG00000270087" "ENSG00000257545" "ENSG00000234377" "ENSG00000224243"
# [6] "ENSG00000257986" "ENSG00000233723" "ENSG00000228956" "ENSG00000145075" "ENSG00000231764"
# [11] "ENSG00000254102"

# LGG data preparation for net building on server -------------------------

# judge data quality
gsg <- goodSamplesGenes(lcpmLGG, verbose = 3)
gsg$allOK

# 哎呀！质量不好，咱们执行以下代码块
if (!gsg$allOK)
{ 
  # 如果goodGenes属性中有不好的基因，打印并移除它们
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(lcpmLGG)[!gsg$goodGenes], collapse = ", "))) 
  
  # 如果goodSamples属性中有不好的样本，打印并移除它们
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(lcpmLGG)[!gsg$goodSamples], collapse = ", "))) 
  
  # 根据质量检查结果对象中的goodSamples和goodGenes属性，过滤数据集exp_tpm
  lcpmLGG <- lcpmLGG[gsg$goodSamples, gsg$goodGenes] 
}

# sample cluster tree
sampleTree <- hclust(dist(lcpmLGG), method = "average")

# Assign different colours to the phenotypic information for each sample and generate colour vectors
traitColors <- numbers2colors(as.numeric(factor(sampleInfoLGG$Type)), 
                              colors = rainbow(length(table(sampleInfoLGG$Type))), 
                              signed = FALSE)

pdf(file = "03-results/46-NCUCR_Genes_WGCNA/01-sampleTreeLGG.pdf", width = 7200/254, height = 1800/254)

plotDendroAndColors(sampleTree, 
                    traitColors,
                    groupLabels = names(sampleInfoLGG),
                    main = "Sample dendrogram and trait heatmap")

dev.off()

# remove SYC189 SYC209 sample
lcpmLGG <- filter(lcpmLGG, !(rownames(lcpmLGG) %in% c("SYC189", "SYC209")))
sampleInfoLGG <- filter(sampleInfoLGG, !(sampleInfoLGG$sample.ID %in% c("SYC189", "SYC209")))
# recluster and plot

# 将层次聚类树中的样本分成不同的簇
clustLGG <- cutreeStatic(sampleTree, cutHeight = 300, minSize = 5)
table(clustLGG)

keepSamples <- clustLGG == 1
lcpmLGG <- lcpmLGG[keepSamples, ]

nGenes <- ncol(lcpmLGG)
nGenes

nSamples <- nrow(lcpmLGG)
nSamples

dim(lcpmLGG)

## 挑选最佳软阈值
# 设置 power 参数选择范围，可以自行修改设定
powersLGG <- c(c(1:10), seq(from = 12, to = 20, by = 2))

# 选择最佳软阈值，获取各个阈值下的 R^2 和平均连接度
sftLGG <- pickSoftThreshold(lcpmLGG, powerVector = powersLGG, verbose = 5)
save(sftLGG, file = "02-analysis/46-NCUCR_Genes_WGCNA/sftLGG.rdata")

# 我们看一下软阈值选择结果
sftLGG
sftLGG$powerEstimate
sftLGG$fitIndices

# 绘制软阈值和拟合指标的关系图
pdf(file = "03-results/46-NCUCR_Genes_WGCNA/scaleIndependence_meanConnectivity_LGG.pdf", width = 900/254, height = 450/254)

par(mfrow = c(1,2), cex = 0.4)  # 将绘图区域分为1行2列

plot(sftLGG$fitIndices[, 1], 
     -sign(sftLGG$fitIndices[, 3]) * sftLGG$fitIndices[, 2],
     xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = paste("Scale independence"))

text(sftLGG$fitIndices[, 1], 
     -sign(sftLGG$fitIndices[, 3]) * sftLGG$fitIndices[, 2],
     labels = powersLGG, col = "red")

# 添加 R^2 水平线，使用 R^2 阈值为 0.90，官网建议最好是0.85或以上
abline(h = 0.90, col = "red")

# 绘制软阈值对平均连接度的影响图
plot(sftLGG$fitIndices[, 1], 
     sftLGG$fitIndices[, 5],
     xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))

text(sftLGG$fitIndices[, 1], 
     sftLGG$fitIndices[, 5], 
     labels = powersLGG, col = "red")

dev.off()

# LGG network analysis ----------------------------------------------------

load("D:/R_project/UCR_project/03-results/46-NCUCR_Genes_WGCNA/GSE146352/minModuleSize100/netLGG100.rdata")
load("D:/R_project/UCR_project/03-results/46-NCUCR_Genes_WGCNA/GSE146352/minModuleSize100/sftLGG.rdata")

# look module status
table(net$colors)

# 0 Indicates genes that are not assigned into any module. 

## module visualization
moduleColors <- labels2colors(net$colors)
table(moduleColors)

# hierarchical clustering tree
pdf(file = "03-results/46-NCUCR_Genes_WGCNA/GSE146352/minModuleSize100/04-clustersLGG.pdf", 
    width = 1350/254, height = 900/254)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, 
                    cex.colorLabels = 0.4, cex.dendroLabels = 0.4, cex.rowText = 0.4)
dev.off()

## corelation and significance

sampleInfoLGG$Type <- factor(sampleInfoLGG$Type)  # 将sampleInfoGBM$Type转换为因子变量
levels(sampleInfoLGG$Type)  # 显示因子变量的水平

lcpmLGG <- lcpmLGG[c(1:31), ]

# 绘制热图
pdf(file = "03-results/46-NCUCR_Genes_WGCNA/GSE146352/minModuleSize100/ModuleTraitRelationshipsGSE147352LGG.pdf", 
    width = 900/254, height = 900/254)

if(T){ 
  nGenes = ncol(lcpmLGG)  # 基因数目
  nSamples = nrow(lcpmLGG)  # 样本数目
  design <- model.matrix(~0+sampleInfoLGG$Type)  # 构建模型矩阵
  colnames(design)= levels(sampleInfoLGG$Type)  # 修改列名，获取分组信息
  MES0 <- moduleEigengenes(lcpmLGG, moduleColors)$eigengenes  # 计算模块特征向量
  MEs = orderMEs(MES0)  # 对模块特征向量排序
  moduleTraitCor <- cor(MEs, design, use = "p")  # 计算模块特征向量与表型的相关系数矩阵
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)  # 计算相关系数矩阵的p值
  textMatrix = paste(signif(moduleTraitCor,2), "(", 
                     signif(moduleTraitPvalue,1), ")", sep = "")  # 构建绘图时用的文本矩阵
  dim(textMatrix) = dim(moduleTraitCor)  # 修改文本矩阵的维度，与相关系数矩阵相同
  par(mar=c(3, 9, 3, 3), cex = 0.4)  # 设置绘图边距
  labeledHeatmap(Matrix = moduleTraitCor,  # 绘制带标签的热图
                 xLabels = colnames(design),  # x轴标签
                 xLabelsAngle = 0, xLabelsAdj = 0, 
                 yLabels = names(MEs),  # y轴标签
                 ySymbols = names(MEs),  # y轴符号
                 colorLabels = FALSE,  # 不显示颜色标签
                 colors = blueWhiteRed(50),  # 颜色范围
                 textMatrix = textMatrix,  # 显示文本矩阵
                 setStdMargins = FALSE,  # 不设置标准边距
                 cex.text = 0.6,  # 文本大小
                 zlim = c(-1,1),  # 颜色映射范围
                 main = paste("Module-trait relationships"))  # 绘图标题
}

dev.off()

# check Green and turquoise module genes
# ncUCR genes to check
ncUCRGenes <- read.table(file = "02-analysis/45-GBM/NCUCR_genes_id.txt", fill = NA)

# green
module="yellow"  # 选择要导出的模块
probes = colnames(lcpmLGG)  # 获取基因名称
inModule = (moduleColors == module)  # 找到属于当前模块的基因
modProbes=probes[inModule]  # 提取属于当前模块的基因名称
head(modProbes)  # 显示基因名称前几行
genesWithinBlueMod <- modProbes
head(genesWithinBlueMod)

print(intersect(ncUCRGenes$V1, genesWithinBlueMod))
# [1] "ENSG00000204929" "ENSG00000245526"

# turquoise
module="turquoise"  # 选择要导出的模块
probes = colnames(lcpmLGG)  # 获取基因名称
inModule = (moduleColors == module)  # 找到属于当前模块的基因
modProbes=probes[inModule]  # 提取属于当前模块的基因名称
head(modProbes)  # 显示基因名称前几行
genesWithinTurquoiseMod <- modProbes
head(genesWithinTurquoiseMod)

intersect(ncUCRGenes$V1, genesWithinTurquoiseMod)


























