# WGCNA identify key gene module
# https://mp.weixin.qq.com/s/fCvLizKQNWDQKeWBuSG3UQ
# https://mp.weixin.qq.com/s/5OUY5KDwgi05MlFrV_7Qjw

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(WGCNA)
library(tidyverse)
library(stringr)
library(forcats)

lwd_pt <- .pt*72.27/96

# data preparation --------------------------------------------------------

## expr matrix tpm
load("D:/R_project/UCR_project/01-data/45-GBM/output_mRNA_lncRNA_expr/TCGA-GBM_lncrna_expr_tpm.rdata")
load("D:/R_project/UCR_project/01-data/45-GBM/output_mRNA_lncRNA_expr/TCGA-GBM_mrna_expr_tpm.rdata")

identical(colnames(lncrna_expr_tpm), colnames(mrna_expr_tpm))

GBM_Expr_Mat <- rbind(mrna_expr_tpm, lncrna_expr_tpm)
GBM_Expr_Mat <- log2(GBM_Expr_Mat + 1)
GBM_Expr_Mat <- data.frame(t(GBM_Expr_Mat))
# rownames_to_column
# GBM_Expr_Mat <- GBM_Expr_Mat %>% 
#   rownames_to_column(var = "ID")

save(GBM_Expr_Mat, file = "02-analysis/46-NCUCR_Genes_WGCNA/GBM_Expr_Mat.rdata")

## clinical info preparation
load("D:/R_project/UCR_project/01-data/45-GBM/output_mRNA_lncRNA_expr/TCGA-GBM_clinicalSE.rdata")

# get barcode and tissue type
GBM_Clinical <- clinicalSE %>%
  select(barcode, tissue_type) %>%  # select column need
  mutate(tissue_type_no = ifelse(tissue_type == "Tumor", 1, 2))  # replace tumor:1 normal:0
GBM_Clinical <- GBM_Clinical[, -1]

save(GBM_Clinical, file = "02-analysis/46-NCUCR_Genes_WGCNA/GBM_Clinical.rdata")

# start WGCNA analysis ----------------------------------------------------

load("02-analysis/46-NCUCR_Genes_WGCNA/GBM_Expr_Mat.rdata")
load("02-analysis/46-NCUCR_Genes_WGCNA/GBM_Clinical.rdata")

# judge data quality
gsg <- goodSamplesGenes(GBM_Expr_Mat, verbose = 3)
gsg$allOK

# 哎呀！质量不好，咱们执行以下代码块
if (!gsg$allOK)
{ 
  # 如果goodGenes属性中有不好的基因，打印并移除它们
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(GBM_Expr_Mat)[!gsg$goodGenes], collapse = ", "))) 
  
  # 如果goodSamples属性中有不好的样本，打印并移除它们
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(GBM_Expr_Mat)[!gsg$goodSamples], collapse = ", "))) 
  
  # 根据质量检查结果对象中的goodSamples和goodGenes属性，过滤数据集exp_tpm
  GBM_Expr_Mat <- GBM_Expr_Mat[gsg$goodSamples, gsg$goodGenes] 
}

# 第一步，咱们就是构建样本的系统聚类树
# 主要为了查看是否有离群样本
sampleTree <- hclust(dist(GBM_Expr_Mat), method = "average")

# Assign different colours to the phenotypic information for each sample and generate colour vectors
traitColors <- numbers2colors(as.numeric(factor(GBM_Clinical$tissue_type)), 
                              colors = rainbow(length(table(GBM_Clinical$tissue_type))), 
                              signed = FALSE)

pdf(file = "02-analysis/46-NCUCR_Genes_WGCNA/SampleTree.pdf", width = 7200/254, height = 1800/254)

plotDendroAndColors(sampleTree, 
                    traitColors,
                    groupLabels = names(GBM_Clinical),
                    main = "Sample dendrogram and trait heatmap")

dev.off()

# 可以在图中画线，剔除离群样本，注意：不需要剔除样本画线就别运行这行！
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 0.4)
abline(h = 150, col = "red")

# 将层次聚类树中的样本分成不同的簇
clust <- cutreeStatic(sampleTree, cutHeight = 300, minSize = 5)

# cutHeight = 200：用于指定在层次聚类树中切割的高度。在树状结构中，高度表示样本之间的相似性或距离。
# 通过指定 cutHeight，你可以控制在哪个高度水平切割树，从而确定最终的簇数。

# minSize = 10：用于指定最小簇的大小。在进行切割时，如果某个簇的大小小于 minSize，
# 则可能会合并到其他簇中，以确保生成的簇都具有足够的样本数。

# 查看详情
table(clust)

keepSamples <- clust == 1 | clust == 5
GBM_datExpr <- GBM_Expr_Mat[keepSamples, ]

nGenes <- ncol(GBM_datExpr)
nGenes

nSamples <- nrow(GBM_datExpr)
nSamples

dim(GBM_datExpr)

GBM_Clinical <- filter(GBM_Clinical, rownames(GBM_Clinical) %in% rownames(GBM_datExpr))
traitColors <- numbers2colors(as.numeric(factor(GBM_Clinical$tissue_type)), 
                              colors = rainbow(length(table(GBM_Clinical$tissue_type))), 
                              signed = FALSE)

# 剔除离群样本后，咱们可以再重新聚类一下
sampleTree_final <- hclust(dist(GBM_datExpr), method = "average")
pdf(file = "02-analysis/46-NCUCR_Genes_WGCNA/SampleTreeFinal.pdf", width = 7200/254, height = 1800/254)
plotDendroAndColors(sampleTree_final, 
                    traitColors,
                    groupLabels = names(GBM_Clinical),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

## 挑选最佳软阈值

# 设置 power 参数选择范围，可以自行修改设定
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

# 选择最佳软阈值，获取各个阈值下的 R^2 和平均连接度
sft <- pickSoftThreshold(GBM_datExpr, powerVector = powers, verbose = 5)
save(sft, file = "02-analysis/46-NCUCR_Genes_WGCNA/sft.rdata")

# 我们看一下软阈值选择结果
sft
sft$powerEstimate
sft$fitIndices

# 绘制软阈值和拟合指标的关系图
pdf(file = "03-results/46-NCUCR_Genes_WGCNA/ScaleIndependence_MeanConnectivity.pdf", width = 900/254, height = 450/254)

par(mfrow = c(1,2), cex = 0.4)  # 将绘图区域分为1行2列

plot(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = paste("Scale independence"))

text(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")

# 添加 R^2 水平线，使用 R^2 阈值为 0.90，官网建议最好是0.85或以上
abline(h = 0.90, col = "red")

# 绘制软阈值对平均连接度的影响图
plot(sft$fitIndices[, 1], 
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))

text(sft$fitIndices[, 1], 
     sft$fitIndices[, 5], 
     labels = powers, col = "red")

dev.off()

# 正式构建加权共表达网络

# run in server
# net <- blockwiseModules(GBM_datExpr, power = sft$powerEstimate,
#                         maxBlockSize = 20000, TOMType = "unsigned",
#                         minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
#                         numericLabels = TRUE, pamRespectsDendro = FALSE,
#                         saveTOMs = F, verbose = 3)

# load net.rdata from nebula server
load("02-analysis/46-NCUCR_Genes_WGCNA/net.rdata")
load("02-analysis/46-NCUCR_Genes_WGCNA/net_20000_100.rdata")

# look module status
table(net$colors)

# 0 Indicates genes that are not assigned into any module. 

## module visualization
moduleColors <- labels2colors(net$colors)
table(moduleColors)

# hierarchical clustering tree
pdf(file = "03-results/46-NCUCR_Genes_WGCNA/clusters.pdf", width = 1350/254, height = 900/254)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, 
                    cex.colorLabels = 0.4, cex.dendroLabels = 0.4, cex.rowText = 0.4)
dev.off()

# correlation and significnce ---------------------------------------------

# Correlation and significance between modules and traits of interest, heatmap of 
# correlation between modules and traits, boxplot of correlation between modules 
# and traits, and scatterplot of correlation between genes, modules and phenotypes

# # 颜色标签
# moduleLables <- net$colors
# moduleColors <- labels2colors(net$colors)
# 
# # ME值，也就是获取eigengenes，每个ME代表一个模块
# MEs <- net$MEs
# head(MEs)[1:5, 1:5]
# #                                      ME59         ME64          ME50          ME98         ME48
# # TCGA-06-0211-01B-01R-1849-01 -0.017687568 -0.015219907 -0.0171138780 -0.0098602255 -0.016872903
# # TCGA-06-0211-02A-02R-2005-01 -0.009734677 -0.015110112 -0.0222683328 -0.0025806605 -0.018511224
# # TCGA-16-1045-01B-01R-1850-01 -0.019190550 -0.012940362  0.0010119309 -0.0226333359 -0.024101413
# # TCGA-06-0125-01A-01R-1849-01 -0.010079356  0.005813045  0.0063616145  0.0002157516 -0.003852471
# # TCGA-06-0125-02A-11R-2005-01 -0.013971219 -0.025807521 -0.0007979277  0.0010568469 -0.018196396
# 
# geneTree <- net$dendrograms[[1]]
# 
# save(moduleLables, moduleColors, MEs, geneTree, 
#      file = "02-analysis/46-NCUCR_Genes_WGCNA//networkConstruction.rdata")
# 
# # 将基因模块与性状进行关联
# 
# # 获取eigengenes，用颜色标签计算ME值
# MEList <-  moduleEigengenes(GBM_datExpr, colors = moduleColors)
# MEs0 <- MEList$eigengenes
# 
# # 查看用颜色标签计算的ME值
# head(MEs0)[1:5, 1:5]
# #                              MEantiquewhite2 MEantiquewhite4    MEbisque4     MEblack       MEblue
# # TCGA-06-0211-01B-01R-1849-01     -0.01873829    -0.006100871  0.154665403 -0.07526370  0.001561828
# # TCGA-06-0211-02A-02R-2005-01     -0.03093754    -0.016835632  0.034280961 -0.08901012 -0.057258089
# # TCGA-16-1045-01B-01R-1850-01      0.02264226    -0.028709629 -0.055032977  0.02062123 -0.001007145
# # TCGA-06-0125-01A-01R-1849-01     -0.01814417    -0.011323409  0.168629645 -0.01683616  0.162279079
# # TCGA-06-0125-02A-11R-2005-01     -0.02234831    -0.007594725 -0.004868453 -0.02511028 -0.093858617
# # 可以看到我们的列名已经变成了颜色，不同的颜色代表不同的模块
# 
# # 排序
# MEs <- orderMEs(MEs0)
# head(MEs)[1:5, 1:5]
# 
# # 计算每个模块和每个性状之间的相关性
# GBM_Clinical <- GBM_Clinical %>% 
#   filter(barcode %in% rownames(MEs))
# moduleTraitCor <- cor(MEs, GBM_Clinical, use = "p")
# head(moduleTraitCor)

GBM_Clinical$tissue_type <- factor(GBM_Clinical$tissue_type)  # 将GBM_Clinical$tissue_type转换为因子变量
levels(GBM_Clinical$tissue_type)  # 显示因子变量的水平

# 绘制热图
pdf(file = "03-results/46-NCUCR_Genes_WGCNA/Module-trait_relationships.pdf", width = 900/254, height = 900/254)

if(T){ 
  nGenes = ncol(GBM_datExpr)  # 基因数目
  nSamples = nrow(GBM_datExpr)  # 样本数目
  design <- model.matrix(~0+GBM_Clinical$tissue_type)  # 构建模型矩阵
  colnames(design)= levels(GBM_Clinical$tissue_type)  # 修改列名，获取分组信息
  MES0 <- moduleEigengenes(GBM_datExpr, moduleColors)$eigengenes  # 计算模块特征向量
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

# MM:Module Membership GS:Gene Significance -------------------------------

# # 模块内基因与表型数据关联
# 
# # 我们发现与 auc 最相关的是 saddlebrown 模块
# # names (colors) of the modules
# 
# # datExpr 表示每个基因在每个样本中的表达量
# # MEs 表示每个模块在每个样本中的模块特征值
# # moduleColors 表示每个基因所属的模块颜色
# 
# # 获取模块名称
# modNames <- substring(names(MEs), 3)
# modNames
# 
# # 计算模块与基因的相关性矩阵
# geneModuleMembership <- as.data.frame(cor(GBM_datExpr, MEs, use = "p"))
# MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
# names(geneModuleMembership) <- paste("MM", modNames, sep = "")
# names(MMPvalue) <- paste("p.MM", modNames, sep = "");
# geneModuleMembership[1:5, 1:5]
# #        MMdarkseagreen4  MMskyblue2 MMmediumpurple1 MMthistle2 MMdarkslateblue
# # MT.CO1      0.04551230 -0.11045144     -0.08858152  0.2704998      0.04299649
# # MT.CO2      0.07884336 -0.09290044     -0.06535712  0.2540856      0.09345345
# # MT.CO3      0.11959701 -0.19091880     -0.11552775  0.1977291      0.05767739
# # MT.ND4      0.04569057 -0.08272758     -0.09549052  0.1894112      0.04365731
# # FTL        -0.04917593 -0.14420409      0.10665293 -0.1859320     -0.02109511
# 
# # 计算性状与基因的相关性矩阵 
# # 只有连续型性状才能进行计算，如果是离散变量，在构建样本表时就转为0-1矩阵。
# geneTraitSignificance <- as.data.frame(cor(GBM_datExpr, GBM_Clinical, use = "p"))
# GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
# names(geneTraitSignificance) <- paste("GS.", names(GBM_Clinical), sep = "")
# names(GSPvalue) <- paste("p.GS.", names(GBM_Clinical), sep = "")
# head(geneTraitSignificance)
# #         GS.barcode GS.tissue_type
# # MT.CO1          NA     -0.2112273
# # MT.CO2          NA     -0.2078630
# # MT.CO3          NA     -0.2720754
# # MT.ND4          NA     -0.1902578
# # FTL             NA      0.3635055
# # MT.ATP6         NA     -0.1244849
# 
# # 最后把两个相关性矩阵联合起来，指定感兴趣模块进行分析
# module = "turquoise"
# pheno = "tissue_type"
# modNames = substring(names(MEs), 3)
# 
# # 获取关注的列
# module_column = match(module, modNames)
# pheno_column = match(pheno, colnames(GBM_Clinical))
# 
# # 获取模块内的基因
# moduleGenes <- moduleColors == module
# 
# # 可视化基因与模块、表型的相关性，绘制散点图（想看所有的咱可以批量作图，咦，后续是不是可以分享一下！）
# pdf(file = "03-results/46-NCUCR_Genes_WGCNA/MM_GS_Cor.pdf", width = 900/254, height = 900/254)
# par(mar = c(1, 3, 3, 3))
# verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
#                    abs(geneTraitSignificance[moduleGenes, 1]),
#                    xlab = paste("Module Membership in", module, "module"),
#                    ylab = paste("Gene significance for LRG"),
#                    main = paste("Module membership vs. gene significance\n"),
#                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
# dev.off()

pdf(file = "03-results/46-NCUCR_Genes_WGCNA/MM_GS_Cor.pdf", width = 480/254, height = 480/254)
if(T){
  modNames = substring(names(MEs), 3)  # 提取模块名称
  geneModuleMembership = as.data.frame(cor(GBM_datExpr, MEs,
                                           use = "p",method = "spearman"))  #计算基因与模块的相关系数矩阵
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))  # 计算相关系数矩阵的p值
  names(geneModuleMembership) = paste("MM", modNames, sep="")  # 修改列名
  names(MMPvalue) = paste("p.MM", modNames, sep="")  # 修改列名
  
  geneTraitSignificance <- as.data.frame(cor(GBM_datExpr, GBM_Clinical$tissue_type_no, use = "p"))  # 计算基因与表型的相关系数矩阵
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))  # 计算相关系数矩阵的p值
  names(geneTraitSignificance)<- paste("GS.",names(GBM_Clinical$tissue_type), sep = "")  # 修改列名
  names(GSPvalue)<-paste("GS.", names(GBM_Clinical$tissue_type), sep = "")  # 修改列名
  
  selectModule<-c("turquoise")  # 选择要绘制的模块（这里只选择了一个）
  #selectModule <- modNames  # 批量作图
  
  #par(mfrow=c(ceiling(length(selectModule)/2),2))  # 设置绘图区域，批量作图开始
  par(mfrow = c(1, 1), cex = 0.4)
  for(module in selectModule){
    column <- match(module,selectModule)  # 找到当前模块在geneModuleMembership中的列号
    print(module)  # 输出当前处理的模块名称
    moduleGenes <- moduleColors==module  # 找到属于当前模块的基因
    
    # 绘制散点图，横轴为基因在当前模块中的模块内连通性，纵轴为基因与表型的相关系数
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = paste("Gene significance for", module, "module"),
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 0.4, cex.lab = 0.4, cex.axis = 0.4, col = module)
  }
}
dev.off()

# view genes within a module ----------------------------------------------

table(net$colors)
moduleLabels3 = net$colors
moduleColors3 = labels2colors(net$colors)
MEs3 = net$MEs
geneTree3 = net$dendrograms[[1]]
gm3 = data.frame(net$colors)
gm3$color = moduleColors
head(gm3)
genes3 = split(rownames(gm3),gm3$color)

# 假设 genes3 是一个包含 110 个元素的列表
# 每个元素都是一个向量或列表

# 使用 lapply 创建一个数据框的列表
df_list <- lapply(names(genes3), function(element_name) {
  data.frame(gene = genes3[[element_name]], 
             element_name = element_name, 
             stringsAsFactors = FALSE)
})

# 将所有数据框合并成一个长型数据框
genes_widthin_module <- do.call(rbind, df_list)

save(genes_widthin_module, file = "03-results/46-NCUCR_Genes_WGCNA/genes_within_module.rdata")

# convert ensembl id in website
converted_NCUCR_genes <- read.table(file = "02-analysis/45-GBM/NCUCR_genes_id.txt", 
                                    sep = "\t", blank.lines.skip = TRUE)
# remove null value
converted_NCUCR_genes <- converted_NCUCR_genes[!(converted_NCUCR_genes$V2 == ""), ]

converted_NCUCR_genes_module <- subset(genes_widthin_module, 
                                       genes_widthin_module$gene %in% converted_NCUCR_genes$V2)



# merge modules -----------------------------------------------------------
# https://www.jianshu.com/p/6d618c748143

# MEList = moduleEigengenes(GBM_datExpr, colors = moduleColors)
# MEs = MEList$eigengenes
# # Calculate dissimilarity of module eigengenes
# MEDiss = 1-cor(MEs);
# # Cluster module eigengenes
# METree = hclust(as.dist(MEDiss), method = "average")
# # Plot the result
# #sizeGrWindow(7, 6)
# pdf(file="6_Clustering of module eigengenes.pdf",width=7,height=6)
# plot(METree, main = "Clustering of module eigengenes",
#      xlab = "", sub = "")
# MEDissThres = 0.94######剪切高度可修改
# # Plot the cut line into the dendrogram
# abline(h=MEDissThres, col = "red")
# dev.off()
