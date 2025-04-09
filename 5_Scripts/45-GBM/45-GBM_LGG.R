# identify the change of NUCR genes in GBM and LGG

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library(easyTCGA)
library(ggplot2)
library(ggprism)
library(sysfonts)
library(showtext)
library(ggview)
library(cowplot)
library(ggrepel)
library(survival)
library(cutoff)

# getmrnaexpr get specific(GBM and LGG) tumor data
# download data to current path dir
setwd(dir = "D:/R_project/UCR_project/01-data/45-GBM/")
getmrnaexpr(project = "TCGA-GBM")
getmrnaexpr(project = "TCGA-LGG")

setwd(dir = "D:/R_project/UCR_project/")

# TCGA GBM ----------------------------------------------------------------

# NCUCR genes in GBM
load("D:/R_project/UCR_project/01-data/45-GBM/output_mRNA_lncRNA_expr/TCGA-GBM_lncrna_expr_count.rdata")

# get NCUCR genes
NCUCR_genes <- read.table(file = "02-analysis/16-New_classification/02-UCR_from_Non_Coding_RNA(66).bed")
NCUCR_genes <- NCUCR_genes[, c(8, 9)]
NCUCR_genes <- unique(NCUCR_genes)
print(noquote(NCUCR_genes[["V8"]]))

# convert ensembl id in website
converted_NCUCR_genes <- read.table(file = "02-analysis/45-GBM/NCUCR_genes_id.txt", 
                                    sep = "\t", blank.lines.skip = TRUE)
# remove null value
converted_NCUCR_genes <- converted_NCUCR_genes[!(converted_NCUCR_genes$V2 == ""), ]

converted_NCUCR_genes_exp <- subset(lncrna_expr_count, 
                                    rownames(lncrna_expr_count) %in% converted_NCUCR_genes$V2)

diff_analysis_results <- diff_analysis(exprset = lncrna_expr_count, 
                                       project = "GBM", 
                                       is_count = TRUE)
diff_analysis_deseq2_results <- diff_analysis_results$deg_deseq2
diff_analysis_deseq2_results <- merge(diff_analysis_deseq2_results, 
                                      converted_NCUCR_genes, 
                                      by.x = "genesymbol", by.y = "V2")
diff_analysis_deseq2_results <- diff_analysis_deseq2_results %>% 
  mutate(type = ifelse(abs(log2FoldChange) < 1 | pvalue >= 0.05, "nosig", 
                       ifelse(log2FoldChange > 1, "up", "down")))

# build for_label
marker <- diff_analysis_deseq2_results[!(diff_analysis_deseq2_results$type %in% "nosig"), ]

lwd_pt <- .pt*72.27/96

font_add(family = "Arial", regular = "C:/Windows/Fonts/arial.ttf")
font_families()

showtext_begin()

NCUCR_genes_volcano_GBM <- ggplot(data = diff_analysis_deseq2_results) +
  
  geom_point(mapping = aes(x = log2FoldChange, y = -log10(pvalue), 
                           fill = type, color = type), size = 0.8) +
  
  scale_fill_manual(values = c("blue", "grey", "red")) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  
  scale_x_continuous(limits = c(-5, 5), expand = c(0, 0)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5/lwd_pt) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.5/lwd_pt) +
  
  geom_text_repel(
    data = subset(diff_analysis_deseq2_results, genesymbol %in% marker$genesymbol), 
    aes(x = log2FoldChange, y = -log10(pvalue), label = genesymbol), 
    alpha = 0.5, 
    max.overlaps = 15, 
    size = unit(1.5, "pt")
  ) +
  
  theme_prism(base_size = 7, 
              base_family = "Arial", 
              base_fontface = "plain", 
              base_line_size = 0.5/lwd_pt, 
              base_rect_size = 0.5/lwd_pt, 
              border = TRUE) +
  
  theme(
    legend.position = "none", 
    legend.key.size = unit(0.25, "cm"), 
    aspect.ratio = 1
  )

NCUCR_genes_volcano_GBM

ggview::ggview(NCUCR_genes_volcano_GBM, 
               width = 4.8, height = 6, units = "cm", dpi = 1200)

showtext_end()

# easyTCGA survival analysis
# https://mp.weixin.qq.com/s/Wj2dmk5hPhkbGIkdUylGOg
load("D:/R_project/UCR_project/01-data/45-GBM/output_mRNA_lncRNA_expr/TCGA-GBM_clinicalSE.rdata")

surv_res <- batch_survival(exprset = lncrna_expr_count, 
                           clin = clinicalSE, 
                           is_count = TRUE, 
                           optimal_cut = T)
save(surv_res, file = "03-results/45-GBM/GBM_surv_res.Rdata")

GBM_log_rank_res <- surv_res$res.logrank
GBM_log_rank_res <- filter(GBM_log_rank_res, gene %in% NCUCR_genes$V9)
print(GBM_log_rank_res)

## survival analysis
load("D:/R_project/UCR_project/01-data/45-GBM/output_mRNA_lncRNA_expr/TCGA-GBM_clinicalSE.rdata")
load("D:/R_project/UCR_project/01-data/45-GBM/output_mRNA_lncRNA_expr/TCGA-GBM_lncrna_expr_count.rdata")

# 准备临床数据
GBM_clinical <- clinicalSE[, c("days_to_last_follow_up", "vital_status")]
names(GBM_clinical) <- c("time", "event")
GBM_clinical$event <- ifelse(GBM_clinical$event=="Dead", 1, 0)

pdf("03-results/45-GBM/surv_LINC01019.pdf", width = 633/254, height = 760/254)

marker = "LINC01117"

res <- plot_KM(exprset = lncrna_expr_count, 
               marker = marker, 
               clin = GBM_clinical, 
               optimal_cut = T, return_data = TRUE)

dev.off()

# # get needed info
# demo <- c("sample_submitter_id","gender","race", "age_at_index", "tissue_type", 
#           "vital_status","days_to_death","days_to_last_follow_up")
# 
# matrix = clinicalSE[, demo] #筛选需要的临床信息
# matrix <- filter(matrix, tissue_type == "Tumor")
# table(matrix$tissue_type)
# 
# colnames(matrix) <- c("ID","Gender","Race","Age","Tissue_type", 
#                       "Status","days_to_death","days_to_last_follow_up")
# 
# matrix = matrix[matrix$Status %in% c('Alive','Dead'),] # 排除结局为"Not Reported"的Sample
# 
# matrix$days_to_last_follow_up[is.na(matrix$days_to_last_follow_up)] = 0
# matrix$days_to_death[is.na(matrix$days_to_death)] = 0
# matrix$days <- ifelse(matrix$Status=='Alive', matrix$days_to_last_follow_up, matrix$days_to_death)
# matrix$days <- as.numeric(unlist(matrix$days))
# matrix$month=round(matrix$days/30, 0) #以month为单位，小数不保留
# 
# # # https://mp.weixin.qq.com/s/0h6O8Gzm_odXECC0h_Ijkw
# # # 用年龄做分组
# # matrix$group = ifelse(matrix$Age > median(matrix$Age),'old','young')
# # meta <- select(matrix,c("month","group","Status"))
# # 
# # survdata = Surv(time = meta$month,            #生存时间数据
# #                 event = meta$Status =='Dead') #判断结局，完全数据/截尾数据
# # survdata[1:10,1:2] #展示分类变量的生存数据
# # KMfit <- survfit(survdata ~ meta$group)   #~ 后是指定分组
# # head(survdata)
# # ggsurvplot(KMfit,                       #拟合对象           
# #            data = matrix,               #变量数据来源           
# #            pval = TRUE,                 #P值           
# #            surv.median.line = "hv",     #中位生存时间线           
# #            conf.int = TRUE,           
# #            risk.table = TRUE,           #风险表           
# #            xlab = "Follow up time(months)",  #x轴标签           
# #            break.x.by = 10)             #x轴刻度间距 
# 
# colnames(lncrna_expr_count) <- str_sub(colnames(lncrna_expr_count), 1, 16)
# colmeans <- function(x){
#   exp_m <- as.matrix(lncrna_expr_count)
#   exp_t <- t(exp_m)
#   exp_t <- limma::avereps(exp_t)
#   t(exp_t)
# }
# lncrna_expr_count <- colmeans(lncrna_expr_count)
# 
# meta <- select(matrix, c("ID","month","Status"))  #matrix中614个样本
# meta <- as.data.frame(meta)
# meta <- limma::avereps(meta, meta$ID) #取平均取出重复临床样本名称的数据，meta中533个样本
# meta <- as.data.frame(meta)
# rownames(meta) <- meta$ID
# 
# # 表达数据和临床数据匹配，找出既有表达数据，又有生存数据的样本
# index <- rownames(meta)[rownames(meta) %in% colnames(lncrna_expr_count)]
# meta <- meta[index, ]
# 
# Gene = "LINC01122"
# meta$Expression <- lncrna_expr_count[Gene, rownames(meta)]
# 
# copy <- meta
# copy <- apply(copy, 2, as.numeric)
# meta <- meta[, c("Status", "Expression")]
# meta2 <- cbind(meta, copy)
# 
# meta <- select(meta2, c("Status","Expression","month"))
# 
# ## 使用cutoff包的logrank函数来寻找最佳阈值
# meta_for_cutoff <- meta
# meta_for_cutoff$Status <- ifelse(meta_for_cutoff$Status == "Dead", 1, 0)
# best_threshold_surv <- logrank(data = meta_for_cutoff, 
#                                time = "month", 
#                                y = "Status",
#                                x = "Expression", 
#                                cut.numb = 1, 
#                                n.per = 0.2, 
#                                y.per = 0.1, 
#                                p.cut = 0.1, 
#                                round = 5)
# print(best_threshold_surv)
# 
# # 根据找到的最佳阈值对数据进行分组
# meta$Expression <- ifelse(lncrna_expr_count[Gene, rownames(meta)] > best_threshold_surv$cut1, "High", "Low")
# 
# 
# # 创建生存对象（Surv对象），包含生存时间和生存状态信息
# survdata <- Surv(time = meta$month, 
#                 event = meta$Status == "Dead")
# 
# survdata[1:10, 1:2]
# KMfit <- survfit(survdata ~ meta$Expression)
# 
# ggsurvplot(KMfit, 
#            data = meta, 
#            pval = TRUE, 
#            conf.int = TRUE, 
#            risk.table = TRUE, 
#            xlab = "Follow up time(months)",  #x轴标签
#            break.x.by = 10,             #x轴刻度间距
#            #legend = c(0.8,0.75),       #图例位置
#            legend.title = "LINC00461",          #图例标题
#            legend.labs = c("High", "Low"))

# TCGA LGG ----------------------------------------------------------------

## survival analysis
load("D:/R_project/UCR_project/01-data/45-GBM/output_mRNA_lncRNA_expr/TCGA-LGG_clinicalSE.rdata")
load("D:/R_project/UCR_project/01-data/45-GBM/output_mRNA_lncRNA_expr/TCGA-LGG_lncrna_expr_count.rdata")

# 准备临床数据
LGG_clinical <- clinicalSE[, c("days_to_last_follow_up", "vital_status")]
names(LGG_clinical) <- c("time", "event")
LGG_clinical$event <- ifelse(LGG_clinical$event=="Dead", 1, 0)

pdf("03-results/45-GBM/surv_LINC01019.pdf", width = 633/254, height = 760/254)

marker = "LINC01117"

res <- plot_KM(exprset = lncrna_expr_count, 
               marker = marker, 
               clin = GBM_clinical, 
               optimal_cut = F, return_data = TRUE)

dev.off()

# GSE147352 LGG -----------------------------------------------------------

# # NCUCR genes in LGG
# load("D:/R_project/UCR_project/01-data/45-GBM/output_mRNA_lncRNA_expr/TCGA-LGG_lncrna_expr_count.rdata")
# 
# # convert ensembl id in website
# converted_NCUCR_genes <- read.table(file = "02-analysis/45-GBM/NCUCR_genes_id.txt", 
#                                     sep = "\t", blank.lines.skip = TRUE)
# # remove null value
# converted_NCUCR_genes <- converted_NCUCR_genes[!(converted_NCUCR_genes$V2 == ""), ]
# 
# converted_NCUCR_genes_exp_LGG <- subset(lncrna_expr_count, 
#                                         rownames(lncrna_expr_count) %in% converted_NCUCR_genes$V2)
# 
# LGG_diff_analysis_results <- diff_analysis(exprset = lncrna_expr_count, 
#                                            project = "LGG", 
#                                            is_count = TRUE)

library(data.table)
library(purrr)
library(DESeq2)
library(openxlsx)

## GSE147352 ReadsCount
# FileNames <- list.files(path = "01-data/45-GBM/GSE147352_RAW/", pattern = "*.txt.gz")
# FileNames <- paste("D:/R_project/UCR_project/01-data/45-GBM/GSE147352_RAW/", FileNames, sep = "")
# DataList <- lapply(FileNames, function(file){
#   # read data
#   data <- fread(file, header = FALSE)
#   
#   # 数据框第一行是基因ID，第二行是表达量
#   colnames(data) <- c("GeneID", gsub("\\.txt\\.gz$", "", file))  # 设置列名
#   return(data)
# })
# 
# # 合并所有数据根据 GeneID
# merged_data <- Reduce(function(x, y) merge(x, y, by = "GeneID", all = TRUE), DataList)
# merged_data <- merged_data[1:58302, ]
# 
# # https://www.jianshu.com/p/bf5bce8a4b1d
# m <- as.matrix(merged_data[, 1])
# v <- as.vector(m)
# merged_data <- merged_data[, -1]
# colnames(merged_data) <- gsub("D:/R_project/UCR_project/01-data/45-GBM/GSE147352_RAW/", "", 
#                               colnames(merged_data))
# merged_data <- as.data.frame(merged_data)
# rownames(merged_data) <- v
# 
# ReadsCount <- merged_data
# save(ReadsCount, file = "01-data/45-GBM/GSE147352_ReadsCount.rdata")

# get ReadsCount
load("D:/R_project/UCR_project/01-data/45-GBM/GSE147352_ReadsCount.rdata")

# get sample information
SampleInfo <- read.xlsx(xlsxFile = "01-data/45-GBM/GSE147352_clininfo.xlsx", sheet = 1)
table(SampleInfo$Type)
# GBM 85 LGG 18 Normal 15
GBMSampleInfo <- filter(SampleInfo, Type == "Normal" | Type == "GBM")
LGGSampleInfo <- filter(SampleInfo, Type == "Normal" | Type == "LGG")

# make ReadsCountTable
LGGReadsCount <- select(ReadsCount, LGGSampleInfo$sample.ID)
LGGReadsCount <- LGGReadsCount[rowSums(LGGReadsCount) > 33, ]
LGGReadsCount <- as.matrix(LGGReadsCount)

# make condition group
Condition <- factor(LGGSampleInfo$Type, levels = c("Normal", "LGG"))
Coldata <- data.frame(row.names = colnames(LGGReadsCount), Condition)

dds <- DESeqDataSetFromMatrix(countData = LGGReadsCount, 
                              colData = Coldata, 
                              design = ~ Condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("Condition", "LGG", "Normal"))
resdata <- as.data.frame(res)

# get NCUCR genes
NCUCR_genes <- read.table(file = "02-analysis/16-New_classification/02-UCR_from_Non_Coding_RNA(66).bed")
NCUCR_genes <- NCUCR_genes[, c(8, 9)]
NCUCR_genes <- unique(NCUCR_genes)
print(noquote(NCUCR_genes[["V8"]]))

LGGDENCUCRGenes <- subset(resdata, rownames(resdata) %in% NCUCR_genes$V8)
LGGDENCUCRGenes$type <- ifelse(LGGDENCUCRGenes$padj > 0.05 | abs(LGGDENCUCRGenes$log2FoldChange) < 1, "nosig", 
                           ifelse(LGGDENCUCRGenes$log2FoldChange > 1, "up", "down"))

LGGMarker <- LGGDENCUCRGenes[!(LGGDENCUCRGenes$type %in% "nosig"), ]
LGGMarker$geneid <- rownames(LGGMarker)
LGGMarker <- merge(LGGMarker, NCUCR_genes, 
                   by.x = "geneid", by.y = "V8")
LGGMarker$V9 <- ifelse(LGGMarker$V9 == ".", LGGMarker$geneid, LGGMarker$V9)

lwd_pt <- .pt*72.27/96

font_add(family = "Arial", regular = "C:/Windows/Fonts/arial.ttf")
font_families()

showtext_begin()

NCUCR_genes_volcano_LGG <- ggplot(data = LGGDENCUCRGenes) +
  
  geom_point(mapping = aes(x = log2FoldChange, y = -log10(padj), 
                           fill = type, color = type), size = 0.8) +
  
  scale_fill_manual(values = c("blue", "grey", "red")) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  
  scale_x_continuous(limits = c(-5, 5), expand = c(0, 0)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5/lwd_pt) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.5/lwd_pt) +
  
  geom_text_repel(
    data = subset(LGGDENCUCRGenes, rownames(LGGDENCUCRGenes) %in% LGGMarker$geneid), 
    aes(x = log2FoldChange, y = -log10(padj), label = LGGMarker$V9), 
    alpha = 0.5, 
    max.overlaps = 15, 
    size = unit(1.5, "pt")
  ) +
  
  theme_prism(base_size = 7, 
              base_family = "Arial", 
              base_fontface = "plain", 
              base_line_size = 0.5/lwd_pt, 
              base_rect_size = 0.5/lwd_pt, 
              border = TRUE) +
  
  theme(
    legend.position = "none", 
    legend.key.size = unit(0.25, "cm"), 
    aspect.ratio = 1
  )

NCUCR_genes_volcano_LGG

ggview::ggview(NCUCR_genes_volcano_LGG, 
               width = 4.8, height = 6, units = "cm", dpi = 1200)

# GSE147352 GBM -----------------------------------------------------------

# get sample information
SampleInfo <- read_xlsx(path = "01-data/45-GBM/GSE147352_clininfo.xlsx", sheet = 1)
table(SampleInfo$Type)
# GBM 85 LGG 18 Normal 15
GBMSampleInfo <- filter(SampleInfo, Type == "Normal" | Type == "GBM")
LGGSampleInfo <- filter(SampleInfo, Type == "Normal" | Type == "LGG")

# make ReadsCountTable
GBMReadsCount <- select(ReadsCount, GBMSampleInfo$`sample ID`)
GBMReadsCount <- GBMReadsCount[rowSums(GBMReadsCount) > 100, ]
GBMReadsCount <- as.matrix(GBMReadsCount)

# make condition group
Condition <- factor(GBMSampleInfo$Type, levels = c("Normal", "GBM"))
Coldata <- data.frame(row.names = colnames(GBMReadsCount), Condition)

dds <- DESeqDataSetFromMatrix(countData = GBMReadsCount, 
                              colData = Coldata, 
                              design = ~ Condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("Condition", "GBM", "Normal"))
resdata <- as.data.frame(res)

# get NCUCR genes
NCUCR_genes <- read.table(file = "02-analysis/16-New_classification/02-UCR_from_Non_Coding_RNA(66).bed")
NCUCR_genes <- NCUCR_genes[, c(8, 9)]
NCUCR_genes <- unique(NCUCR_genes)
print(noquote(NCUCR_genes[["V8"]]))

GBMDENCUCRGenes <- subset(resdata, rownames(resdata) %in% NCUCR_genes$V8)
GBMDENCUCRGenes$type <- ifelse(GBMDENCUCRGenes$padj > 0.05 | abs(GBMDENCUCRGenes$log2FoldChange) < 1, "nosig", 
                               ifelse(GBMDENCUCRGenes$log2FoldChange > 1, "up", "down"))
write.table(GBMDENCUCRGenes, file = "03-results/45-GBM/GSE147352_GBM_DEG.txt", 
            sep = "\t", quote = FALSE)

GBMMarker <- GBMDENCUCRGenes[!(GBMDENCUCRGenes$type %in% "nosig"), ]
GBMMarker$geneid <- rownames(GBMMarker)
GBMMarker <- merge(GBMMarker, NCUCR_genes, 
                   by.x = "geneid", by.y = "V8")
GBMMarker$V9 <- ifelse(GBMMarker$V9 == ".", GBMMarker$geneid, GBMMarker$V9)

lwd_pt <- .pt*72.27/96

font_add(family = "Arial", regular = "C:/Windows/Fonts/arial.ttf")
font_families()

showtext_begin()

NCUCR_genes_volcano_GEOGBM <- ggplot(data = GBMDENCUCRGenes) +
  
  geom_point(mapping = aes(x = log2FoldChange, y = -log10(padj), 
                           fill = type, color = type), size = 0.8) +
  
  scale_fill_manual(values = c("blue", "grey", "red")) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  
  scale_x_continuous(limits = c(-5, 5), expand = c(0, 0)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5/lwd_pt) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.5/lwd_pt) +
  
  geom_text_repel(
    data = subset(GBMDENCUCRGenes, rownames(GBMDENCUCRGenes) %in% GBMMarker$geneid), 
    aes(x = log2FoldChange, y = -log10(padj), label = GBMMarker$V9), 
    alpha = 0.5, 
    max.overlaps = 15, 
    size = unit(1.5, "pt")
  ) +
  
  theme_prism(base_size = 7, 
              base_family = "Arial", 
              base_fontface = "plain", 
              base_line_size = 0.5/lwd_pt, 
              base_rect_size = 0.5/lwd_pt, 
              border = TRUE) +
  
  theme(
    legend.position = "none", 
    legend.key.size = unit(0.25, "cm"), 
    aspect.ratio = 1
  )

NCUCR_genes_volcano_GEOGBM

showtext_end()

# merge 3 volcano ---------------------------------------------------------

pdf(file = "03-results/45-GBM/NCUCR_genes_volcano_TCGA_GSE147352.pdf", width = 1500/254, height = 600/254)
showtext_begin()

cowplot::plot_grid(NCUCR_genes_volcano_GBM, NCUCR_genes_volcano_GEOGBM, NCUCR_genes_volcano_LGG, 
                   align = c("hv"), axis = c("tb"), 
                   nrow = 1, ncol = 3, rel_widths = 1, rel_heights = 1, 
                   labels = c("a", "b", "c"), label_size = 8)

showtext_end()
dev.off()












