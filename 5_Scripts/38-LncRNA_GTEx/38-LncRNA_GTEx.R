# 先在GTEx里查看lncRNA在哪些组织里表达：可能是大脑和生殖系统

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(pheatmap)

# ncRNA UCR重叠的lncRNA
ncRNA_UCR_lncRNA <- read.table(file = "02-analysis/16-New_classification/ncRNA_UCR(66)_ensembl_id.txt", sep = "\t")
ncRNA_UCR_lncRNA <- unique(ncRNA_UCR_lncRNA)

# GTEx data
GTEx_exp <- read.table(file = "01-data/22-Tissue_specificity_across_different_organs/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", 
                       sep = "\t", skip = 2, header = TRUE)
GTEx_exp$Name <- str_sub(GTEx_exp$Name, 1, 15)
GTEx_exp <- GTEx_exp[, -2]

ncRNA_UCR_lncRNA_exp <- merge(ncRNA_UCR_lncRNA, 
                              GTEx_exp, 
                              by.x = "V1", 
                              by.y = "Name")
rownames(ncRNA_UCR_lncRNA_exp) <- ncRNA_UCR_lncRNA_exp$V1
ncRNA_UCR_lncRNA_exp <- ncRNA_UCR_lncRNA_exp[, -1]

# 过滤低表达部分
zero_count <- rowSums(ncRNA_UCR_lncRNA_exp == 0)
ncRNA_UCR_lncRNA_exp <- ncRNA_UCR_lncRNA_exp %>% 
  filter(zero_count < ncol(ncRNA_UCR_lncRNA_exp) / 2)

ncRNA_UCR_lncRNA_exp <- log2(ncRNA_UCR_lncRNA_exp + 1)
ncRNA_UCR_lncRNA_exp <- round(ncRNA_UCR_lncRNA_exp, 2)

colnames(ncRNA_UCR_lncRNA_exp) <- factor(x = colnames(ncRNA_UCR_lncRNA_exp), 
                                         levels = colnames(ncRNA_UCR_lncRNA_exp))

max(ncRNA_UCR_lncRNA_exp)
# 4.66
min(ncRNA_UCR_lncRNA_exp)
# 0

# breaks
bk <- c(seq(0, 2.4, by = 0.01), seq(2.5, 5, by = 0.01))

# 绘图
pdf(file = "03-results/38-LncRNA_GTEx/lncRNA_GTEx_heatmap.pdf", width = 10, height = 15)
pheatmap(mat = ncRNA_UCR_lncRNA_exp, 
         cellwidth = 10, cellheight = 10, 
         cluster_rows = T, cluster_cols = T, 
         angle_col = 90, main = "log2(exp+1)", 
         color = c(colorRampPalette(colors = c("#083490","#f1f1f1"))(length(bk)/2), 
                   colorRampPalette(colors = c("#f1f1f1","#E60212"))(length(bk)/2)),
         legend_breaks=seq(0, 5, 1), 
         breaks = bk, 
         fontsize = 8, 
         display_numbers = FALSE)
dev.off()

install.packages("LncPath")
library(LncPath)

# get lncRNA-mRNA interaction network
NetLncPath <- getNet()

print(head(NetLncPath), row.names = FALSE)

# get lncRNA sets
# 使用paste函数将一列数据框中的值连接成一个字符串
ncRNA_UCR_lncRNA_strings <- as.character(unlist(ncRNA_UCR_lncRNA))
print(head(ncRNA_UCR_lncRNA_strings))

# get example lncRNA-mRNA interaction network
ExampleNet <- getExampleData("ExampleNet")

Result <- lncPath(LncRNAList = ncRNA_UCR_lncRNA_strings, 
                  Network = NetLncPath, 
                  Weighted = TRUE, 
                  PathwayDataSet = "KEGG", 
                  nperm = 1000, 
                  minPathSize = 15, 
                  maxPathSize = 500)

save(Result, file = "03-results/38-LncRNA_GTEx/LncPath_result.Rdata")

# Now start the random walking...
# [1] "ENSG00000228956" "ENSG00000224243" "ENSG00000245526" "ENSG00000247828" "ENSG00000233723" "ENSG00000271860"
# [7] "ENSG00000234377" "ENSG00000224592" "ENSG00000248079" "ENSG00000234350" "ENSG00000234323" "ENSG00000228835"
# [13] "ENSG00000224577" "ENSG00000243620" "ENSG00000231764" "ENSG00000248118" "ENSG00000258038" "ENSG00000235724"
# [19] "ENSG00000227290" "ENSG00000257585" "ENSG00000249483" "ENSG00000233891" "ENSG00000257523" "ENSG00000254946"
# [25] "ENSG00000258107" "ENSG00000257126"

# Generate a table of the summary of each pathway
PathwaySummaryTable <- lncPath2Table(Result)
PathwaySummaryTable <- filter(PathwaySummaryTable, PathwaySummaryTable$`False Discovery Rate` < 0.05)

write.csv(PathwaySummaryTable, 
          file = "03-results/38-LncRNA_GTEx/PathwaySummaryTable.csv", row.names = FALSE)

# MAPK signaling pathway
# Calcium signaling pathway
# Chemokine signaling pathway
# Neuroactive ligand-receptor interaction
# Cardiac muscle contraction
# Vascular smooth muscle contraction
# Adherens junction
# Gap junction
# Neurotrophin signaling pathway
# Taste transduction
# Regulation of actin cytoskeleton
# GnRH signaling pathway
# Pathways in cancer
# Hypertrophic cardiomyopathy (HCM)
# Arrhythmogenic right ventricular cardiomyopathy (ARVC)
# Dilated cardiomyopathy
# ErbB signaling pathway
# Wnt signaling pathway
# Long-term potentiation
# Long-term depression
# Melanogenesis
# Vasopressin-regulated water reabsorption
# Focal adhesion
# Leukocyte transendothelial migration
# Ribosome
# Phosphatidylinositol signaling system
# Fc epsilon RI signaling pathway
# Glioma
# Chronic myeloid leukemia
# Prostate cancer

ENSG00000228956
ENSG00000224243
ENSG00000245526
ENSG00000247828
ENSG00000233723
ENSG00000271860
ENSG00000234377
ENSG00000224592
ENSG00000248079
ENSG00000234350
ENSG00000234323
ENSG00000228835
ENSG00000224577
ENSG00000243620
ENSG00000231764
ENSG00000248118
ENSG00000258038
ENSG00000235724
ENSG00000227290
ENSG00000257585
ENSG00000249483
ENSG00000233891
ENSG00000257523
ENSG00000254946
ENSG00000258107
ENSG00000257126

















