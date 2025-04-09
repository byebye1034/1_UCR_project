# 绘制染色体核型图，展示479个UCR的具体位置和类型

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(karyoploteR)
library(GenomicRanges) # 将数据框转换成GRanges对象
library(tidyverse) # 没有这个不能使用lwd的设置
library(chromoMap)

UCR_location <- read.table(file = "01-data/UCR_raw/UCR_location.txt", sep = "\t", 
                           header = T)
# UCR_location$chr_name <- paste("chr", UCR_location$chr_name, sep = "")

# UCR分类信息
# coding UCR -----------------------------------------------------------------

coding_UCR <- read.table(file = "02-analysis/16-New_classification/01-UCR_from_protein_coding_gene.bed", 
                         sep = "\t", header = FALSE)
coding_UCR$V9 <- ifelse(coding_UCR$V9 == ".", 
                        coding_UCR$V8, 
                        coding_UCR$V9)

coding_UCR_intergrated <- aggregate(V9 ~ V4, 
                                    data = coding_UCR, 
                                    FUN = paste, 
                                    collapse = "/")

coding_UCR_intergrated <- coding_UCR_intergrated %>% 
  mutate(type = "coding_UCR")
colnames(coding_UCR_intergrated) <- c("UCR_name", "Gene_name", "UCR_type")

# ncRNA UCR ------------------------------------------------------------------

ncRNA_UCR <- read.table(file = "02-analysis/16-New_classification/02-UCR_from_Non_Coding_RNA(66).bed", 
                        sep = "\t", header = FALSE)
ncRNA_UCR$V9 <- ifelse(ncRNA_UCR$V9 == ".", 
                    ncRNA_UCR$V8, 
                    ncRNA_UCR$V9)

ncRNA_UCR_intergrated <- aggregate(V9 ~ V4,  # 对前面的进行分类，分类的原则是后面的
                                   data = ncRNA_UCR, 
                                   FUN = paste, 
                                   collapse = "/")

ncRNA_UCR_intergrated <- ncRNA_UCR_intergrated %>% 
  mutate(type = "ncRNA_UCR")
colnames(ncRNA_UCR_intergrated) <- c("UCR_name", "Gene_name", "UCR_type")

# intergenic UCR ----------------------------------------------------------

# intergenic
intergenic_UCR <- read.table(file = "02-analysis/16-New_classification/03-UCR_intergenic.bed", 
                             sep = "\t", header = FALSE)
intergenic_UCR <- intergenic_UCR %>% 
  mutate(Gene_name = "intergenic") %>% 
  mutate(UCR_type = "intergenic")
intergenic_UCR <- intergenic_UCR[, c(4:6)]
colnames(intergenic_UCR)[1] <- "UCR_name"

# UCR_type_gene -----------------------------------------------------------

# coding UCR 314
# ncRNA UCR 66
# intergenic UCR 99

UCR_type_gene <- rbind(coding_UCR_intergrated, ncRNA_UCR_intergrated, intergenic_UCR)
save(UCR_type_gene, file = "02-analysis/12-karyoplote/UCR_type_gene.Rdata")

# UCR_location_plot -------------------------------------------------------

rm(list = ls())

load("D:/R_project/UCR_project/02-analysis/12-karyoplote/UCR_type_gene.Rdata")
UCR_location <- read.table(file = "01-data/UCR_raw/UCR_location.txt", sep = "\t", 
                           header = T)
UCR_location <- UCR_location[, c(1:4)]

UCR_feature <- merge(UCR_location[, c(4, 1:3)], UCR_type_gene, by = "UCR_name")
UCR_feature$chr_name <- paste("chr", UCR_feature$chr_name, sep = "")

UCR_GRanges <- GRanges(seqnames = UCR_feature$chr_name, 
                       ranges = IRanges(start = UCR_feature$start, end = UCR_feature$end))

# 将类别添加到GRanges对象的元数据中
mcols(UCR_GRanges)$UCR_type <- UCR_feature$UCR_type

coding_UCR_gene_label <- UCR_GRanges[mcols(UCR_GRanges)$UCR_type == "coding_UCR"]
ncRNA_UCR_gene_label <- UCR_GRanges[mcols(UCR_GRanges)$UCR_type == "ncRNA_UCR"]

# 将基因添加到GRanges对象的元数据中
mcols(UCR_GRanges)$gene <- UCR_feature$Gene_symbol

# 创建一个颜色映射，将类别映射到颜色
col_mapping <- c("coding_UCR" = "#E64B35FF", 
                 "ncRNA_UCR" = "#4DBBD5FF",
                 "intergenic" = "#3C5488FF")

lwd_pt <- .pt*72.27/96

# 根据需要调整边距和字符大小
params$mar <- c(3, 3, 3, 3)  # 这里可以根据需要调整数值
params$cex <- 0.8  # 字符大小也可以根据需要调整

# 绘制kp
kp <- plotKaryotype(genome = "hg38", plot.type = 1)

# 在karyoplot上绘制区域
kpPlotRegions(kp, data = UCR_GRanges, 
              col = col_mapping[mcols(UCR_GRanges)$UCR_type], 
              lwd = 0.5/lwd_pt)

# 添加基因名
kpText(kp, data = coding_UCR_gene_label, labels = mcols(coding_UCR_gene_label)$gene, r0 = 0.1, r1 = 0.2, cex = 0.4, y = 1)
kpText(kp, data = ncRNA_UCR_gene_label, labels = mcols(ncRNA_UCR_gene_label)$gene, r0 = 0.2, r1 = 0.4, cex = 0.4, y = 1)

pdf(file = "03-results/12-karyoploeplot/UCR_location_plot.pdf", width = 17, height = 14, unit = "cm")

dev.off()



























