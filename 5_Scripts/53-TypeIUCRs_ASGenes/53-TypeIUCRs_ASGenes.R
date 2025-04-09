# the overlap relationship of type I UCRs and as genes
# get the bed file of exonic UCRs, intronic UCRs, AS genes and chr length
# plot circle plot by TB-tools

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library(readxl)
library(circlize)

lwd_pt <- .pt*72.27/96

## exonic UCRs bed
exonicUCRs <- read.table(
  file = "02-analysis/50-TypeIUCRs_PCG_Location/exonicUCRs.bed", 
  sep = "\t"
)
exonicUCRs <- exonicUCRs[, c(1:4)]
exonicUCRs <- unique(exonicUCRs)
write.table(exonicUCRs, 
            file = "03-results/53-TypeIUCRs_ASGenes/exonicUCRs.bed", 
            col.names = FALSE, row.names = FALSE, 
            sep = "\t", quote = FALSE)

## intronic UCRs bed
typeIUCRs <- read.table(
  file = "02-analysis/16-New_classification/01-UCR_from_protein_coding_gene.bed"
)
typeIUCRs <- unique(typeIUCRs[, c(1:4)])
intronicUCRs <- anti_join(typeIUCRs, exonicUCRs, by = "V4")
write.table(intronicUCRs, 
            file = "03-results/53-TypeIUCRs_ASGenes/intronicUCRs.bed", 
            col.names = FALSE, row.names = FALSE, 
            sep = "\t", quote = FALSE)

## alternative splicing related genes bed
typeIUCRs <- read.table(
  file = "02-analysis/16-New_classification/01-UCR_from_protein_coding_gene.bed"
)
typeIGenes <- typeIUCRs[, c(5:9)]
typeIGenes <- unique(typeIGenes)
write.table(typeIGenes, 
            file = "03-results/53-TypeIUCRs_ASGenes/typeIGenes.bed", 
            col.names = FALSE, row.names = FALSE, 
            sep = "\t", quote = FALSE)
write.table(typeIGenes[, 4], 
            file = "03-results/53-TypeIUCRs_ASGenes/typeIGenesID.bed", 
            col.names = FALSE, row.names = FALSE, 
            sep = "\t", quote = FALSE)

typeIGenesGO <- read_xlsx(
  path = "03-results/53-TypeIUCRs_ASGenes/typeIGenesGO/metascape_result.xlsx", 
  sheet = 2)
typeIGenesGO <- typeIGenesGO[str_detect(typeIGenesGO$GroupID, ".*Summary"), ]
AS1 <- unlist(str_split(typeIGenesGO[1, 8], ","))
AS2 <- unlist(str_split(typeIGenesGO[4, 8], ","))
ASGenes <- unique(c(AS1, AS2))
ASGenes[36] <- "FAM172A"
ASGenesbed <- typeIGenes[typeIGenes$V9 %in% ASGenes, ]
write.table(ASGenesbed, 
            file = "03-results/53-TypeIUCRs_ASGenes/ASGenes.bed", 
            col.names = FALSE, row.names = FALSE, 
            sep = "\t", quote = FALSE)

# get all AS Genes and coresponding UCR
ASOverlapUCR <- typeIUCRs[typeIUCRs$V9 %in% ASGenes, ]
write.csv(ASOverlapUCR, 
          file = "03-results/53-TypeIUCRs_ASGenes/ASGenesOverlapUCR.csv", 
          row.names = FALSE)

## chromosome length
chromosomeInfo <- read.table(
  file = "02-analysis/03-SNP_INFO_Retrieval/GRCh38_p14_sequence_report.tsv", 
  sep = "\t", header = TRUE
)
chromosomeInfo <- chromosomeInfo[c(1:24), c(4, 12)]
chromosomeInfo$Seq.length <- as.numeric(chromosomeInfo$Seq.length)
write.table(chromosomeInfo, 
            file = "03-results/53-TypeIUCRs_ASGenes/chrLength.txt", 
            col.names = FALSE, row.names = FALSE, 
            sep = "\t", quote = FALSE)

######################################################
#                STEP 1 - 数据预处理                 #
######################################################
# 处理染色体信息 (确保列名对应)
genome <- chromosomeInfo %>%
  rename(chr = Chromosome.name, length = Seq.length)

# 处理exonicUCRs (列名标准化)
exonic_ucr <- exonicUCRs %>%
  select(chr = V1, start = V2, end = V3, name = V4)

# 处理intronicUCRs (列名标准化)
intronic_ucr <- intronicUCRs %>%
  select(chr = V1, start = V2, end = V3, name = V4)

# 处理ASGenesbed (提取关键列)
as_genes <- ASGenesbed %>%
  select(chr = V5, start = V6, end = V7, name = V9) %>%
  mutate(chr = as.character(chr))  # 确保染色体为字符型

# 设置绘图参数
circos.par(
  gap.degree = 2,
  cell.padding = c(0, 0, 0, 0),
  track.margin = c(0, 0.02)
)

# 初始化染色体
circos.initialize(
  factors = genome$chr,
  xlim = cbind(rep(0, nrow(genome)), genome$length)
)

# 绘制染色体基线（中线）
circos.track(
  ylim = c(0, 1),
  panel.fun = function(x, y) {
    # 绘制染色体中线
    circos.lines(CELL_META$cell.xlim, c(0.5, 0.5), 
                 col = "gray30", lwd = 1.5)
  },
  track.height = 0.05,
  bg.border = NA
)

# 自定义绘制垂直线段的函数
draw_segments <- function(data, color, y_center = 0.5, 
                          segment_length = 0.3) {
  circos.genomicTrack(
    data, bg.border = NA, 
    ylim = c(0, 1),
    panel.fun = function(region, value, ...) {
      # 计算中点位置
      x_center <- (region$start + region$end) / 2
      # 绘制垂直线段
      circos.segments(
        x0 = x_center, y0 = y_center - segment_length/2,
        x1 = x_center, y1 = y_center + segment_length/2,
        col = color,
        lwd = 1.5/lwd_pt
      )
    },
    track.height = 0.15
  )
}

# 需要标注的基因
highlight_genes <- c("SRSF1", "SRSF3", "SRSF6", "SRSF7", "SRSF11", 
                     "HNRNPU", "HNRNPA1L3", "HNRNPM", "HNRNPH1", "HNRNPK")

# 筛选出需要标注的基因
highlight_as_genes <- as_genes %>% filter(name %in% highlight_genes)

# 绘制SR和HNRNP蛋白家族基因
draw_segments(highlight_as_genes, "#DC0000FF", y_center = 1.6, segment_length = 1)

# 绘制AS基因（上方线段）
draw_segments(as_genes, "#3C5488FF", y_center = 1.2, segment_length = 0.4)

# 绘制exonic UCRs（中线位置）
draw_segments(exonic_ucr, "#E64B35FF", y_center = 1.1, segment_length = 0.4)

# 绘制intronic UCRs（下方线段）
draw_segments(intronic_ucr, "#4DBBD5FF", y_center = 1.0, segment_length = 0.4)

# # 在圈图上添加基因名称
# circos.genomicTrack(
#   highlight_as_genes, bg.border = NA,
#   ylim = c(0, 1),
#   panel.fun = function(region, value, ...) {
#     x_center <- (region$start + region$end) / 2
#     circos.text(x = x_center, y = 0.8, labels = value$name,
#                 facing = "clockwise", adj = c(0, 0.5),
#                 cex = 0.7, col = "black", font = 2/lwd_pt)
#   }
# )

# 添加图例
legend(
  x = "bottomleft",
  legend = c("AS Genes", "Exonic UCRs", "Intronic UCRs"),
  col = c("#3C5488FF", "#E64B35FF", "#4DBBD5FF"),
  lwd = 2,
  bty = "n",
  cex = 1.2/lwd_pt,
  title = NULL
)

circos.clear()

pdf(file = "03-results/53-TypeIUCRs_ASGenes/circle.pdf", 
    width = 900/254, height = 900/254)

dev.off()
































