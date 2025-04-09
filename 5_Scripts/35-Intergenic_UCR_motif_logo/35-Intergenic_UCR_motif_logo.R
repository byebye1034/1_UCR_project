# 绘制MEME在intergenic UCR当中找到的富集的motif logo

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(seqLogo)
library(tidyverse)

# CTCF <- read.table(file = "01-data/11.1-Get_intergenic_UCR_closest_pc_gene/MA0139.1.jaspar", 
#                    header = FALSE)
# CTCF <- CTCF[, c(3:21)]
# rownames(CTCF) <- c("A", "C", "G", "T")
# 
# # 转为PWM矩阵
# # 使用apply()函数对每列进行操作
# CTCF_1 <- apply(CTCF, 2, FUN = function(x) x / sum(x))
# CTCF_1 <- data.frame(CTCF_1)
# 
# trans_ctcf <- makePWM(CTCF_1)
# seqLogo(trans_ctcf)

# MEME found a motif enriched in intergenic UCR
intergenic_motif <- read.table(
  file = "01-data/11.1-Get_intergenic_UCR_closest_pc_gene/motif_intregenic_UCR_MEME_found.txt", 
  sep = " "
)
colnames(intergenic_motif) <- c("A", "C", "G", "T")
intergenic_motif <- data.frame(t(intergenic_motif))

intergenic_motif <- makePWM(intergenic_motif)
seqLogo(intergenic_motif, 
        xfontsize = 7, yfontsize = 7)

pdf(
  file = "03-results/35-Intergenic_UCR_motif_logo/intergenic_motif.pdf", 
  width = 1100/254, height = 500/254
)

seqLogo(intergenic_motif, 
        xfontsize = 7, yfontsize = 7)

dev.off()





