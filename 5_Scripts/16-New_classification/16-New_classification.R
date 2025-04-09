# 新分类和分类饼图

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(dbplyr)
library(tidyverse)

UCR_coding <- read.table(file = "02-analysis/16-New_classification/01-UCR_from_protein_coding_gene.bed", 
                         sep = "\t")
write.table(UCR_coding$V8, file = "02-analysis/16-New_classification/coding_UCR_ensembl_id.txt", 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# 整理01.1-exonic UCR
exonic_UCR <- read.table(file = "02-analysis/16-New_classification/01.1-UCR_from_protein_coding_gene_exon.bed", 
                         sep = "\t")
write.table(unique(exonic_UCR$V8), 
            file = "02-analysis/16-New_classification/exonic_UCR_ensembl_id.txt", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE) # 126个exonic UCR对应111个coding gene
exonic_UCR <- exonic_UCR[, c(1:4)]
exonic_UCR <- unique(exonic_UCR)
# exonic UCR:126

# 获得01.2-intronic UCR
temp <- UCR_coding[, c(1:4)]
intronic_UCR <- setdiff(temp, exonic_UCR)

# intronic UCR:188
temp <- merge(intronic_UCR, UCR_coding, by = "V4")
temp <- temp[, c(2:4, 1, 11, 12)]
temp <- unique(temp)
write.table(unique(temp$V8), 
            file = "02-analysis/16-New_classification/intronic_UCR_ensembl_id.txt", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE) # 188个intronic UCR对应99个coding gene
write.table(temp, 
            file = "02-analysis/16-New_classification/01.2-UCR_from_protein_coding_gene_intron.bed", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# 将ncRNA当中和coding基因有交集的部分去除
UCR_ncRNA <- read.table(file = "02-analysis/16-New_classification/02-UCR_from_Non_Coding_RNA.bed", 
                        sep = "\t")
write.table(UCR_ncRNA$V8, file = "02-analysis/16-New_classification/ncRNA_UCR(100)_ensembl_id.txt", 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

UCR_ncRNA <- UCR_ncRNA[!(UCR_ncRNA$V4 %in% UCR_coding$V4), ]
write.table(UCR_ncRNA$V8, file = "02-analysis/16-New_classification/ncRNA_UCR(66)_ensembl_id.txt", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(UCR_ncRNA, file = "02-analysis/16-New_classification/02-UCR_from_Non_Coding_RNA(66).bed", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# 绘制分类饼图 ------------------------------------------------------------------

# 按照新思路对UCR进行分类
# 1、coding_UCR 314 exonic_UCR 126 intronic_UCR 188
# 2、ncRNA_UCR 100(去除和1重叠的还有66个)
# 3、intergenic_UCR 99
# 1和2重叠的有34个
# [1] "uc.43"  "uc.325" "uc.326" "uc.327" "uc.338" "uc.343" "uc.344" "uc.368" "uc.389" "uc.408" "uc.415" "uc.416"
# [13] "uc.417" "uc.420" "uc.427" "uc.64"  "uc.69"  "uc.70"  "uc.71"  "uc.72"  "uc.100" "uc.136" "uc.142" "uc.167"
# [25] "uc.169" "uc.213" "uc.217" "uc.222" "uc.223" "uc.232" "uc.263" "uc.264" "uc.472" "uc.477"
# note:在intergenic当中uc.329和一个假基因重叠
# 11	32176446	32176752	uc.329	11	32112049	32343409	ENSG00000227160	THEM7P	transcribed_unitary_pseudogene

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(dbplyr)
library(ggplot2)
library(ggsci)
library(ggview)
library(ggprism)
library(sysfonts)
library(showtext)

lwd_pt <- .pt*72.27/96

New_classification <- data.frame(type = c("coding UCR", "ncRNA UCR", "other UCR"), 
                                 num = c(314, 66, 99))
New_classification$type <- factor(x= New_classification$type, 
                                  levels = c("other UCR", 
                                             "ncRNA UCR", 
                                             "coding UCR"))

p1 <- ggplot(data = New_classification) +
  geom_bar(stat = "identity", 
           mapping = aes(x = "", y = num, fill = type), 
           alpha = 0.8, width = 1) +   # width >= 1去除中心的杂点，但是好像本来就没有杂点
  
  geom_text(aes(y = num/2 + c(0, cumsum(num)[-(length(num))]), 
                x = 1.1, 
                label = c("coding
  UCR(314)", "ncRNA
  UCR(66)", "other
  UCR(99)")), 
            size = 7, size.unit = "pt") +
  
  coord_polar("y", start = 0) +
  
  scale_color_npg() +
  
  labs(x = "", 
       y = "", 
       title = "") +
  
  theme(
    panel.background = element_blank(), 
    panel.grid = element_blank(), 
    text = element_text(size = 7), 
    line = element_line(linewidth = 0.5/lwd_pt), 
    legend.position = "none", 
    
    aspect.ratio = 1:1, 
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    
    axis.ticks = element_blank(),  # 去除左上角的点
    axis.text.x = element_blank()   # 去掉白框的数字
  )
p1

ggview(p1, width = 4.8, height = 4.8, units = "cm", dpi = 1200)

pdf(file = "03-results/16-New_classification/newClassificationPieChart.pdf", width = 480/254, height = 480/254)
p1
dev.off()

# 获得03-intergenic -------------------------------------------------------

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(dbplyr)
library(tidyverse)

coding_UCR <- read.table(file = "02-analysis/16-New_classification/01-UCR_from_protein_coding_gene.bed", 
                         sep = "\t")
ncRNA_UCR <- read.table(file = "02-analysis/16-New_classification/02-UCR_from_Non_Coding_RNA.bed", 
                        sep = "\t")
temp_UCR <- rbind(coding_UCR[, c(1:4)], ncRNA_UCR[, c(1:4)])
temp_UCR <- unique(temp_UCR)
colnames(temp_UCR) <- colnames(UCR)[1:4]

UCR <- read.table(file = "01-data/UCR_raw/UCR_location.txt", sep = "\t", header = TRUE)
UCR_bed <- UCR[, c(1:4)]

intergenic_UCR <- anti_join(UCR_bed, temp_UCR, by = "UCR_name")
write.table(intergenic_UCR, file = "02-analysis/16-New_classification/03-UCR_intergenic.bed", 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

