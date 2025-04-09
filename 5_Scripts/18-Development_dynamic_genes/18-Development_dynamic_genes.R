# 看在不同器官里DDG的数目

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(readxl)
library(ggplot2)

development_dynamic_organ <- read_xlsx(path = "01-data/18-Development_dynamic_genes/development_dynamic_organ.xlsx", 
                                        sheet = 7, skip = 2)
development_dynamic_organ <- development_dynamic_organ[, c(1, 10:16)]

# coding UCR --------------------------------------------------------------

coding_UCR <- read.table(file = "02-analysis/16-New_classification/01-UCR_from_protein_coding_gene.bed", 
                         sep = "\t")
coding_UCR <- coding_UCR[, c(8, 9)]
coding_UCR <- unique(coding_UCR)
coding_UCR <- merge(coding_UCR, development_dynamic_organ, by.x = "V8", by.y = "Human_ID")

# 假设数据框名为coding_UCR，将第4到第10列转换为数值型，然后计算它们的和
cols_to_sum <- coding_UCR[, 3:9]

# 将数据转换为数值型（如果有非数值型数据，则会被转换为NA）
cols_to_sum <- apply(cols_to_sum, 2, as.numeric)

# 计算和并将结果存储在新列中
coding_UCR$OrganDDGNum <- rowSums(cols_to_sum, na.rm = TRUE)
save(coding_UCR, file = "02-analysis/18-Development_dynamic_genes/coding_UCR_OrganDDGNum.Rdata")

print(table(coding_UCR$OrganDDGNum))
# 0  1  2  3  4  5  6  7 
# 6 14 13 27 36 43 36 24

# random fragments --------------------------------------------------------

coding_rf <- read.table(file = "02-analysis/16-New_classification/01-rf_form_protein_coding_gene.bed", sep = "\t")
coding_rf <- coding_rf[, c(8, 9)]
coding_rf <- unique(coding_rf)

coding_rf <- merge(coding_rf, development_dynamic_organ, by.x = "V8", by.y = "Human_ID")

# 假设数据框名为coding_UCR，将第4到第10列转换为数值型，然后计算它们的和
cols_to_sum <- coding_rf[, 3:9]

# 将数据转换为数值型（如果有非数值型数据，则会被转换为NA）
cols_to_sum <- apply(cols_to_sum, 2, as.numeric)

# 计算和并将结果存储在新列中
coding_rf$OrganDDGNum <- rowSums(cols_to_sum, na.rm = TRUE)
save(coding_rf, file = "02-analysis/18-Development_dynamic_genes/coding_rf_OrganDDGNum.Rdata")

print(table(coding_rf$OrganDDGNum))
# 0  1  2  3  4  5  6  7 
# 13  9 23 36 48 66 31 15

# lollipop 棒棒糖图 -----------------------------------------------------------

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(reshape2)
library(ggsci)

development_dynamic_genes <- data.frame(num = c("DNDG", "1", "2", "3", "4", "5", "6", "7"), 
                                        coding_UCR = c("6", "14", "13", "27", "36", "43", "36", "24"), 
                                        coding_rf = c("13", "9", "23", "36", "48", "66", "31", "15"))
development_dynamic_genes[, 2] <- as.numeric(development_dynamic_genes[, 2])
development_dynamic_genes[, 3] <- as.numeric(development_dynamic_genes[, 3])

development_dynamic_genes$coding_UCR <- ((development_dynamic_genes$coding_UCR)/199)*100
development_dynamic_genes$coding_rf <- ((development_dynamic_genes$coding_rf)/241)*100

development_dynamic_genes_L <- melt(development_dynamic_genes, id.vars = "num", variable.name = "DDG_proportion", value.name = "Num")

lwd_pt <- .pt*72.27/96

development_dynamic_genes_L$num <- factor(x = development_dynamic_genes_L$num, 
                                          levels = c("DNDG", "1", "2", "3", "4", "5", "6", "7"))

# 点图 + 线图
p1 <- ggplot() +
  
  geom_point(data = development_dynamic_genes_L[c(1:8), ], 
             mapping = aes(x = num, y = Num), 
             position = position_nudge(x = 0.1, y = 0), 
             color = "#E64B35FF", size = 0.5) +
  
  geom_point(data = development_dynamic_genes_L[c(9:16), ], 
             mapping = aes(x = num, y = Num), 
             position = position_nudge(x = -0.1, y = 0), 
             color = "#4DBBD5FF", size = 0.5) +
  
  geom_segment(data = development_dynamic_genes_L[c(1:8), ], 
               mapping = aes(x = num, xend = num, y = 0, yend = Num), 
               linewidth = 0.5/lwd_pt, color = "#E64B35FF", 
               position = position_nudge(x = 0.1, y = 0)) +
  
  geom_segment(data = development_dynamic_genes_L[c(9:16), ], 
               mapping = aes(x = num, xend = num, y = 0, yend = Num), 
               linewidth = 0.5/lwd_pt, color = "#4DBBD5FF", 
               position = position_nudge(x = -0.1, y = 0)) +
  
  scale_y_continuous(limits = c(0, 30), expand = c(0, 0)) +
  
  theme(
    panel.background = element_blank(), 
    panel.grid = element_blank(), 
    line = element_line(linewidth = 0.5/lwd_pt), 
    text = element_text(size = 7, color = "black"), 
    legend.position = "NONE", 
    
    axis.line = element_line(linewidth = 0.5/lwd_pt), 
    axis.text = element_text(size = 7, color = "black"), 
    
    aspect.ratio = 1:1
  ) +
  
  labs(y = "Proportion of developmentally
            dynamic genes(DDGs)(%)", 
       x = "No. of organs")
p1

pdf(file = "03-results/18-Development_dynamic_genes/development_dynamic_genes_proportion.pdf", width = 480/254, height = 480/254)
p1
dev.off()

# 进行Fisher's exact test
data <- development_dynamic_genes[, c(2:3)]
data <- round(data)
result <- fisher.test(data)
print(result)

# mcnemar.test
coding_UCR <- development_dynamic_genes$coding_UCR
coding_rf <- development_dynamic_genes$coding_rf
result <- mcnemar.test(coding_UCR, coding_rf)

# wilcox.test
wilcox.test(coding_UCR, coding_rf, paired = TRUE)



