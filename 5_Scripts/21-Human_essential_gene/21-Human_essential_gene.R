# 整理所有的human essential gene

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)

# 获得human deg -------------------------------------------------------------

deg_eukaryotes <- read.csv(file = "01-data/21-Human_essential_gene/deg_eukaryotes.csv", 
                           sep = ";", header = FALSE)
human_deg_list <- deg_eukaryotes[c(grep("Homo sapiens", deg_eukaryotes$V1)), ]

human_deg <- read.csv(file = "01-data/21-Human_essential_gene/deg_annotation_e.csv", 
                      sep = ";", header = FALSE)
human_deg <- human_deg[c(human_deg$V1 %in% human_deg_list$V13), ]
human_deg <- human_deg[, c(3, 7)]
human_deg <- unique(human_deg)
human_deg <- human_deg %>% 
  group_by(V3) %>% 
  summarise(gene_description = paste(V7, collapse = "/"))
colnames(human_deg)[1] <- "gene_symbol"

human_deg$gene_symbol[1] <- "MARCHF5"
human_deg$gene_symbol[2] <- "SEPTIN5"
human_deg$gene_symbol[3] <- "MARCHF6"

# 统计UCR相关基因是否在human deg当中富集 -----------------------------------------------

library(biomaRt)
library(curl)

# 利用useMart链接到人类的数据库，其他数据库使用listDatasets(ensembl)查看
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# attr为基因id目标类型，listAttributes(ensembl)展示其他id输出类型
attr <- c("ensembl_gene_id", "hgnc_symbol")

# listFilters(ensembl)展示其他id输入类型
ids <- getBM(attributes = attr, 
             filters = "hgnc_symbol", 
             values = human_deg$gene_symbol, 
             mart = ensembl)
human_deg <- merge(human_deg, ids, by.x = "gene_symbol", by.y = "hgnc_symbol")
save(human_deg, file = "01-data/21-Human_essential_gene/human_deg.Rdata")

# UCR related gene
coding_UCR_genes <- read.table(file = "02-analysis/16-New_classification/coding_UCR_ensembl_id.txt", 
                               sep = "\t")
ncRNA_UCR_genes <- read.table(file = "02-analysis/16-New_classification/ncRNA_UCR_ensembl_id.txt", 
                              sep = "\t")
UCR_genes <- rbind(coding_UCR_gene, ncRNA_UCR_genes)
UCR_genes_in_EG <- merge(UCR_genes, human_deg, by.x = "V1", by.y = "ensembl_gene_id")

# rf related gene
coding_rf_genes <- read.table(file = "02-analysis/16-New_classification/01-rf_form_protein_coding_gene.bed", 
                              sep = "\t", header = FALSE)
coding_rf_genes <- coding_rf_genes[, c(8, 10)]

ncRNA_rf_genes <- read.table(file = "02-analysis/16-New_classification/02-rf_from_Non_Coding_RNA.bed", 
                             sep = "\t", header = FALSE)
ncRNA_rf_genes <- ncRNA_rf_genes[, c(8, 10)]
rf_genes <- rbind(coding_rf_genes, ncRNA_rf_genes)
rf_genes_in_EG <- merge(rf_genes, human_deg, by.x = "V8", by.y = "ensembl_gene_id")

# 统计分析 --------------------------------------------------------------------

# UCR genes:445
# UCR genes in EG:272
# rf genes:319
# rf genes in EG:138

# 现在我有两组基因，组名分别是UCR和rf，其中分别有445和319个基因，
# 在UCR组中有272个基因是EG，其余的是NEG，在rf组中有138个是EG，其余的是NEG，
# 我想用统计学方法检验UCR组中的基因是不是更多的会位于EG当中。我可以使用哪些检验方法？

# 创建一个2x2的列联表
gene_data <- matrix(c(272, 173, 138, 181), nrow = 2, byrow = TRUE, 
                    dimnames = list(c("UCR", "rf"), c("EG", "NEG")))

# 执行卡方检验
chisq_result <- chisq.test(gene_data)

# 输出检验结果
print(chisq_result)

# 绘制统计图 -------------------------------------------------------------------

library(reshape2)
library(ggplot2)
library(ggsci)
library(scales) # 将Y轴标签改成百分比
library(gginnards)
library(ggprism)
library(sysfonts)
library(showtext)

font_add(family = "Arial", regular = "C:/Windows/Fonts/arial.ttf")
font_families()

lwd_pt <- .pt*72.27/96

# 创建数据框

gene_data <- data.frame(
  category = c("EG", "NEG"), 
  UCR = c(61.1236, 38.8764), 
  rf = c(43.26019, 56.73981)
)

# 将数据从宽格式转换为长格式
gene_data_L <- melt(gene_data, id.vars = "category", variable.name = "type", value.name = "proportion")
gene_data_L$category <- factor(gene_data_L$category, 
                               levels = c("NEG", "EG"), 
                               ordered = TRUE)

# 绘制堆叠柱形图

p1 <- ggplot(data = gene_data_L) +
  geom_bar(mapping = aes(x = type, y = proportion, fill = category), 
           stat = "identity", position = "fill", width = 0.6) + 
  
  scale_y_continuous(labels = label_percent(), expand = c(0, 0)) +
  
  scale_fill_manual(values = c("EG" = "#00A087FF", "NEG" = "#3C5488FF"), 
                    name = "") +
  
  scale_x_discrete(labels = c("UCR" = "UCR", "rf" = "RF")) +
  
  geom_text(aes(x = type, y = proportion/2,   # 在内部添加具体数值
                label = paste0(round(proportion, 1), "%")), 
            position = position_fill(vjust = 0.5), 
            size = 7, size.unit = "pt") +
  
  labs(x = "", 
       y = "Percent (%)", 
       title = "p-value = 1.513e-06") +
  
  theme_prism(base_size = 7, 
              base_family = "Arial", 
              base_fontface = "plain", 
              base_line_size = 0.5/lwd_pt, 
              base_rect_size = 0.5/lwd_pt, 
              border = TRUE) +
  
  theme(
    legend.position = "top", 
    legend.key.size = unit(0.25, "cm")
  )

p1

plotHumanEssentialGene <- p1

ggview::ggview(plot = ggplot2::last_plot(), 
               width = 4.8, height = 6, units = "cm", dpi = 1200)

ggsave(filename = "03-results/21-Human_essential_gene/proportion_UCR_gene_in_EG.tiff", 
       width = 4.8, height = 6, units = "cm", dpi = 1200)

pdf(file = "03-results/21-Human_essential_gene/proportion_UCR_gene_in_EG.pdf", 
    width = 480/254, height = 600/254)

# 使用showtext渲染字体
showtext_begin()

p1

# 关闭showtext渲染字体
showtext_end()

dev.off()




