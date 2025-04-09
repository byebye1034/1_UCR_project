# protein-protein/protein-link community/protein complex
# condensates

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(dbplyr)
library(tidyverse)
library(biomaRt) # 在使用biomaRt的时候dbplyr的版本不能太高
library(curl)

rf_neighboring_protein_coding <- read.table(file = "02-analysis/13-Network_analysis/rf_neighboring_protein_coding.bed")
rf_neighboring_protein_coding$V9[rf_neighboring_protein_coding$V9 == "."] <- 
  rf_neighboring_protein_coding$V8[rf_neighboring_protein_coding$V9 == "."]
rf_neighboring_protein_coding <- rf_neighboring_protein_coding[, c(4, 8, 9)]
colnames(rf_neighboring_protein_coding) <- c("rf_name", "ensembl_id", "neighboring_gene")

# obtain uniprot id，即RF相关的蛋白编码基因的ENTRY ID-------------------------------------------------------

ensembl_id <- rf_neighboring_protein_coding$ensembl_id
ensembl_id <- as.data.frame(ensembl_id)

# 利用useMart链接到人类的数据库，其他数据库使用listDatasets(ensembl)查看
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# attr为基因id目标类型，listAttributes(ensembl)展示其他id输出类型
attr <- c("ensembl_gene_id", "uniprot_gn_symbol", "uniprot_gn_id")

# listFilters(ensembl)展示其他id输入类型
ids <- getBM(attributes = attr, 
             filters = "ensembl_gene_id", 
             values = ensembl_id, 
             mart = ensembl)

ids <- unique(ids)
colnames(ids)[1] <- "ensembl_id"

# 将rf_neighboring_protein_coding和ids合并
rf_neighboring_protein_coding <- merge(rf_neighboring_protein_coding, ids, by = "ensembl_id")

# 以ensembl id为唯一标识，删除uniprot_gn_symbol这一列，基因名以neighboring_gene为准
# 去除没有uniprot_gn_id的行
rf_neighboring_protein_coding <- rf_neighboring_protein_coding[, c(-4)]
rf_neighboring_protein_coding <- filter(rf_neighboring_protein_coding, rf_neighboring_protein_coding$uniprot_gn_id != "")

save(rf_neighboring_protein_coding, file = "02-analysis/13-Network_analysis/rf_protein_uniprot_id.Rdata")

# 人类所有蛋白的reviewed之后的ENTRY ID信息
uniprot_reviewed_protein <- read.table(file = "01-data/13-Network_analysis/simple_uniprotkb_reviewed_true_AND_model_organ_2024_03_19.tsv", 
                                       sep = "\t", header = T)

rf_neighboring_protein_coding_reviewed <- 
  rf_neighboring_protein_coding[rf_neighboring_protein_coding$uniprot_gn_id %in% uniprot_reviewed_protein$Entry, ]
save(rf_neighboring_protein_coding_reviewed, file = "02-analysis/13-Network_analysis/rf_neighboring_protein_coding_reviewed.Rdata")

# 154个rf，rf.39位置上有两个蛋白编码基因，HBE1和OR51B5

# UCR相关的蛋白编码基因的ENRTY ID ---------------------------------------------------

library(biomaRt)

load("D:/R_project/UCR_project/02-analysis/07-UCR_Classification/UCR classification/UCR_noninter_nonintronic_classification.Rdata")
intergenic_neighboring_gene <- read.table(file = "02-analysis/12-karyoplote/intergenic_UCR_neighboring_gene.bed", 
                                          sep = "\t", header = FALSE)
intronic_neighboring_gene <- read.table(file = "02-analysis/12-karyoplote/intronic_UCR_neighboring_gene.bed", 
                                        sep = "\t", header = FALSE)

typeII_neighboring_gene <- UCR_classification[, c(1:2, 6)]
colnames(typeII_neighboring_gene) <- c("UCR_name", "UCR_type", "ensembl_id")
typeII_neighboring_gene <- typeII_neighboring_gene[, c(1, 3, 2)]

intergenic_neighboring_gene <- intergenic_neighboring_gene[, c(4, 8)]
intergenic_neighboring_gene[, 3] <- "intergenic"
colnames(intergenic_neighboring_gene) <- c("UCR_name", "ensembl_id", "UCR_type")

intronic_neighboring_gene <- intronic_neighboring_gene[, c(4, 8)]
intronic_neighboring_gene[, 3] <- "intronic"
colnames(intronic_neighboring_gene) <- c("UCR_name", "ensembl_id", "UCR_type")

all_UCR_neighboring_gene <- rbind(intergenic_neighboring_gene, 
                                  intronic_neighboring_gene, 
                                  typeII_neighboring_gene)
all_UCR_neighboring_gene <- all_UCR_neighboring_gene[!(all_UCR_neighboring_gene$UCR_name == "uc.18"), ]
all_UCR_neighboring_gene <- all_UCR_neighboring_gene[!(all_UCR_neighboring_gene$UCR_name == "uc.304"), ]

# 利用useMart链接到人类的数据库，其他数据库使用listDatasets(ensembl)查看
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# attr为基因id目标类型，listAttributes(ensembl)展示其他id输出类型
attr <- c("ensembl_gene_id", "uniprot_gn_symbol", "uniprot_gn_id")

# listFilters(ensembl)展示其他id输入类型
ids <- getBM(attributes = attr, 
             filters = "ensembl_gene_id", 
             values = all_UCR_neighboring_gene$ensembl_id, 
             mart = ensembl)
colnames(ids)[1] <- "ensembl_id"

all_UCR_neighboring_gene <- merge(all_UCR_neighboring_gene, ids, by = "ensembl_id")
all_UCR_neighboring_gene <- all_UCR_neighboring_gene[, -4]

all_UCR_neighboring_gene <- all_UCR_neighboring_gene[all_UCR_neighboring_gene$uniprot_gn_id %in% uniprot_reviewed_protein$Entry, ]
save(all_UCR_neighboring_gene, file = "02-analysis/13-Network_analysis/all_UCR_neighboring_protein_coding_reviewed.Rdata")

all_UCR_neighboring_gene_reviewed <- all_UCR_neighboring_gene

# protein-protein interaction ---------------------------------------------

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(dbplyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(gginnards)

load("D:/R_project/UCR_project/02-analysis/13-Network_analysis/all_UCR_neighboring_protein_coding_reviewed.Rdata")
load("D:/R_project/UCR_project/02-analysis/13-Network_analysis/rf_neighboring_protein_coding_reviewed.Rdata")

global_PPI_network <- read.table(file = "01-data/13-Network_analysis/Global Protein-Protein Interaction Network.txt", 
                                 sep = "\t", 
                                 header = TRUE)
global_PPI_network$interaction <- rowSums(global_PPI_network[, 3:ncol(global_PPI_network)])
global_PPI_network <- global_PPI_network[, c(1:2, 76)]

# 将不为0的值赋值为1
global_PPI_network$interaction[global_PPI_network$interaction != 0] <- 1

all_UCR_reviewed_protein_interaction <- merge(all_UCR_neighboring_gene_reviewed, 
                                              global_PPI_network, 
                                              by.x = "uniprot_gn_id", 
                                              by.y = "InteractorA")
all_UCR_reviewed_protein_interaction <- filter(all_UCR_reviewed_protein_interaction, 
                                               interaction > 0)

# 分开都和rf相比，这样可以获得UCR和rf的差异，还能获得UCR内部的一致性或者差异性

# UCR ---------------------------------------------------------------------
intergenic_reviewed_protein_interaction <- filter(all_UCR_reviewed_protein_interaction, 
                                                  UCR_type == "intergenic")
intronic_reviewed_protein_interaction <- filter(all_UCR_reviewed_protein_interaction, 
                                                UCR_type == "intronic")
typeII_reviewed_protein_interaction <- filter(all_UCR_reviewed_protein_interaction, 
                                              UCR_type != "intergenic" & UCR_type != "intronic")

# 使用group_by()函数按照第一列进行分组，然后使用summarize()函数对每个组中的第六列进行求和
intergenic_reviewed_protein_interaction_num <- intergenic_reviewed_protein_interaction %>% 
  group_by(uniprot_gn_id) %>% 
  summarise(interaction_num = sum(interaction)) %>% 
  mutate(UCR_type = "intergenic")

intronic_reviewed_protein_interaction_num <- intronic_reviewed_protein_interaction %>% 
  group_by(uniprot_gn_id) %>% 
  summarise(interaction_num = sum(interaction)) %>% 
  mutate(UCR_type = "intronic")

typeII_reviewed_protein_interaction_num <- typeII_reviewed_protein_interaction %>% 
  group_by(uniprot_gn_id) %>% 
  summarise(interaction_num = sum(interaction)) %>% 
  mutate(UCR_type = "typeII")

# rf ----------------------------------------------------------------------

load("D:/R_project/UCR_project/02-analysis/13-Network_analysis/rf_neighboring_protein_coding_reviewed.Rdata")

all_rf_reviewed_protein_interaction <- merge(rf_neighboring_protein_coding_reviewed, 
                                             global_PPI_network, 
                                             by.x = "uniprot_gn_id", 
                                             by.y = "InteractorA")
all_rf_reviewed_protein_interaction <- filter(all_rf_reviewed_protein_interaction, 
                                              interaction > 0) 
all_rf_reviewed_protein_interaction_num <- all_rf_reviewed_protein_interaction %>% 
  group_by(uniprot_gn_id) %>% 
  summarise(interaction_num = sum(interaction)) %>% 
  mutate(UCR_type = "ranfrag")

# PPI绘图 -------------------------------------------------------------------

UCR_pro_interaction_num <- rbind(intergenic_reviewed_protein_interaction_num, 
                                 intronic_reviewed_protein_interaction_num, 
                                 typeII_reviewed_protein_interaction_num, 
                                 all_rf_reviewed_protein_interaction_num)
UCR_pro_interaction_num$UCR_type <- factor(UCR_pro_interaction_num$UCR_type, 
                                           levels = c("ranfrag", 
                                                      "intergenic",
                                                      "intronic", 
                                                      "typeII"))

lwd_pt <- .pt*72.27/96
compare_list <- list(c("ranfrag", "intergenic"), 
                     c("ranfrag", "intronic"), 
                     c("ranfrag", "typeII"))

p1 <- ggplot(data = UCR_pro_interaction_num, mapping = aes(x = UCR_type, y = interaction_num)) +
  geom_boxplot(mapping = aes(color = UCR_type), 
               outlier.shape = NA, 
               linewidth = 0.5/lwd_pt)  +
  
  scale_y_continuous(limits = c(0, 90), expand = c(0, 0)) +
  
  scale_color_npg() +
  
  theme(
    panel.grid = element_blank(), 
    panel.background = element_blank(), 
    text = element_text(size = 7), 
    line = element_line(linewidth = 0.5/lwd_pt), 
    legend.position = "none",
    
    axis.text = element_text(size = 7), 
    axis.line = element_line(linewidth = 0.5/lwd_pt), 
    
    aspect.ratio = 1
  )
p1

pdf(file = "03-results/13-Network_analysis/PPI_num.pdf", width = 480/254, height = 480/254)
p1
dev.off()

# 用下面这张图拿到显著性标记

p2 <- ggplot(data = UCR_pro_interaction_num, mapping = aes(x = UCR_type, y = interaction_num)) +
  geom_boxplot(mapping = aes(color = UCR_type), 
               outlier.shape = NA, 
               linewidth = 0.5/lwd_pt)  +
  
  scale_y_continuous(limits = c(0, 2000)) +
  
  scale_color_npg() +
  
  theme(
    panel.grid = element_blank(), 
    panel.background = element_blank(), 
    text = element_text(size = 7), 
    line = element_line(linewidth = 0.5/lwd_pt), 
    legend.position = "none",
    
    axis.text = element_text(size = 7), 
    axis.line = element_line(linewidth = 0.5/lwd_pt), 
    
    aspect.ratio = 1
  )
p2

p3 <- p2 + stat_compare_means(
  comparisons = compare_list, 
  method = "wilcox.test", 
  label = "p.signif", 
  bracket.size = 0.5/lwd_pt, 
  step.increase = 0.25
)

p3$layers[[which_layers(p3, "GeomSignif")]]$aes_params$textsize <- 5/lwd_pt
p3

pdf(file = "03-results/13-Network_analysis/PPI_num_signif.pdf", width = 480/254, height = 480/254)
p3
dev.off()


# 