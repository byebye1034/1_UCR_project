# 利用ggtree绘制进化树

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(ggtree)
library(ggplot2)
library(ggsci)
library(patchwork)

# 进化树 ---------------------------------------------------------------------

tree <- read.tree(file = "01-data/26-Evolution_newly_emerging_UCR_distribution/phyliptree_no_HRM.phy")
tree[["tip.label"]] <- gsub("'", "", tree[["tip.label"]])

tree.a <- full_join(tree, group, by.x = "", by.y = "species")

p1 <- ggtree(tr = tree, linewidth = 0.5/lwd_pt) +
  coord_flip()
p1

# 柱形图 ---------------------------------------------------------------------

UCR_num <- data.frame(species = c("Pan troglodytes", "Homo sapiens", "Macaca mulatta", 
                            "Bos taurus", "Ovis aries", "Sus scrofa", 
                            "Canis lupus familiaris", "Felis catus", "Mus musculus", 
                            "Rattus", "Equus caballus", "Ochotona princeps", 
                            "Danio rerio", "Takifugu rubripes", "Gallus gallus", 
                            "Meleagris gallopavo", "Neoceratodus forsteri", "Xenopus laevis", 
                            "Caenorhabditis elegans", "Drosophila melanogaster", "Saccharomyces cerevisiae"), 
                UCR_num = c(416, 479, 386, 
                            235, 234, 263, 
                            261, 281, 479, 
                            479, 257, 155, 
                            0, 0, 31, 
                            30, 0, 0, 
                            0, 0, 0), 
                )

UCR_num <- data.frame(species = c("Pantr", "Homo", "Macaca", 
                                  "Bosta", "Ovis", "Sussc", 
                                  "Canis", "Felis", "Equus", "Ochpri", 
                                  "Danio", "fugu", "Gallus", 
                                  "Mele", "Neofo", "Xenopus", 
                                  "elegans", "Dromel", "Saccer"), 
                      UCR_num = c(416, 479, 386, 
                                  235, 234, 263, 
                                  261, 281, 257, 155, 
                                  0, 0, 31, 
                                  30, 0, 0, 
                                  0, 0, 0), 
                      group = c(rep("Primates", 3),    # 灵长目
                                rep("Artiodactyla", 3),    # 偶蹄目
                                rep("Carnivora", 2),    # 食肉目
                                rep("Perissodactyla", 1),    # 奇蹄类
                                rep("Lagomorpha", 1),   # 兔形目
                                rep("Actinopteri", 2),   # 辐鳍鱼纲
                                rep("Phasianidae", 2),   # 稚科
                                rep("Ceratodontiformes", 1),   # 肉鳍鱼
                                rep("Amphibia", 1),   # 两栖类
                                rep("Eukaryota", 3)  # 真核生物
                      ))

UCR_num$species <- factor(x = UCR_num$species, 
                          levels = UCR_num$species)
UCR_num$species <- fct_rev(UCR_num$species)

lwd_pt <- .pt*72.27/96

npg_palette <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF",
                 "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF")
reversed_npg_plette <- rev(npg_palette)

p2 <- ggplot(data = UCR_num, aes(x = species, y = UCR_num, group = group, fill = group)) +
  geom_bar(stat = "identity", linewidth = 0.5/lwd_pt) +
  
  scale_y_continuous(limits = c(0, 500), expand = c(0, 0)) +
  
  scale_fill_manual(values = reversed_npg_plette) +  # 设置颜色映射
  
  labs(
    x = "", 
    y = "No. of UCR", 
    title = ""
  ) +
  
  theme(
    panel.background = element_blank(), 
    panel.grid = element_blank(), 
    text = element_text(size = 7, color = "black"), 
    line = element_line(linewidth = 0.5/lwd_pt), 
    legend.position = "none", 
    
    axis.text.x = element_text(size = 7,
                               angle = 90, 
                               vjust = 0.5, 
                               color = "black"), 
    axis.text.y = element_text(size = 7, 
                               color = "black"), 
    axis.line = element_line(linewidth = 0.5/lwd_pt)
  )
p2

pdf(file = "03-results/26-Evolution_newly_emerging_UCR_distribution/evolution_UCR_num.pdf", width = 960/254, height = 720/254)
(p2 / p1) +
  plot_layout(nrow = 2, 
              heights = c(3, 0.75))
dev.off()


UCR_num <- data.frame(species = c("Pan troglodytes", "Homo sapiens", "Macaca mulatta", 
                                  "Bos taurus", "Ovis aries", "Sus scrofa", 
                                  "Canis lupus familiaris", "Felis catus", "Equus caballus", "Ochotona princeps", 
                                  "Danio rerio", "Takifugu rubripes", "Gallus gallus", 
                                  "Meleagris gallopavo", "Neoceratodus forsteri", "Xenopus laevis", 
                                  "Caenorhabditis elegans", "Drosophila melanogaster", "Saccharomyces cerevisiae"), 
                      UCR_num = c(416, 479, 386, 
                                  235, 234, 263, 
                                  261, 281, 257, 155, 
                                  0, 0, 31, 
                                  30, 0, 0, 
                                  0, 0, 0), 
                      group = c(rep("Primates", 3),    # 灵长目
                                rep("Artiodactyla", 3),    # 偶蹄目
                                rep("Carnivora", 2),    # 食肉目
                                rep("Perissodactyla", 1),    # 奇蹄类
                                rep("Lagomorpha", 1),   # 兔形目
                                rep("Actinopteri", 2),   # 辐鳍鱼纲
                                rep("Phasianidae", 2),   # 稚科
                                rep("Ceratodontiformes", 1),   # 肉鳍鱼
                                rep("Amphibia", 1),   # 两栖类
                                rep("Eukaryota", 3)  # 真核生物
                      ))

# human 479
# rat 479
# mouse 479
# chimpanzee 416
# rhemonkey 386
# cat 281
# pig 263
# dog 261
# horse 257
# cow 235
# sheep 234
# mochpri 155
# chicken 31
# turkey 30
# fugu 0
# zebrafish 0
# x.tropicalis 0
# lizard 0
# lungfish 0

# 20240528 绘制折线图 ----------------------------------------------------------

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(ggview)

load("D:/R_project/UCR_project/03-results/26-Evolution_newly_emerging_UCR_distribution/summary_1_before.Rdata")
load("D:/R_project/UCR_project/03-results/26-Evolution_newly_emerging_UCR_distribution/summary_2_added.Rdata")

evolution_summary <- rbind(
  summary_1, 
  summary_2
)

evolution_summary$species <- gsub("c.elegans", "Caenorhabditis elegans", evolution_summary$species)
evolution_summary$species <- gsub("d.melanogaster", "Drosophila melanogaster", evolution_summary$species)
evolution_summary$species <- gsub("qisaiman", "Petromyzon marinus", evolution_summary$species)
evolution_summary$species <- gsub("ribenhou", "Tachypleus tridentatus", evolution_summary$species)
evolution_summary$species <- gsub("wenchangyu", "Branchiostoma lanceolatum", evolution_summary$species)
evolution_summary$species <- gsub("daxiyangxueyu", "Gadus morhua", evolution_summary$species)
evolution_summary$species <- gsub("hongzun", "Oncorhynchus mykiss", evolution_summary$species)
evolution_summary$species <- gsub("fugu", "Takifugu rubripes", evolution_summary$species)
evolution_summary$species <- gsub("niluoluofeiyu", "Oreochromis niloticus", evolution_summary$species)
evolution_summary$species <- gsub("lungfish", "Neoceratodus forsteri", evolution_summary$species)
evolution_summary$species <- gsub("xiaoxinggouyu", "Scyliorhinus canicula", evolution_summary$species)
evolution_summary$species <- gsub("dabaisha", "Carcharodon carcharias", evolution_summary$species)
evolution_summary$species <- gsub("qingwa", "Rana temporaria", evolution_summary$species)
evolution_summary$species <- gsub("lizard", "Iguana", evolution_summary$species)
evolution_summary$species <- gsub("miandianmang", "Python bivittatus", evolution_summary$species)
evolution_summary$species <- gsub("turkey", "Meleagris gallopavo", evolution_summary$species)
evolution_summary$species <- gsub("banxiongcaoque", "Taeniopygia guttata", evolution_summary$species)
evolution_summary$species <- gsub("chicken", "Gallus gallus", evolution_summary$species)
evolution_summary$species <- gsub("jiage", "Columba livia", evolution_summary$species)
evolution_summary$species <- gsub("pingtadaoxianggui", "Chelonoidis abingdonii", evolution_summary$species)
evolution_summary$species <- gsub("dog", "Canis lupus familiaris", evolution_summary$species)
evolution_summary$species <- gsub("mochpri", "Ochotona princeps", evolution_summary$species)
evolution_summary$species <- gsub("sheep", "Ovis aries", evolution_summary$species)
evolution_summary$species <- gsub("cow", "Bos taurus", evolution_summary$species)
evolution_summary$species <- gsub("horse", "Equus caballus", evolution_summary$species)
evolution_summary$species <- gsub("pig", "Sus scrofa", evolution_summary$species)
evolution_summary$species <- gsub("cat", "Felis catus", evolution_summary$species)
evolution_summary$species <- gsub("haitun", "Tursiops truncatus", evolution_summary$species)
evolution_summary$species <- gsub("human", "Homo sapiens", evolution_summary$species)
evolution_summary$species <- gsub("chimpanzee", "Pan troglodytes", evolution_summary$species)
evolution_summary$species <- gsub("rat", "Rattus", evolution_summary$species)
evolution_summary$species <- gsub("mouse", "Mus musculus", evolution_summary$species)
evolution_summary$species <- gsub("x.tropicalis", "Xenopus laevis", evolution_summary$species)
evolution_summary$species <- gsub("zebrafish", "Danio rerio", evolution_summary$species)

# 删除酵母
evolution_summary <- evolution_summary[!(evolution_summary$species == "s.cerevisiae"), ]

lwd_pt <- .pt*72.27/96
windowsFonts(Arial = windowsFont('Arial'), TNR = windowsFont('Times New Roman'))

p1 <- ggplot(data = evolution_summary) +
  
  geom_line(mapping = aes(x = species, y = number, group = identity, color = identity)) +
  
  scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF")) +
  
  theme(
    panel.background = element_blank(), 
    panel.grid = element_blank(), 
    text = element_text(size = 7, color = "black"), 
    line = element_line(linewidth = 0.5/lwd_pt), 
    legend.position = "none", 
    
    axis.text.x = element_text(size = 7,
                               angle = 90, 
                               vjust = 0.5, 
                               color = "black"), 
    axis.text.y = element_text(size = 7, 
                               color = "black"), 
    axis.line = element_line(linewidth = 0.5/lwd_pt)
  )

p1

# 物种的进化树暂时画不出来，我觉得就先按照>95%的数目大小从小到大排列好了
evolution_summary <- evolution_summary %>% 
  arrange(desc(identity), number)

evolution_summary$species <- factor(x = evolution_summary$species[1:34], 
                  levels = evolution_summary$species[1:34])

p3 <- ggplot(data = evolution_summary) +
  
  geom_line(mapping = aes(x = species, y = number, group = identity, color = identity)) +
  
  scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF")) +
  scale_y_continuous(limits = c(0, 500), expand = c(0, 0)) +
  
  labs(x = "", 
       y = "Number of UCR") +
  
  theme(
    panel.background = element_blank(), 
    panel.grid = element_blank(), 
    text = element_text(size = 7, color = "black"), 
    line = element_line(linewidth = 0.5/lwd_pt), 
    legend.justification = c(0.05, 1), 
    legend.position = c(0.05, 1), 
    legend.key.size = unit(0.25, "cm"), 
    legend.text = element_text(size = 7), 
    legend.title = element_text(size = 7),
    
    axis.text.x = element_text(size = 7,
                               angle = 90, 
                               vjust = 0.5, 
                               color = "black"), 
    axis.text.y = element_text(size = 7, 
                               color = "black"), 
    axis.line = element_line(linewidth = 0.5/lwd_pt)
  )

p3

ggview(last_plot(), width = 9.6, height = 7, units = "cm", dpi = 1200)
ggsave(filename = "03-results/26-Evolution_newly_emerging_UCR_distribution/evolution_UCR_num_lineplot.tiff", 
       width = 9.6, height = 6, units = "cm", dpi = 1200)

pdf(file = "03-results/26-Evolution_newly_emerging_UCR_distribution/evolution_UCR_num_lineplot.pdf", 
    width = 960/254, height = 700/254)
p3
dev.off()

# Pan troglodytes
# Homo sapiens
# Macaca mulatta
# Bos taurus
# Ovis aries
# Sus scrofa
# Canis lupus familiaris
# Felis catus
# Mus musculus
# Rattus
# Equus caballus
# Ochotona princeps
# Danio rerio
# Takifugu rubripes
# Gallus gallus
# Meleagris gallopavo
# Neoceratodus forsteri
# Xenopus laevis
# Caenorhabditis elegans
# Drosophila melanogaster
# Saccharomyces cerevisiae
# Branchiostoma lanceolatum
# Petromyzonmarinus
# Scyliorhinus canicula
# Carcharodon carcharias
# Tachypleus tridentatus
# Oncorhynchus mykiss
# Oreochromis niloticus
# Gadus morhua
# Rana temporaria
# Iguana
# Python bivittatus
# Columba livia
# Melopsittacus undulatus
# Taeniopygia guttata
# Tursiops truncatus

tree <- read.tree(file = "01-data/26-Evolution_newly_emerging_UCR_distribution/phyliptree_3.phy")
tree[["tip.label"]] <- gsub("'", "", tree[["tip.label"]])

p2 <- ggtree(tr = tree, linewidth = 0.5/lwd_pt) +
  
  geom_text(mapping = aes(label = label), size = 7, size.unit = "pt", 
            hjust = -0.1) +
  
  theme(
    text = element_text(size = 7, color = "#000000")
  )
p2

# 从timetree网站得到的进化树
tree <- read.tree(file = "01-data/26-Evolution_newly_emerging_UCR_distribution/spelices_for_evo_tree.nwk")
tree[["tip.label"]] <- gsub("'", "", tree[["tip.label"]])

p3 <- ggtree(tr = tree, linewidth = 0.5/lwd_pt) +
  
  geom_text(mapping = aes(label = label), size = 7, size.unit = "pt", 
            hjust = -0.1) +
  
  theme(
    text = element_text(size = 7, color = "#000000")
  )
p3


old <- read.table(file = "C:/Users/baiyun6/Desktop/Taxonomy ID (1).txt", sep = "\t")
new <- read.table(file = "C:/Users/baiyun6/Desktop/Taxonomy ID.txt", sep = "\t")
new$V1 <- gsub(" - .*", "", new$V1)



