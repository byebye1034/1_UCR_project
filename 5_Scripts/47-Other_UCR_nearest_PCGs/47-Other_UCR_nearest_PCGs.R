# get other UCRs nearest upstream and downstream protein coding genes

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library(readxl)
library(stringr)
library(ggplot2)
library(ggview)
library(cowplot)

otherUCRsNearestPCGsINFO <- read.table(
  file = "01-data/47-Other_UCR_nearest_PCGs/UCR_closest_genes.bed", 
  fill = TRUE
)
otherUCRsNearestPCGsINFO <- subset(otherUCRsNearestPCGsINFO, V6 != "NA")
otherUCRsNearestPCGsINFO <- otherUCRsNearestPCGsINFO[, c(-11:-14)]
otherUCRsNearestPCGsINFO <- otherUCRsNearestPCGsINFO %>% 
  mutate(nearestGene = ifelse(abs(V10) > V20, V19, V9))

# 47 other UCRs maybe functions as enhancers
otherUCRsOverlapcCREs <- read.table(
  file = "02-analysis/36-Encode_screen_cCREs/cCREs_overlap_intergenic.bed"
)
length(table(otherUCRsOverlapcCREs$V4)) # 47

nearestGeneOtherUCRsOverlapcCREs <- otherUCRsNearestPCGsINFO %>% 
  filter(V4 %in% otherUCRsOverlapcCREs$V4)
write.csv(nearestGeneOtherUCRsOverlapcCREs, 
          file = "03-results/47-Other_UCR_nearest_PCGs/nearestPCGsOtherUCRsOverlapcCREs.csv", 
          row.names = FALSE)

nearestGeneOtherUCRsNotOverlapcCREs <- anti_join(otherUCRsNearestPCGsINFO, 
                                                 nearestGeneOtherUCRsOverlapcCREs, 
                                                 by= "V4")
write.csv(nearestGeneOtherUCRsNotOverlapcCREs, 
          file = "03-results/47-Other_UCR_nearest_PCGs/nearestPCGsOtherUCRsNotOverlapcCREs.csv", 
          row.names = FALSE)

# goOfOtherUCRsAdjacentPCGs -----------------------------------------------

## cCREs Overlap
goOfAdjacentPCGsOfOtherUCRsOfcCREs <- read_xlsx(
  path = "03-results/47-Other_UCR_nearest_PCGs/nearestPCGsOtherUCRsOverlapcCREsGO/metascape_result.xlsx", 
  sheet = 2
)
goOfAdjacentPCGsOfOtherUCRsOfcCREs <- 
  goOfAdjacentPCGsOfOtherUCRsOfcCREs[str_detect(goOfAdjacentPCGsOfOtherUCRsOfcCREs$GroupID, ".*Summary"), ]

# order
goOfAdjacentPCGsOfOtherUCRsOfcCREs <- 
  goOfAdjacentPCGsOfOtherUCRsOfcCREs[order(goOfAdjacentPCGsOfOtherUCRsOfcCREs$LogP, decreasing = TRUE), ]
goOfAdjacentPCGsOfOtherUCRsOfcCREs$Description <- 
  factor(x = goOfAdjacentPCGsOfOtherUCRsOfcCREs$Description, 
         levels = goOfAdjacentPCGsOfOtherUCRsOfcCREs$Description, 
         ordered = TRUE)

lwd_pt <- .pt*72.27/96

plotGOOtherUCRsOverlapcCREs <- ggplot(data = goOfAdjacentPCGsOfOtherUCRsOfcCREs) +
  
  geom_bar(mapping = aes(x = -(LogP), y = Description), 
           stat = "identity", width = 0.8, fill = "#E3738B", alpha = 0.5) +
  
  geom_text(aes(x = 0.1, y = Description, label = Description), 
            size = 6, size.unit = "pt", hjust = 0) +
  
  scale_x_continuous(limits = c(0, 10), expand = c(0, 0)) +
  
  labs(x = "-Log10 pvalue", 
       y = "GO biological process 
       (adjacent PCGs of other UCRs overlapping cCREs)", 
       title = "") +
  
  theme(
    plot.title = element_text(size = 6),
    legend.key.height = unit(100/254, "cm"), # 设置图例的高度为2厘米
    legend.key.width = unit(100/254, "cm"),  # 设置图例的宽度为1厘米
    
    panel.grid = element_blank(), 
    panel.background = element_blank(), 
    text = element_text(size = 6), 
    line = element_line(linewidth = 0.5/lwd_pt), 
    
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(), 
    axis.line.y = element_line(linewidth = 0.5/lwd_pt, color = "black"), 
    
    axis.text.x = element_text(size = 6, color = "black"), 
    axis.line = element_line(linewidth = 0.5/lwd_pt, color = "black")
  )

plotGOOtherUCRsOverlapcCREs

ggview(plotGOOtherUCRsOverlapcCREs, width = 8.25, height = 6, units = "cm", dpi = 1200)

## No cCREs Overlap
goOfAdjacentPCGsOfOtherUCRsNotOverlapcCREs <- read_xlsx(
  path = "03-results/47-Other_UCR_nearest_PCGs/neatestPCGsOtherUCRsNotOverlapcCREs/metascape_result.xlsx", 
  sheet = 2
)
goOfAdjacentPCGsOfOtherUCRsNotOverlapcCREs <- 
  goOfAdjacentPCGsOfOtherUCRsNotOverlapcCREs[str_detect(goOfAdjacentPCGsOfOtherUCRsNotOverlapcCREs$GroupID, ".*Summary"), ]

# order
goOfAdjacentPCGsOfOtherUCRsNotOverlapcCREs <- 
  goOfAdjacentPCGsOfOtherUCRsNotOverlapcCREs[order(goOfAdjacentPCGsOfOtherUCRsNotOverlapcCREs$LogP, decreasing = TRUE), ]
goOfAdjacentPCGsOfOtherUCRsNotOverlapcCREs$Description <- 
  factor(
    x = goOfAdjacentPCGsOfOtherUCRsNotOverlapcCREs$Description, 
    levels = goOfAdjacentPCGsOfOtherUCRsNotOverlapcCREs$Description, 
    ordered = TRUE
  )

plotGOOtherUCRsNotOverlapcCREs <- ggplot(data = goOfAdjacentPCGsOfOtherUCRsNotOverlapcCREs) +
  
  geom_bar(mapping = aes(x = -(LogP), y = Description), 
           stat = "identity", width = 0.8, fill = "#F9D5DD", alpha = 0.5) +
  
  geom_text(aes(x = 0.1, y = Description, label = Description), 
            size = 6, size.unit = "pt", hjust = 0) +
  
  scale_x_continuous(limits = c(0, 8), expand = c(0, 0)) +
  
  labs(x = "-Log10 pvalue", 
       y = "GO biological process 
       (adjacent PCGs of other UCRs not overlapping cCREs)", 
       title = "") +
  
  theme(
    plot.title = element_text(size = 6),
    legend.key.height = unit(100/254, "cm"), # 设置图例的高度为2厘米
    legend.key.width = unit(100/254, "cm"),  # 设置图例的宽度为1厘米
    
    panel.grid = element_blank(), 
    panel.background = element_blank(), 
    text = element_text(size = 6), 
    line = element_line(linewidth = 0.5/lwd_pt), 
    
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(), 
    axis.line.y = element_line(linewidth = 0.5/lwd_pt, color = "black"), 
    
    axis.text.x = element_text(size = 6, color = "black"), 
    axis.line = element_line(linewidth = 0.5/lwd_pt, color = "black")
  )

plotGOOtherUCRsNotOverlapcCREs

cowplot::plot_grid(plotGOOtherUCRsOverlapcCREs, 
                   plotGOOtherUCRsNotOverlapcCREs, 
                   align = c("hv"), axis = c("tb"), 
                   nrow = 1, ncol = 2, rel_widths = 1, rel_heights = 1, 
                   labels = c("a", "b"), label_size = 8)

pdf(file = "03-results/47-Other_UCR_nearest_PCGs/barplotGOPCGsOfOtherUCRsOverlapcCREsOrNot.pdf", 
    width = 1650/254, height = 600/254)

dev.off()

library(readxl)

typeIIIUCRs <- read_xlsx(
  path = "D:/C_英文论文/Figures/Additional File2_Supplementary Tables.xlsx",
  sheet = 3, skip = 2
)

lwd_pt <- .pt*72.27/96

p1 <- ggplot(
  data = typeIIIUCRs, 
  mapping = aes(x = log10(distance...16))
) +
  
  geom_density(linewidth = 0.5/lwd_pt, color = "#E64B35FF") +
  geom_rug(linewidth = 0.5/lwd_pt, color = "#E64B35FF") +
  
  geom_vline(xintercept = log10(2000), 
             linetype = "dashed", 
             color = "gray", 
             linewidth = 0.5/lwd_pt) +
  
  scale_y_continuous(expand = c(0, 0)) +
  
  labs(
    x = "log10(distance between type III
    UCR and protein coding genes)", 
    y = "density", 
    title = "Human"
  ) +
  
  theme(
    panel.background = element_blank(), 
    panel.grid = element_blank(), 
    line = element_line(linewidth = 0.5/lwd_pt), 
    
    axis.text = element_text(size = 5, color = "#000000"), 
    axis.ticks = element_line(linewidth = 0.5/lwd_pt), 
    axis.line = element_line(linewidth = 0.5/lwd_pt), 
    axis.title = element_text(size = 5, color = "#000000"), 
    
    plot.title = element_text(size = 5),
    
    aspect.ratio = 1:1
  )
p1   # 距离主要在10000，一万bp以上

p1 + canvas(width = 4.8, height = 6, units = "cm", dpi = 1200)

## mouse type III UCRs nearest PCGs
mouseIIIUp <- read.table(
  file = "02-analysis/47-Other_UCR_nearest_PCGs/mouseTypeIIIUCRsClosestUpstreamPCGs.bed", 
  sep = "\t"
)
mouseIIIUp <- unique(mouseIIIUp)
mouseIIIUp <- mouseIIIUp[, c(4, 9, 10)]

mouseIIIDown <- read.table(
  file = "02-analysis/47-Other_UCR_nearest_PCGs/mouseTypeIIIUCRsClosestDownstreamPCGs.bed", 
  sep = "\t"
)
mouseIIIDown <- unique(mouseIIIDown)
mouseIIIDown <- mouseIIIDown[, c(4, 9, 10)]

mouseIIIGenes <- merge(mouseIIIUp, mouseIIIDown, 
                       by.x = "V4", by.y = "V4")
mouseIIIGenes$Genes <- if_else(
  mouseIIIGenes$V10.x + mouseIIIGenes$V10.y < 0, mouseIIIGenes$V9.y, mouseIIIGenes$V9.x
)

write.csv(mouseIIIGenes, 
          file = "03-results/47-Other_UCR_nearest_PCGs/mouseTypeIIIUCRsNearestGenes.csv", 
          row.names = FALSE)
