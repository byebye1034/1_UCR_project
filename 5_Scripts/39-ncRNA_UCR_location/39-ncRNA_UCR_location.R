# 查看ncRNA UCR在lncRNA上的位置，是exon还是intron

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library(ggprism)

# 66个ncRNA UCR，65个ncRNA UCR和56个lncRNA重叠，uc.15和miRNA ENSG00000224592重叠
ncRNA_UCR <- read.table(file = "02-analysis/16-New_classification/02-UCR_from_Non_Coding_RNA(66).bed")
length(unique(ncRNA_UCR$V8))

# pie plot
ncUCRGene <- tibble(
  Location = c("lncRNA", "miRNA"), 
  Number = c(65, 1)
)

ncUCRGene$Location <- factor(x = ncUCRGene$Location, 
                             levels = c("lncRNA", "miRNA"), 
                             ordered = TRUE)

lwd_pt <- .pt*72.27/96

PiePlot <- function(data){
  ggplot(data = data) +
    
    geom_bar(mapping = aes(x = "", y = Number, fill = Location), 
             stat = "identity", alpha = 0.8, width = 1) +
    
    coord_polar("y", start = 0, direction = -1) +
    
    scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF")) +
    
    theme_prism(
      base_size = 7, 
      base_line_size = 0.5/lwd_pt, 
      base_rect_size = 0.5/lwd_pt, 
      border = FALSE
    ) +
    
    theme(
      panel.background = element_blank(), 
      panel.grid = element_blank(), 
      legend.position = "top", 
      legend.key.size = unit(0.25, "cm"), 
      
      aspect.ratio = 1:1, 
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      
      axis.ticks = element_blank(),  # 去除左上角的点
      axis.text.x = element_blank()   # 去掉白框的数字
    )
}

PiePlot(data = ncUCRGene)

pdf(file = "03-results/39-ncRNA_UCR_location/lncRNA_proportion.pdf", 
    width = 360/254, height = 360/254)
PiePlot(data = ncUCRGene)
dev.off()






























