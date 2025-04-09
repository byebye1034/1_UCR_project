# annotate the distribution of type I UCRs: exonic & intronic?

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library(ggprism)
library(ggsci)
library(ggview)
library(ggalluvial)
library(readxl)
library(cowplot)

lwd_pt <- .pt*72.27/96

exonicUCRs <- read.table(
  file = "02-analysis/50-TypeIUCRs_PCG_Location/exonicUCRs.bed", 
  sep = "\t"
)
length(unique(exonicUCRs$V4)) # 126 UCRs
length(unique(exonicUCRs$V9)) # 124 genes

write.csv(unique(exonicUCRs[, c(8:9)]), 
          file = "03-results/50-TypeIUCRs_PCG_Location/exonicUCRsGenes.csv")

exonicUCRs <- exonicUCRs[, c(4, 9)] %>% 
  mutate(Type = "exon")
exonicUCRs <- unique(exonicUCRs)

# prime utr ---------------------------------------------------------------

primeUTRUCRs <- read.table(
  file = "02-analysis/50-TypeIUCRs_PCG_Location/primeutrUCRs.bed", sep = "\t"
)
primeUTRUCRs <- unique(primeUTRUCRs)

length(unique(primeUTRUCRs$V4)) # 91 UCRs
length(unique(primeUTRUCRs$V18)) # 81 genes

ud <- primeUTRUCRs[, c("V4", "V14")]
ud <- unique(ud)
table(ud$V14)
# five prime utr 22
# three prime utr 72


# intronic utr ------------------------------------------------------------

typeIUCRs <- read.table(
  file = "02-analysis/16-New_classification/01-UCR_from_protein_coding_gene.bed"
)
typeIUCRs <- unique(typeIUCRs[, c(4, 8)])
exonicUCRs <- unique(exonicUCRs[, c(4, 8)])
intronicUCRs <- setdiff(typeIUCRs, exonicUCRs)
write.csv(intronicUCRs, 
          file = "03-results/50-TypeIUCRs_PCG_Location/intronicUCRsGenes.csv")

# visualization -----------------------------------------------------------

typeIUCRsLocation <- tibble(
  Location = c("exon", "intron"), 
  Number = c(126, 188)
)

typeIUCRsLocation$Location <- factor(x = typeIUCRsLocation$Location, 
                                     levels = c("exon", "intron"), 
                                     ordered = TRUE)

PiePlot <- function(data){
  ggplot(data = data) +
    
    geom_bar(mapping = aes(x = "", y = Number, fill = Location), 
             stat = "identity", alpha = 0.8, width = 1) +
    
    coord_polar("y", start = 0, direction = -1) +
    
    scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF")) +
    
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

PiePlot(data = typeIUCRsLocation) + canvas(width = 4.8, height = 4.8, units = "cm", dpi = 1200)

pdf(file = "03-results/50-TypeIUCRs_PCG_Location/tyepIUCRsLocation.pdf", 
    width = 480/254, height = 480/254)
PiePlot(data = typeIUCRsLocation)
dev.off()

# exonic and intronic GO visualization ------------------------------------

## exonic
exonicMetascape <- read_xlsx(
  path = "03-results/50-TypeIUCRs_PCG_Location/exonicUCRsGenesGO/metascape_result.xlsx", 
  sheet = 2
  )
exonicMetascape <- exonicMetascape[str_detect(exonicMetascape$GroupID, ".*Summary"), ]
exonicMetascape <- exonicMetascape[c(1:10), ]

# order
exonicMetascape <- exonicMetascape[order(exonicMetascape$LogP, decreasing = TRUE), ]
exonicMetascape$Description <- factor(
  x = exonicMetascape$Description, 
  levels = exonicMetascape$Description, 
  ordered = TRUE
)

plotGO <- function(GOdata, maxPvalue){
  ggplot(data = GOdata) +
    
    geom_bar(mapping = aes(x = -(LogP), y = Description), 
             stat = "identity", width = 0.8, 
             fill = "#3C5488FF", alpha = 0.5) +
    
    geom_text(aes(x = 0.1, y = Description, label = Description), 
              size = 5, size.unit = "pt", hjust = 0) +
    
    scale_x_continuous(limits = c(0, maxPvalue), expand = c(0, 0)) +
    
    labs(x = "-Log10 pvalue", 
         y = "GO biological process", 
         title = "") +
    
    theme(
      plot.title = element_text(size = 5),
      
      panel.grid = element_blank(), 
      panel.background = element_blank(), 
      text = element_text(size = 5), 
      line = element_line(linewidth = 0.5/lwd_pt), 
      
      axis.text.y = element_blank(), 
      axis.ticks.y = element_blank(), 
      axis.line.y = element_line(linewidth = 0.5/lwd_pt, color = "black"), 
      
      axis.text.x = element_text(size = 5, color = "black"), 
      axis.line = element_line(linewidth = 0.5/lwd_pt, color = "black"), 
      
      plot.margin = unit(c(0,0.2,0,0.2), "cm")
    )
}

plotGO(GOdata = exonicMetascape, maxPvalue = 25) + canvas(width = 9.6, height = 6, units = "cm", dpi = 1200)

pdf(file = "03-results/50-TypeIUCRs_PCG_Location/exonicMetascape.pdf", 
    width = 960/254, height = 600/254)
plotGO(GOdata = exonicMetascape, maxPvalue = 25)
dev.off()

## intronic
intronicMetascape <- read_xlsx(
  path = "03-results/50-TypeIUCRs_PCG_Location/intronicUCRsGenesGO/metascape_result.xlsx", 
  sheet = 2
)
intronicMetascape <- intronicMetascape[str_detect(intronicMetascape$GroupID, ".*Summary"), ]

# order
intronicMetascape <- intronicMetascape[order(intronicMetascape$LogP, decreasing = TRUE), ]
intronicMetascape$Description <- factor(
  x = intronicMetascape$Description, 
  levels = intronicMetascape$Description, 
  ordered = TRUE
)

plotGO(intronicMetascape, maxPvalue = 8)

pdf(file = "03-results/50-TypeIUCRs_PCG_Location/intronicMetascape.pdf", 
    width = 960/254, height = 600/254)
plotGO(GOdata = intronicMetascape, maxPvalue = 8)
dev.off()

p1 <- plotGO(exonicMetascape, maxPvalue = 25)
p2 <- plotGO(intronicMetascape, maxPvalue = 8)

cowplot::plot_grid(p1, p2, 
                   align = c("hv"), axis = c("tb"), 
                   nrow = 1, ncol = 2, rel_widths = 1, rel_heights = 1, 
                   labels = c("a", "b"), label_size = 8)

pdf(file = "03-results/50-TypeIUCRs_PCG_Location/mergedMetascape.pdf", 
    width = 1800/254, height = 600/254)
cowplot::plot_grid(p1, p2, 
                   align = c("hv"), axis = c("tb"), 
                   nrow = 1, ncol = 2, rel_widths = 1, rel_heights = 1, 
                   labels = c("a", "b"), label_size = 8)
dev.off()
