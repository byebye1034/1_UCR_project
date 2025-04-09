# get human/mouse type III UCRs nearest upstream and downstream protein coding genes

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library(readxl)
library(stringr)
library(ggplot2)
library(ggview)
library(cowplot)

mouseUCRs <- read.table(
  file = "02-analysis/47-Other_UCR_nearest_PCGs/ucr_mouse.bed", sep = "\t"
)

# type I UCRs
mouseTypeIUCRs <- read.table(
  file = "02-analysis/47-Other_UCR_nearest_PCGs/mousetypeIUCRs.bed", sep = "\t"
)
length(unique(mouseTypeIUCRs$V4)) # 245 type I UCRs
length(unique(mouseTypeIUCRs$V8)) # 157 PCGs

# type II UCRs
mouseTypeIIUCRs <- read.table(
  file = "02-analysis/47-Other_UCR_nearest_PCGs/mousetypeIIUCRs.bed", sep = "\t"
)
length(unique(mouseTypeIIUCRs$V4)) # 35 type I UCRs
length(unique(mouseTypeIIUCRs$V8)) # 28 ncRNA
table(mouseTypeIIUCRs$V10) # 15 antisense_RNA 20 lincRNA

nonTypeIIIUCRs <- c(unique(mouseTypeIUCRs$V4), 
                    unique(mouseTypeIIUCRs$V4))
nonTypeIIIUCRs <- unique(nonTypeIIIUCRs) # 7 mouse UCR overlap with PCGs and ncRNA

mouseTypeIIIUCRs <- mouseUCRs[!(mouseUCRs$V4 %in% nonTypeIIIUCRs), ] # 206 mouse type III UCRs
write.table(mouseTypeIIIUCRs, 
            file = "02-analysis/47-Other_UCR_nearest_PCGs/mouseTypeIIIUCRs.bed", 
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

# get nearest up and down PCGs on server
mouseUpPCGs <- read.table(
  file = "02-analysis/47-Other_UCR_nearest_PCGs/mouseTypeIIIUCRsUpPCGs.bed", sep = "\t"
)
mouseUpPCGs <- mouseUpPCGs[, c(4, 8, 9, 11)]
mouseDownPCGs <- read.table(
  file = "02-analysis/47-Other_UCR_nearest_PCGs/mouseTypeIIIUCRsDownPCGs.bed", sep = "\t"
)
mouseDownPCGs <- mouseDownPCGs[, c(4, 8, 9, 11)]
mouseIIIPCGs <- merge(
  mouseUpPCGs, mouseDownPCGs, 
  by.x = "V4", by.y = "V4"
)
mouseIIIPCGs <- mouseIIIPCGs %>% 
  mutate(Nearest = if_else(V11.x + V11.y < 0, V9.y, V9.x))
mouseIIIPCGs <- mouseIIIPCGs %>% 
  mutate(Distance = if_else(V11.x + V11.y < 0, V11.y, abs(V11.x)))

write.csv(mouseIIIPCGs, 
          file = "02-analysis/47-Other_UCR_nearest_PCGs/mouseTypeIIIUCRsPCGs.csv", 
          row.names = FALSE)

lwd_pt <- .pt*72.27/96

p1 <- ggplot(
  data = mouseIIIPCGs, 
  mapping = aes(x = log10(Distance))
) +
  
  geom_density(linewidth = 0.5/lwd_pt, color = "#4DBBD5FF") +
  geom_rug(linewidth = 0.5/lwd_pt, color = "#4DBBD5FF") +
  
  geom_vline(xintercept = log10(2000), 
             linetype = "dashed", 
             color = "gray", 
             linewidth = 0.5/lwd_pt) +
  
  scale_y_continuous(expand = c(0, 0)) +
  
  labs(
    x = "log10(distance between
    type III UCRs and PCGs)", 
    y = "Density", 
    title = "Mouse"
  ) +
  
  theme(
    panel.background = element_blank(), 
    panel.grid = element_blank(), 
    line = element_line(linewidth = 0.5/lwd_pt), 
    
    axis.text = element_text(size = 7, color = "#000000"), 
    axis.ticks = element_line(linewidth = 0.5/lwd_pt), 
    axis.line = element_line(linewidth = 0.5/lwd_pt), 
    axis.title = element_text(size = 7, color = "#000000"), 
    
    plot.title = element_text(size = 7),
    
    aspect.ratio = 1:1
  )
p1   # 距离主要在10000，一万bp以上

p1 + canvas(width = 4.8, height = 6, units = "cm", dpi = 1200)



humanTypeIIIUCRs <- read.table(
  file = "02-analysis/16-New_classification/03-UCR_intergenic.bed", sep = "\t"
)
setdiff(humanTypeIIIUCRs$V4, mouseTypeIIIUCRs$V4)

humanTypeIUCRs <- read.table(
  file = "02-analysis/47-Other_UCR_nearest_PCGs/human/humanTypeIUCRs.bed", sep = "\t"
)
length(unique(humanTypeIUCRs$V4))

humanTypeIUCRs1 <- read.table(
  file = "02-analysis/16-New_classification/01-UCR_from_protein_coding_gene.bed", sep = "\t"
)

humanTypeIIUCRs1 <- read.table(
  file = "02-analysis/47-Other_UCR_nearest_PCGs/human/humanTypeIIUCRs1.bed", sep = "\t"
)
humanTypeIIUCRs2 <- read.table(
  file = "02-analysis/47-Other_UCR_nearest_PCGs/human/humanTypeIIUCRs2.bed", sep = "\t"
)
humanTypeIIUCRsNew <- c(humanTypeIIUCRs2$V4, humanTypeIIUCRs1$V4)
humanTypeIIUCRsNew <- unique(humanTypeIIUCRsNew)
humanTypeIIUCRsNew1 <- setdiff(humanTypeIIUCRsNew, humanTypeIUCRs$V4)

humanTypeIIUCRs <- read.table(
  file = "02-analysis/16-New_classification/02-UCR_from_Non_Coding_RNA(66).bed", sep = "\t"
)
setdiff(humanTypeIIUCRsNew1, humanTypeIIUCRs$V4) # uc.304

## 确认human type III UCRs nearest PCGs没错
humanIIIUp <- read.table(
  file = "02-analysis/47-Other_UCR_nearest_PCGs/human/humanTypeIIIUCRsUpPCGs.bed", sep = "\t"
)
humanIIIUp <- humanIIIUp[, c(4, 8, 9, 11)]
humanIIIDown <- read.table(
  file = "02-analysis/47-Other_UCR_nearest_PCGs/human/humanTypeIIIUCRsDownPCGs.bed", sep = "\t"
)
humanIIIDown <- humanIIIDown[, c(4, 8, 9, 11)]
humanIIIPCGs <- merge(
  humanIIIUp, humanIIIDown, 
  by.x = "V4", by.y = "V4"
)
humanIIIPCGs <- humanIIIPCGs %>% 
  mutate(Nearest = if_else(V11.x + V11.y < 0, V9.y, V9.x))
humanIIIPCGs <- humanIIIPCGs %>% 
  mutate(Distance = if_else(V11.x + V11.y < 0, V11.y, abs(V11.x)))

write.csv(humanIIIPCGs, 
          file = "humanIIIPCGs.csv", row.names = FALSE)

lwd_pt <- .pt*72.27/96

p2 <- ggplot(
  data = humanIIIPCGs, 
  mapping = aes(x = log10(Distance))
) +
  
  geom_density(linewidth = 0.5/lwd_pt, color = "#E64B35FF") +
  geom_rug(linewidth = 0.5/lwd_pt, color = "#E64B35FF") +
  
  geom_vline(xintercept = log10(2000), 
             linetype = "dashed", 
             color = "gray", 
             linewidth = 0.5/lwd_pt) +
  
  scale_y_continuous(expand = c(0, 0)) +
  
  labs(
    x = "log10(distance between
    type III UCRs and PCGs)", 
    y = "Density", 
    title = "Human"
  ) +
  
  theme(
    panel.background = element_blank(), 
    panel.grid = element_blank(), 
    line = element_line(linewidth = 0.5/lwd_pt), 
    
    axis.text = element_text(size = 7, color = "#000000"), 
    axis.ticks = element_line(linewidth = 0.5/lwd_pt), 
    axis.line = element_line(linewidth = 0.5/lwd_pt), 
    axis.title = element_text(size = 7, color = "#000000"), 
    
    plot.title = element_text(size = 7),
    
    aspect.ratio = 1:1
  )
p2   # 距离主要在10000，一万bp以上

p2 + canvas(width = 4.8, height = 6, units = "cm", dpi = 1200)

cowplot::plot_grid(p2, p1, 
                   align = "hv", axis = c("tb"), 
                   ncol = 2, nrow = 1, rel_widths = 1, rel_heights = 1, 
                   labels = c("b", "c"), label_size = 8)

pdf(file = "03-results/47-Other_UCR_nearest_PCGs/distanceBtTypeIIIUCRsPCGs.pdf", 
    width = 960/254, height = 600/254)

dev.off()
