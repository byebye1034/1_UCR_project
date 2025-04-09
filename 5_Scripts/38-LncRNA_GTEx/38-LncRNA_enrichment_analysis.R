# Immune-related 3-lncRNA signature with prognostic connotation in a multi-cancer setting
# 文章里使用LncRNAs2Pathways进行富集分析，进行了一些调整，我看看我的富集分析结果是否需要调整

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(ggplot2)
library(plyr)
library(LncPath)
library(stringr)

# Set parameters
version = "GRCh38"

# get lncRNA-mRNA interaction network
NetLncPath <- getNet()

length(unique(NetLncPath$V1))
length(unique(NetLncPath$V2))

# ncRNA UCR重叠的lncRNA
ncRNA_UCR_lncRNA <- read.table(file = "02-analysis/16-New_classification/ncRNA_UCR(66)_ensembl_id.txt", sep = "\t")
ncRNA_UCR_lncRNA <- unique(ncRNA_UCR_lncRNA)

SigLncs <- as.character(unlist(ncRNA_UCR_lncRNA))
print(SigLncs)

#Perform the Walkscore to get score for each gene (pcg/mRNA)
if (length(SigLncs) == 0) 
  stop("The list is empty.")
cat("Now start the random walking...\n")
if (!exists("LncPathEnvir")) 
  LncPathEnvir <- initializeLncPathEnvir()
Network <- NetLncPath
NetLncPath <- graph_from_edgelist(as.matrix(Network), directed = FALSE)
VertexWeight <- rep(0, length(V(NetLncPath)))
names(VertexWeight) <- V(NetLncPath)$name

VertexWeight[V(NetLncPath)$name %in% SigLncs] <- 1
print(names(which(VertexWeight > 0)))
cat(paste(length(which(VertexWeight > 0)), "of", length(SigLncs), 
          "were mapped to the huge net with", length(VertexWeight), 
          "Nodes.", "\n"))

WalkRes <- RandomWalk2igraph(NetLncPath, VertexWeight, EdgeWeight = FALSE)

# save(WalkRes, file = "02-analysis/38-LncRNA_GTEx/WalkRes.Rdata")
# save(VertexWeight, file = "02-analysis/38-LncRNA_GTEx/VertexWeight.Rdata")

WalkScore <- as.data.frame(cbind(V(NetLncPath)$name, WalkRes))
colnames(WalkScore) <- c("Gene","WalkRes")
WalkScore$Gene <- as.character(as.vector(WalkScore$Gene))
WalkScore$WalkRes <- as.numeric(as.vector(WalkScore$WalkRes))

if (length(which(as.numeric(WalkScore[[2]]) == 0)) > 0) 
  WalkScore <- WalkScore[-(which(as.numeric(WalkScore[[2]]) == 
                                   0)), ]
WalkScore[[2]] <- sqrt(as.numeric(WalkScore[[2]]))
WalkScore <- WalkScore[order(WalkScore[[2]], decreasing = TRUE), 
]
PCEnsem2Sym <- get("PCEnsem2Sym", envir = LncPathEnvir)
PCWalkScore <- merge(PCEnsem2Sym, WalkScore, by.x = 2, by.y = 1)
PCWalkScore <- PCWalkScore[order(PCWalkScore[[3]], decreasing = TRUE), 
]

differntial_PCWalkScore <- PCWalkScore[PCWalkScore$WalkRes>=0.01,]















