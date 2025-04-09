# http://enhanceratlas.org/
# 需要intergenic的bed文件

rm(list = ls())

library(tidyverse)

intergenic <- read.table(file = "02-analysis/16-New_classification/03-UCR_intergenic.bed", 
                         sep = "\t", header = FALSE)
intergenic$V1 <- paste("chr", intergenic$V1, sep = "")
write.table(intergenic, file = "02-analysis/11-intergenic_region_enhancer_site/01-intergenic(chr).bed", sep = "\t", 
            col.names = F, row.names = F, quote = F)





























