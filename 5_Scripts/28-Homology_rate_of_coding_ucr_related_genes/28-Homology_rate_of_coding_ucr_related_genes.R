# 检查一下coding UCR相关基因的在人类，大鼠和小鼠当中是否是同源基因
# 是：UCR确保了基因的功能
# 否：可能是UCR本身的功能

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(dbplyr)
library(stringi)

# mouse -------------------------------------------------------------------

# mouse：mouse genome blastn to ucr(sseq是ucr序列)
mouse_to_ucr <- read.table(file = "01-data/26-Evolution_newly_emerging_UCR_distribution/mouse_results.xls", 
                    sep = "\t")
colnames(mouse_to_ucr) <- c("qseqid", "sseqid", "pident", "length", 
                     "mismatch", "gapopen", 
                     "qstart", "qend", "sstart", "send", 
                     "qseq", "sseq", "evalue", "bitscore")
mouse_to_ucr <- mouse_to_ucr[c(mouse_to_ucr$sstart == 1 | mouse_to_ucr$send == 1), ]   # 和下面一行一起组成：能匹配到ucr全长
mouse_to_ucr <- mouse_to_ucr[c(mouse_to_ucr$length == (abs(mouse_to_ucr$send - mouse_to_ucr$sstart)) + 1), ]
mouse_to_ucr <- mouse_to_ucr[c(mouse_to_ucr$pident == 100), ]
mouse_to_ucr$sseqid <- gsub("::.*", "", mouse_to_ucr$sseqid)
mouse_to_ucr <- mouse_to_ucr[!(c(mouse_to_ucr$sseqid == "uc.18" | mouse_to_ucr$sseqid == "uc.304")), ]
mouse_to_ucr_bed <- mouse_to_ucr[, c(2, 9, 10, 12)]
mouse_to_ucr_bed <- mouse_to_ucr_bed %>% 
  mutate(sseq = ifelse(sstart > send, stri_reverse(sseq), sseq))

# Code:
for (i in 1:nrow(mouse_to_ucr_bed)) {
  if (mouse_to_ucr_bed[i, 2] > mouse_to_ucr_bed[i, 3]) {
    temp <- gsub("A", "X", mouse_to_ucr_bed[i, 4])
    temp <- gsub("T", "A", temp)
    temp <- gsub("G", "Y", temp)
    temp <- gsub("C", "G", temp)
    temp <- gsub("X", "T", temp)
    temp <- gsub("Y", "C", temp)
    mouse_to_ucr_bed[i, 4] <- temp
  }
}

mouse_to_ucr_bed <- mouse_to_ucr_bed[, c(1, 4)]

# ucr blastn to mouse(qseq是ucr序列)
ucr_to_mouse <- 
  read.table(file = "01-data/28-Homology_rate_of_coding_ucr_related_genes/UCR_balstn_to_genomes/ucr_blastn_to_mouse_results.xls")

colnames(ucr_to_mouse) <- c("qseqid", "sseqid", "pident", "length", 
                            "mismatch", "gapopen", 
                            "qstart", "qend", "sstart", "send", 
                            "qseq", "sseq", "evalue", "bitscore")
ucr_to_mouse <- ucr_to_mouse[c(ucr_to_mouse$qstart == 1), ]   # 和下面一行一起组成：能匹配到ucr全长
ucr_to_mouse <- ucr_to_mouse[c(ucr_to_mouse$length == (abs(ucr_to_mouse$qend - ucr_to_mouse$qstart)) + 1), ]
ucr_to_mouse <- ucr_to_mouse[c(ucr_to_mouse$pident == 100), ]
ucr_to_mouse_bed <- ucr_to_mouse[, c(2, 9, 10, 11)]

# UCR mouse bed
ucr_mouse_bed <- merge(ucr_to_mouse_bed, mouse_to_ucr_bed, 
                       by.x = "qseq", by.y = "sseq")
ucr_mouse_bed <- ucr_mouse_bed[, c(2:5)]
colnames(ucr_mouse_bed) <- c("chr", "start", "end", "ucr_name")

ucr_mouse_bed$chr <- gsub("gb", "", ucr_mouse_bed$chr)
ucr_mouse_bed$chr <- str_sub(ucr_mouse_bed$chr, 2, 11)
save(ucr_mouse_bed, 
     file = "01-data/28-Homology_rate_of_coding_ucr_related_genes/ucr_mouse_bed.Rdata")

write.table(ucr_mouse_bed, 
            file = "01-data/28-Homology_rate_of_coding_ucr_related_genes/ucr_mouse_.bed")

# 替换染色体表示
rm(list = ls())

load("D:/R_project/UCR_project/01-data/28-Homology_rate_of_coding_ucr_related_genes/ucr_mouse_bed.Rdata")

mouse_chr_core <- read.table(
  file = "01-data/28-Homology_rate_of_coding_ucr_related_genes/mouse_chr_correspondence.tsv", 
  sep = "\t", header = TRUE)
mouse_chr_core <- mouse_chr_core[c(1:21), ]
mouse_chr_core <- mouse_chr_core[, c(7, 4)]

ucr_mouse_bed <- merge(ucr_mouse_bed, mouse_chr_core, 
                       by.x = "chr", by.y = "GenBank.seq.accession")
ucr_mouse_bed <- ucr_mouse_bed[, c(5, 2:4)]
save(ucr_mouse_bed, 
     file = "01-data/28-Homology_rate_of_coding_ucr_related_genes/ucr_mouse_simple_chrname_bed.Rdata")

# Code:
for (i in 1:nrow(ucr_mouse_bed)) {
  if (ucr_mouse_bed[i, 2] > ucr_mouse_bed[i, 3]) {
    temp <- ucr_mouse_bed[i, 2]
    ucr_mouse_bed[i, 2] <- ucr_mouse_bed[i, 3]
    ucr_mouse_bed[i, 3] <- temp
  }
}

write.table(ucr_mouse_bed, 
            file = "01-data/28-Homology_rate_of_coding_ucr_related_genes/ucr_mouse.bed", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# human mouse homology ----------------------------------------------------------------

# 在linux使用bedtools intersect获得mouse_coding_ucr

rm(list = ls())

library(homologene)

mouse_coding_UCR <- read.table(
  file = "02-analysis/28-Homology_rate_of_coding_ucr_related_genes/mouse_coding_ucr.bed")
mouse_coding_UCR <- mouse_coding_UCR[, c(4, 8:9)]

human_coding_UCR <- read.table(
  file = "02-analysis/16-New_classification/01-UCR_from_protein_coding_gene.bed")
human_coding_UCR <- human_coding_UCR[, c(4, 8:9)]

# 用homologene获得小鼠的所有的同源基因，查看多少基因的同源基因在human_coding_ucr_genes里
homologous_genes_of_mouse_genes <- homologene(mouse_coding_UCR$V9, 
                                              inTax = 10090, 
                                              outTax = 9606)

intersect(human_coding_UCR$V9, homologous_genes_of_mouse_genes$'9606')   #181个
length(unique(human_coding_UCR$V9))   # 201个
# 也就是说201个基因里面有181个基因和小鼠是同源基因

# rat ---------------------------------------------------------------------

rm(list = ls())

rat_to_ucr <- read.table(
  file = "01-data/26-Evolution_newly_emerging_UCR_distribution/rat_results.xls", 
  sep = "\t"
)

colnames(rat_to_ucr) <- c("qseqid", "sseqid", "pident", "length", 
                            "mismatch", "gapopen", 
                            "qstart", "qend", "sstart", "send", 
                            "qseq", "sseq", "evalue", "bitscore")
rat_to_ucr <- rat_to_ucr[c(rat_to_ucr$sstart == 1 | rat_to_ucr$send == 1), ]   # 和下面一行一起组成：能匹配到ucr全长
rat_to_ucr <- rat_to_ucr[c(rat_to_ucr$length == (abs(rat_to_ucr$send - rat_to_ucr$sstart)) + 1), ]
rat_to_ucr <- rat_to_ucr[c(rat_to_ucr$pident == 100), ]
rat_to_ucr$sseqid <- gsub("::.*", "", rat_to_ucr$sseqid)
rat_to_ucr <- rat_to_ucr[!(c(rat_to_ucr$sseqid == "uc.18" | rat_to_ucr$sseqid == "uc.304")), ]

rat_to_ucr <- rat_to_ucr[!(c(rat_to_ucr$sseqid == "uc.101" & rat_to_ucr$sstart == 60)), ]
rat_to_ucr <- rat_to_ucr[!(c(rat_to_ucr$sseqid == "uc.373" & rat_to_ucr$sstart == 73)), ]

rat_to_ucr_bed <- rat_to_ucr[, c(2, 9, 10, 12)]
rat_to_ucr_bed <- rat_to_ucr_bed %>% 
  mutate(sseq = ifelse(sstart > send, stri_reverse(sseq), sseq))

# Code:
for (i in 1:nrow(rat_to_ucr_bed)) {
  if (rat_to_ucr_bed[i, 2] > rat_to_ucr_bed[i, 3]) {
    temp <- gsub("A", "X", rat_to_ucr_bed[i, 4])
    temp <- gsub("T", "A", temp)
    temp <- gsub("G", "Y", temp)
    temp <- gsub("C", "G", temp)
    temp <- gsub("X", "T", temp)
    temp <- gsub("Y", "C", temp)
    rat_to_ucr_bed[i, 4] <- temp
  }
}

rat_to_ucr_bed <- rat_to_ucr_bed[, c(1, 4)]

# ucr blastn to rat(qseq是ucr序列)
ucr_to_rat <- 
  read.table(file = "01-data/28-Homology_rate_of_coding_ucr_related_genes/UCR_balstn_to_genomes/ucr_blastn_to_rat_results.xls")

colnames(ucr_to_rat) <- c("qseqid", "sseqid", "pident", "length", 
                            "mismatch", "gapopen", 
                            "qstart", "qend", "sstart", "send", 
                            "qseq", "sseq", "evalue", "bitscore")
ucr_to_rat <- ucr_to_rat[c(ucr_to_rat$qstart == 1), ]   # 和下面一行一起组成：能匹配到ucr全长
ucr_to_rat <- ucr_to_rat[c(ucr_to_rat$length == (abs(ucr_to_rat$qend - ucr_to_rat$qstart)) + 1), ]
ucr_to_rat <- ucr_to_rat[c(ucr_to_rat$pident == 100), ]
ucr_to_rat_bed <- ucr_to_rat[, c(2, 9, 10, 11)]

# UCR rat bed
ucr_rat_bed <- merge(ucr_to_rat_bed, rat_to_ucr_bed, 
                       by.x = "qseq", by.y = "sseq")
ucr_rat_bed <- ucr_rat_bed[, c(2:5)]
colnames(ucr_rat_bed) <- c("chr", "start", "end", "ucr_name")

ucr_rat_bed$chr <- gsub("gb", "", ucr_rat_bed$chr)
ucr_rat_bed$chr <- str_sub(ucr_rat_bed$chr, 2, 11)
save(ucr_rat_bed, 
     file = "01-data/28-Homology_rate_of_coding_ucr_related_genes/ucr_rat_bed.Rdata")

write.table(ucr_rat_bed, 
            file = "01-data/28-Homology_rate_of_coding_ucr_related_genes/ucr_rat_.bed")

# 替换染色体表示
rm(list = ls())

load("D:/R_project/UCR_project/01-data/28-Homology_rate_of_coding_ucr_related_genes/ucr_rat_bed.Rdata")

rat_chr_core <- read.table(
  file = "01-data/28-Homology_rate_of_coding_ucr_related_genes/rat_chr_correspondence.tsv", 
  sep = "\t", header = TRUE)
rat_chr_core <- rat_chr_core[c(1:22), ]
rat_chr_core <- rat_chr_core[, c(7, 4)]

ucr_rat_bed <- merge(ucr_rat_bed, rat_chr_core, 
                       by.x = "chr", by.y = "GenBank.seq.accession")
ucr_rat_bed <- ucr_rat_bed[, c(5, 2:4)]
save(ucr_rat_bed, 
     file = "01-data/28-Homology_rate_of_coding_ucr_related_genes/ucr_rat_simple_chrname_bed.Rdata")

# Code:
for (i in 1:nrow(ucr_rat_bed)) {
  if (ucr_rat_bed[i, 2] > ucr_rat_bed[i, 3]) {
    temp <- ucr_rat_bed[i, 2]
    ucr_rat_bed[i, 2] <- ucr_rat_bed[i, 3]
    ucr_rat_bed[i, 3] <- temp
  }
}

write.table(ucr_rat_bed, 
            file = "01-data/28-Homology_rate_of_coding_ucr_related_genes/ucr_rat.bed", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# human rat homology ------------------------------------------------------

# 在Linux使用bedtools intersect获得rat_coding_ucr

rm(list = ls())

library(homologene)

rat_coding_UCR <- read.table(
  file = "02-analysis/28-Homology_rate_of_coding_ucr_related_genes/rat_coding_ucr.bed")
rat_coding_UCR <- rat_coding_UCR[, c(4, 8:9)]

human_coding_UCR <- read.table(
  file = "02-analysis/16-New_classification/01-UCR_from_protein_coding_gene.bed")
human_coding_UCR <- human_coding_UCR[, c(4, 8:9)]

# 用homologene获得大鼠的所有的同源基因，查看多少基因的同源基因在human_coding_ucr_genes里
homologous_genes_of_rat_genes <- homologene(rat_coding_UCR$V9, 
                                              inTax = 10116, 
                                              outTax = 9606)

intersect(human_coding_UCR$V9, homologous_genes_of_rat_genes$'9606')   #155个
length(unique(human_coding_UCR$V9))   # 201个
# 也就是说201个基因里面有155个基因和小鼠是同源基因

# 绘制一张韦恩图，展示human rat mouse同源的部分 ------------------------------------------





