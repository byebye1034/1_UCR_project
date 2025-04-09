# 统计在不同物种当中新出现的UCR

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(stringr)

# human -------------------------------------------------------------------

human <- read.table(file = "01-data/26-Evolution_newly_emerging_UCR_distribution/human_results.xls", 
                    sep = "\t")
colnames(human) <- c("qseqid", "sseqid", "pident", "length", 
                     "mismatch", "gapopen", 
                     "qstart", "qend", "sstart", "send", 
                     "qseq", "sseq", "evalue", "bitscore")
human <- human[c(human$sstart == 1), ]   # 和下面一行一起组成：能匹配到ucr全长
human <- human[c(human$length == human$send), ]
human <- human[c(human$pident == 100), ]
human$sseqid <- gsub("::.*", "", human$sseqid)
human <- human[!(c(human$sseqid == "uc.18" | human$sseqid == "uc.304")), ]

# 使用table函数统计每个数值的频次
freq_sseqid <- table(human$sseqid)

# 筛选出出现超过两次的数值
result <- names(freq_sseqid[freq_sseqid > 1])
print(result)
# [1] "uc.186" "uc.410" "uc.411" "uc.412" "uc.43"  "uc.458" "uc.96"  "uc.97"  "uc.98"  "uc.99"

human <- human[, c(2, 9, 10, 12)]   # 现在不需要统计匹配到哪里，只要有100%匹配的位置就可以
human <- unique(human)

save(human, file = "02-analysis/26-Evolution_newly_emerging_UCR_distribution/human_UCR.Rdata")

# mouse -------------------------------------------------------------------

rm(list = ls())

mouse <- read.table(file = "01-data/26-Evolution_newly_emerging_UCR_distribution/mouse_results.xls", 
                        sep = "\t")
colnames(mouse) <- c("qseqid", "sseqid", "pident", "length", 
                         "mismatch", "gapopen", 
                         "qstart", "qend", "sstart", "send", 
                         "qseq", "sseq", "evalue", "bitscore")
mouse <- mouse[c(mouse$sstart == 1 | mouse$send == 1), ]   # 和下面一行一起组成：能匹配到ucr全长
mouse <- mouse[c(mouse$length == (abs(mouse$send - mouse$sstart)) + 1), ]
mouse <- mouse[c(mouse$pident == 100), ]
mouse$sseqid <- gsub("::.*", "", mouse$sseqid)
mouse <- mouse[!(c(mouse$sseqid == "uc.18" | mouse$sseqid == "uc.304")), ]

# 使用table函数统计每个数值的频次
freq_sseqid <- table(mouse$sseqid)

# 筛选出出现超过两次的数值
result <- names(freq_sseqid[freq_sseqid > 1])
print(result)
# [1] "uc.333" "uc.480"

mouse <- mouse[, c(2, 9, 10, 12)]   # 现在不需要统计匹配到哪里，只要有100%匹配的位置就可以
mouse <- unique(mouse)
length(unique(mouse$sseqid))
# 479 UCR

save(mouse, file = "02-analysis/26-Evolution_newly_emerging_UCR_distribution/mouse_UCR.Rdata")

# rat ---------------------------------------------------------------------

rm(list = ls())

rat <- read.table(file = "01-data/26-Evolution_newly_emerging_UCR_distribution/rat_results.xls", 
                    sep = "\t")
colnames(rat) <- c("qseqid", "sseqid", "pident", "length", 
                     "mismatch", "gapopen", 
                     "qstart", "qend", "sstart", "send", 
                     "qseq", "sseq", "evalue", "bitscore")
rat <- rat[c(rat$sstart == 1 | rat$send == 1), ]   # 和下面一行一起组成：能匹配到ucr全长
rat <- rat[c(rat$length == (abs(rat$send - rat$sstart)) + 1), ]
rat <- rat[c(rat$pident == 100), ]
rat$sseqid <- gsub("::.*", "", rat$sseqid)
rat <- rat[!(c(rat$sseqid == "uc.18" | rat$sseqid == "uc.304")), ]

# 使用table函数统计每个数值的频次
freq_sseqid <- table(rat$sseqid)

# 筛选出出现超过两次的数值
result <- names(freq_sseqid[freq_sseqid > 1])
print(result)
# [1] "uc.101" "uc.373"

rat <- rat[, c(2, 9, 10, 12)]   # 现在不需要统计匹配到哪里，只要有100%匹配的位置就可以
rat <- unique(rat)
length(unique(rat$sseqid))
# 479 UCR

rat <- rat[!(rat$sseqid == "uc.101" & rat$sstart == 60), ]
rat <- rat[!(rat$sseqid == "uc.373" & rat$sstart == 73), ]

save(rat, file = "02-analysis/26-Evolution_newly_emerging_UCR_distribution/rat_UCR.Rdata")

# rhemonkey ---------------------------------------------------------------

rm(list = ls())

rhemonkey <- read.table(file = "01-data/26-Evolution_newly_emerging_UCR_distribution/rhemonkey_results.xls", 
                    sep = "\t")
colnames(rhemonkey) <- c("qseqid", "sseqid", "pident", "length", 
                     "mismatch", "gapopen", 
                     "qstart", "qend", "sstart", "send", 
                     "qseq", "sseq", "evalue", "bitscore")
rhemonkey <- rhemonkey[c(rhemonkey$sstart == 1 | rhemonkey$send == 1), ]   # 和下面一行一起组成：能匹配到ucr全长
rhemonkey <- rhemonkey[c(rhemonkey$length == (abs(rhemonkey$send - rhemonkey$sstart)) + 1), ]
rhemonkey <- rhemonkey[c(rhemonkey$pident == 100), ]
rhemonkey$sseqid <- gsub("::.*", "", rhemonkey$sseqid)
rhemonkey <- rhemonkey[!(c(rhemonkey$sseqid == "uc.18" | rhemonkey$sseqid == "uc.304")), ]

# 使用table函数统计每个数值的频次
freq_sseqid <- table(rhemonkey$sseqid)

# 筛选出出现超过两次的数值
result <- names(freq_sseqid[freq_sseqid > 1])
print(result)
# [1] "uc.333" "uc.480"

rhemonkey <- rhemonkey[, c(2, 9, 10, 12)]   # 现在不需要统计匹配到哪里，只要有100%匹配的位置就可以
rhemonkey <- unique(rhemonkey)
length(unique(rhemonkey$sseqid))
# 386 UCR

save(rhemonkey, file = "02-analysis/26-Evolution_newly_emerging_UCR_distribution/rhemonkey_UCR.Rdata")

# chimpanzee --------------------------------------------------------------

rm(list = ls())

chimpanzee <- read.table(file = "01-data/26-Evolution_newly_emerging_UCR_distribution/chimpanzee_results.xls", 
                  sep = "\t")
colnames(chimpanzee) <- c("qseqid", "sseqid", "pident", "length", 
                   "mismatch", "gapopen", 
                   "qstart", "qend", "sstart", "send", 
                   "qseq", "sseq", "evalue", "bitscore")
chimpanzee <- chimpanzee[c(chimpanzee$sstart == 1 | chimpanzee$send == 1), ]   # 和下面一行一起组成：能匹配到ucr全长
chimpanzee <- chimpanzee[c(chimpanzee$length == (abs(chimpanzee$send - chimpanzee$sstart)) + 1), ]
chimpanzee <- chimpanzee[c(chimpanzee$pident == 100), ]
chimpanzee$sseqid <- gsub("::.*", "", chimpanzee$sseqid)
chimpanzee <- chimpanzee[!(c(chimpanzee$sseqid == "uc.18" | chimpanzee$sseqid == "uc.304")), ]

# 使用table函数统计每个数值的频次
freq_sseqid <- table(chimpanzee$sseqid)

# 筛选出出现超过两次的数值
result <- names(freq_sseqid[freq_sseqid > 1])
print(result)
# [1] "uc.30"

chimpanzee <- chimpanzee[, c(2, 9, 10, 12)]   # 现在不需要统计匹配到哪里，只要有100%匹配的位置就可以
chimpanzee <- unique(chimpanzee)
length(unique(chimpanzee$sseqid))
# 416 UCR

chimpanzee <- chimpanzee[!(chimpanzee$sseqid == "uc.30" & chimpanzee$sstart == 243), ]

save(chimpanzee, file = "02-analysis/26-Evolution_newly_emerging_UCR_distribution/chimpanzee_UCR.Rdata")

# sheep -------------------------------------------------------------------

rm(list = ls())

sheep <- read.table(file = "01-data/26-Evolution_newly_emerging_UCR_distribution/sheep_results.xls", 
                         sep = "\t")
colnames(sheep) <- c("qseqid", "sseqid", "pident", "length", 
                          "mismatch", "gapopen", 
                          "qstart", "qend", "sstart", "send", 
                          "qseq", "sseq", "evalue", "bitscore")
sheep <- sheep[c(sheep$sstart == 1 | sheep$send == 1), ]   # 和下面一行一起组成：能匹配到ucr全长
sheep <- sheep[c(sheep$length == (abs(sheep$send - sheep$sstart)) + 1), ]
sheep <- sheep[c(sheep$pident == 100), ]
sheep$sseqid <- gsub("::.*", "", sheep$sseqid)
sheep <- sheep[!(c(sheep$sseqid == "uc.18" | sheep$sseqid == "uc.304")), ]

# 使用table函数统计每个数值的频次
freq_sseqid <- table(sheep$sseqid)

# 筛选出出现超过两次的数值
result <- names(freq_sseqid[freq_sseqid > 1])
print(result)
# [1] ""

sheep <- sheep[, c(2, 9, 10, 12)]   # 现在不需要统计匹配到哪里，只要有100%匹配的位置就可以
sheep <- unique(sheep)
length(unique(sheep$sseqid))
# 234 UCR

save(sheep, file = "02-analysis/26-Evolution_newly_emerging_UCR_distribution/sheep_UCR.Rdata")

# horse -------------------------------------------------------------------

rm(list = ls())

horse <- read.table(file = "01-data/26-Evolution_newly_emerging_UCR_distribution/horse_results.xls", 
                         sep = "\t")
colnames(horse) <- c("qseqid", "sseqid", "pident", "length", 
                          "mismatch", "gapopen", 
                          "qstart", "qend", "sstart", "send", 
                          "qseq", "sseq", "evalue", "bitscore")
horse <- horse[c(horse$sstart == 1 | horse$send == 1), ]   # 和下面一行一起组成：能匹配到ucr全长
horse <- horse[c(horse$length == (abs(horse$send - horse$sstart)) + 1), ]
horse <- horse[c(horse$pident == 100), ]
horse$sseqid <- gsub("::.*", "", horse$sseqid)
horse <- horse[!(c(horse$sseqid == "uc.18" | horse$sseqid == "uc.304")), ]

# 使用table函数统计每个数值的频次
freq_sseqid <- table(horse$sseqid)

# 筛选出出现超过两次的数值
result <- names(freq_sseqid[freq_sseqid > 1])
print(result)
# [1] "uc.242" "uc.401"

horse <- horse[, c(2, 9, 10, 12)]   # 现在不需要统计匹配到哪里，只要有100%匹配的位置就可以
horse <- unique(horse)
length(unique(horse$sseqid))
# 257
# horse里uc.242可以匹配到两个位置，uc.401可以匹配到4个位置

save(horse, file = "02-analysis/26-Evolution_newly_emerging_UCR_distribution/horse_UCR.Rdata")

# mochpri -----------------------------------------------------------------

rm(list = ls())

mochpri <- read.table(file = "01-data/26-Evolution_newly_emerging_UCR_distribution/mochpri_results.xls", 
                    sep = "\t")
colnames(mochpri) <- c("qseqid", "sseqid", "pident", "length", 
                     "mismatch", "gapopen", 
                     "qstart", "qend", "sstart", "send", 
                     "qseq", "sseq", "evalue", "bitscore")
mochpri <- mochpri[c(mochpri$sstart == 1 | mochpri$send == 1), ]   # 和下面一行一起组成：能匹配到ucr全长
mochpri <- mochpri[c(mochpri$length == (abs(mochpri$send - mochpri$sstart)) + 1), ]
mochpri <- mochpri[c(mochpri$pident == 100), ]
mochpri$sseqid <- gsub("::.*", "", mochpri$sseqid)
mochpri <- mochpri[!(c(mochpri$sseqid == "uc.18" | mochpri$sseqid == "uc.304")), ]

# 使用table函数统计每个数值的频次
freq_sseqid <- table(mochpri$sseqid)

# 筛选出出现超过两次的数值
result <- names(freq_sseqid[freq_sseqid > 1])
print(result)
# [1] "uc.242" "uc.401"

mochpri <- mochpri[, c(2, 9, 10, 12)]   # 现在不需要统计匹配到哪里，只要有100%匹配的位置就可以
mochpri <- unique(mochpri)
length(unique(mochpri$sseqid))
# 155 UCR

save(mochpri, file = "02-analysis/26-Evolution_newly_emerging_UCR_distribution/mochpri_UCR.Rdata")

# cow ---------------------------------------------------------------------

rm(list = ls())

cow <- read.table(file = "01-data/26-Evolution_newly_emerging_UCR_distribution/cow_results.xls", 
                      sep = "\t")
colnames(cow) <- c("qseqid", "sseqid", "pident", "length", 
                       "mismatch", "gapopen", 
                       "qstart", "qend", "sstart", "send", 
                       "qseq", "sseq", "evalue", "bitscore")
cow <- cow[c(cow$sstart == 1 | cow$send == 1), ]   # 和下面一行一起组成：能匹配到ucr全长
cow <- cow[c(cow$length == (abs(cow$send - cow$sstart)) + 1), ]
cow <- cow[c(cow$pident == 100), ]
cow$sseqid <- gsub("::.*", "", cow$sseqid)
cow <- cow[!(c(cow$sseqid == "uc.18" | cow$sseqid == "uc.304")), ]

# 使用table函数统计每个数值的频次
freq_sseqid <- table(cow$sseqid)

# 筛选出出现超过两次的数值
result <- names(freq_sseqid[freq_sseqid > 1])
print(result)
# [1] "uc.471"

cow <- cow[, c(2, 9, 10, 12)]   # 现在不需要统计匹配到哪里，只要有100%匹配的位置就可以
cow <- unique(cow)
length(unique(cow$sseqid))
# 235 UCR

cow <- cow[!(cow$sseqid == "uc.471" & cow$sstart == 205), ]

save(cow, file = "02-analysis/26-Evolution_newly_emerging_UCR_distribution/cow_UCR.Rdata")

# pig ---------------------------------------------------------------------

rm(list = ls())

pig <- read.table(file = "01-data/26-Evolution_newly_emerging_UCR_distribution/pig_results.xls", 
                  sep = "\t")
colnames(pig) <- c("qseqid", "sseqid", "pident", "length", 
                   "mismatch", "gapopen", 
                   "qstart", "qend", "sstart", "send", 
                   "qseq", "sseq", "evalue", "bitscore")
pig <- pig[c(pig$sstart == 1 | pig$send == 1), ]   # 和下面一行一起组成：能匹配到ucr全长
pig <- pig[c(pig$length == (abs(pig$send - pig$sstart)) + 1), ]
pig <- pig[c(pig$pident == 100), ]
pig$sseqid <- gsub("::.*", "", pig$sseqid)
pig <- pig[!(c(pig$sseqid == "uc.18" | pig$sseqid == "uc.304")), ]

# 使用table函数统计每个数值的频次
freq_sseqid <- table(pig$sseqid)

# 筛选出出现超过两次的数值
result <- names(freq_sseqid[freq_sseqid > 1])
print(result)
# [1] ""

pig <- pig[, c(2, 9, 10, 12)]   # 现在不需要统计匹配到哪里，只要有100%匹配的位置就可以
pig <- unique(pig)
length(unique(pig$sseqid))
# 263 UCR

# pig <- pig[!(pig$sseqid == "uc.471" & pig$sstart == 205), ]

save(pig, file = "02-analysis/26-Evolution_newly_emerging_UCR_distribution/pig_UCR.Rdata")

# cat ---------------------------------------------------------------------

rm(list = ls())

cat <- read.table(file = "01-data/26-Evolution_newly_emerging_UCR_distribution/cat_results.xls", 
                  sep = "\t")
colnames(cat) <- c("qseqid", "sseqid", "pident", "length", 
                   "mismatch", "gapopen", 
                   "qstart", "qend", "sstart", "send", 
                   "qseq", "sseq", "evalue", "bitscore")
cat <- cat[c(cat$sstart == 1 | cat$send == 1), ]   # 和下面一行一起组成：能匹配到ucr全长
cat <- cat[c(cat$length == (abs(cat$send - cat$sstart)) + 1), ]
cat <- cat[c(cat$pident == 100), ]
cat$sseqid <- gsub("::.*", "", cat$sseqid)
cat <- cat[!(c(cat$sseqid == "uc.18" | cat$sseqid == "uc.304")), ]

# 使用table函数统计每个数值的频次
freq_sseqid <- table(cat$sseqid)

# 筛选出出现超过两次的数值
result <- names(freq_sseqid[freq_sseqid > 1])
print(result)
# [1] "“

cat <- cat[, c(2, 9, 10, 12)]   # 现在不需要统计匹配到哪里，只要有100%匹配的位置就可以
cat <- unique(cat)
length(unique(cat$sseqid))
# 281 UCR

# cat <- cat[!(cat$sseqid == "uc.471" & cat$sstart == 205), ]

save(cat, file = "02-analysis/26-Evolution_newly_emerging_UCR_distribution/cat_UCR.Rdata")

# dog ---------------------------------------------------------------------

rm(list = ls())

dog <- read.table(file = "01-data/26-Evolution_newly_emerging_UCR_distribution/dog_results.xls", 
                  sep = "\t")
colnames(dog) <- c("qseqid", "sseqid", "pident", "length", 
                   "mismatch", "gapopen", 
                   "qstart", "qend", "sstart", "send", 
                   "qseq", "sseq", "evalue", "bitscore")
dog <- dog[c(dog$sstart == 1 | dog$send == 1), ]   # 和下面一行一起组成：能匹配到ucr全长
dog <- dog[c(dog$length == (abs(dog$send - dog$sstart)) + 1), ]
dog <- dog[c(dog$pident == 100), ]
dog$sseqid <- gsub("::.*", "", dog$sseqid)
dog <- dog[!(c(dog$sseqid == "uc.18" | dog$sseqid == "uc.304")), ]

# 使用table函数统计每个数值的频次
freq_sseqid <- table(dog$sseqid)

# 筛选出出现超过两次的数值
result <- names(freq_sseqid[freq_sseqid > 1])
print(result)
# [1] ""

dog <- dog[, c(2, 9, 10, 12)]   # 现在不需要统计匹配到哪里，只要有100%匹配的位置就可以
dog <- unique(dog)
length(unique(dog$sseqid))
# 261 UCR

# dog <- dog[!(dog$sseqid == "uc.471" & dog$sstart == 205), ]

save(dog, file = "02-analysis/26-Evolution_newly_emerging_UCR_distribution/dog_UCR.Rdata")

# chicken -----------------------------------------------------------------

rm(list = ls())

chicken <- read.table(file = "01-data/26-Evolution_newly_emerging_UCR_distribution/chicken_results.xls", 
                  sep = "\t")
colnames(chicken) <- c("qseqid", "sseqid", "pident", "length", 
                   "mismatch", "gapopen", 
                   "qstart", "qend", "sstart", "send", 
                   "qseq", "sseq", "evalue", "bitscore")
chicken <- chicken[c(chicken$sstart == 1 | chicken$send == 1), ]   # 和下面一行一起组成：能匹配到ucr全长
chicken <- chicken[c(chicken$length == (abs(chicken$send - chicken$sstart)) + 1), ]
chicken <- chicken[c(chicken$pident == 100), ]
chicken$sseqid <- gsub("::.*", "", chicken$sseqid)
chicken <- chicken[!(c(chicken$sseqid == "uc.18" | chicken$sseqid == "uc.304")), ]

# 使用table函数统计每个数值的频次
freq_sseqid <- table(chicken$sseqid)

# 筛选出出现超过两次的数值
result <- names(freq_sseqid[freq_sseqid > 1])
print(result)
# [1] ""

chicken <- chicken[, c(2, 9, 10, 12)]   # 现在不需要统计匹配到哪里，只要有100%匹配的位置就可以
chicken <- unique(chicken)
length(unique(chicken$sseqid))
# 31 UCR

# chicken <- chicken[!(chicken$sseqid == "uc.471" & chicken$sstart == 205), ]

save(chicken, file = "02-analysis/26-Evolution_newly_emerging_UCR_distribution/chicken_UCR.Rdata")

# turkey ------------------------------------------------------------------

rm(list = ls())

turkey <- read.table(file = "01-data/26-Evolution_newly_emerging_UCR_distribution/turkey_results.xls", 
                  sep = "\t")
colnames(turkey) <- c("qseqid", "sseqid", "pident", "length", 
                   "mismatch", "gapopen", 
                   "qstart", "qend", "sstart", "send", 
                   "qseq", "sseq", "evalue", "bitscore")
turkey <- turkey[c(turkey$sstart == 1 | turkey$send == 1), ]   # 和下面一行一起组成：能匹配到ucr全长
turkey <- turkey[c(turkey$length == (abs(turkey$send - turkey$sstart)) + 1), ]
turkey <- turkey[c(turkey$pident == 100), ]
turkey$sseqid <- gsub("::.*", "", turkey$sseqid)
turkey <- turkey[!(c(turkey$sseqid == "uc.18" | turkey$sseqid == "uc.304")), ]

# 使用table函数统计每个数值的频次
freq_sseqid <- table(turkey$sseqid)

# 筛选出出现超过两次的数值
result <- names(freq_sseqid[freq_sseqid > 1])
print(result)
# [1] ""

turkey <- turkey[, c(2, 9, 10, 12)]   # 现在不需要统计匹配到哪里，只要有100%匹配的位置就可以
turkey <- unique(turkey)
length(unique(turkey$sseqid))
# 30 UCR

# turkey <- turkey[!(turkey$sseqid == "uc.471" & turkey$sstart == 205), ]

save(turkey, file = "02-analysis/26-Evolution_newly_emerging_UCR_distribution/turkey_UCR.Rdata")

# lizard ------------------------------------------------------------------

rm(list = ls())

lizard <- read.table(file = "01-data/26-Evolution_newly_emerging_UCR_distribution/lizard_results.xls", 
                  sep = "\t")
colnames(lizard) <- c("qseqid", "sseqid", "pident", "length", 
                   "mismatch", "gapopen", 
                   "qstart", "qend", "sstart", "send", 
                   "qseq", "sseq", "evalue", "bitscore")
lizard <- lizard[c(lizard$sstart == 1 | lizard$send == 1), ]   # 和下面一行一起组成：能匹配到ucr全长
lizard <- lizard[c(lizard$length == (abs(lizard$send - lizard$sstart)) + 1), ]
lizard <- lizard[c(lizard$pident == 100), ]
lizard$sseqid <- gsub("::.*", "", lizard$sseqid)
lizard <- lizard[!(c(lizard$sseqid == "uc.18" | lizard$sseqid == "uc.304")), ]

# 使用table函数统计每个数值的频次
freq_sseqid <- table(lizard$sseqid)

# 筛选出出现超过两次的数值
result <- names(freq_sseqid[freq_sseqid > 1])
print(result)
# [1] ""

lizard <- lizard[, c(2, 9, 10, 12)]   # 现在不需要统计匹配到哪里，只要有100%匹配的位置就可以
lizard <- unique(lizard)
length(unique(lizard$sseqid))
# 11 UCR

# lizard <- lizard[!(lizard$sseqid == "uc.471" & lizard$sstart == 205), ]

save(lizard, file = "02-analysis/26-Evolution_newly_emerging_UCR_distribution/lizard_UCR.Rdata")

# x.tropicalis ------------------------------------------------------------

rm(list = ls())

x.tropicalis <- read.table(file = "01-data/26-Evolution_newly_emerging_UCR_distribution/x.tropicalis_results.xls", 
                  sep = "\t")
colnames(x.tropicalis) <- c("qseqid", "sseqid", "pident", "length", 
                   "mismatch", "gapopen", 
                   "qstart", "qend", "sstart", "send", 
                   "qseq", "sseq", "evalue", "bitscore")
x.tropicalis <- x.tropicalis[c(x.tropicalis$sstart == 1 | x.tropicalis$send == 1), ]   # 和下面一行一起组成：能匹配到ucr全长
x.tropicalis <- x.tropicalis[c(x.tropicalis$length == (abs(x.tropicalis$send - x.tropicalis$sstart)) + 1), ]
x.tropicalis <- x.tropicalis[c(x.tropicalis$pident == 100), ]
x.tropicalis$sseqid <- gsub("::.*", "", x.tropicalis$sseqid)
x.tropicalis <- x.tropicalis[!(c(x.tropicalis$sseqid == "uc.18" | x.tropicalis$sseqid == "uc.304")), ]

# 使用table函数统计每个数值的频次
freq_sseqid <- table(x.tropicalis$sseqid)

# 筛选出出现超过两次的数值
result <- names(freq_sseqid[freq_sseqid > 1])
print(result)
# [1] ""

x.tropicalis <- x.tropicalis[, c(2, 9, 10, 12)]   # 现在不需要统计匹配到哪里，只要有100%匹配的位置就可以
x.tropicalis <- unique(x.tropicalis)
length(unique(x.tropicalis$sseqid))
# 0 UCR

# x.tropicalis <- x.tropicalis[!(x.tropicalis$sseqid == "uc.471" & x.tropicalis$sstart == 205), ]

# save(x.tropicalis, file = "02-analysis/26-Evolution_newly_emerging_UCR_distribution/x.tropicalis_UCR.Rdata")

# zebrafish ---------------------------------------------------------------

rm(list = ls())

zebrafish <- read.table(file = "01-data/26-Evolution_newly_emerging_UCR_distribution/zebrafish_results.xls", 
                  sep = "\t")
colnames(zebrafish) <- c("qseqid", "sseqid", "pident", "length", 
                   "mismatch", "gapopen", 
                   "qstart", "qend", "sstart", "send", 
                   "qseq", "sseq", "evalue", "bitscore")
zebrafish <- zebrafish[c(zebrafish$sstart == 1 | zebrafish$send == 1), ]   # 和下面一行一起组成：能匹配到ucr全长
zebrafish <- zebrafish[c(zebrafish$length == (abs(zebrafish$send - zebrafish$sstart)) + 1), ]
zebrafish <- zebrafish[c(zebrafish$pident == 100), ]
# 0 UCR

# fugu --------------------------------------------------------------------

rm(list = ls())

fugu <- read.table(file = "01-data/26-Evolution_newly_emerging_UCR_distribution/fugu_results.xls", 
                        sep = "\t")
colnames(fugu) <- c("qseqid", "sseqid", "pident", "length", 
                         "mismatch", "gapopen", 
                         "qstart", "qend", "sstart", "send", 
                         "qseq", "sseq", "evalue", "bitscore")
fugu <- fugu[c(fugu$sstart == 1 | fugu$send == 1), ]   # 和下面一行一起组成：能匹配到ucr全长
fugu <- fugu[c(fugu$length == (abs(fugu$send - fugu$sstart)) + 1), ]
fugu <- fugu[c(fugu$pident == 100), ]
# 0 UCR

# lungfish ----------------------------------------------------------------

rm(list = ls())

lungfish <- read.table(file = "01-data/26-Evolution_newly_emerging_UCR_distribution/lungfish_results.xls", 
                        sep = "\t")
colnames(lungfish) <- c("qseqid", "sseqid", "pident", "length", 
                         "mismatch", "gapopen", 
                         "qstart", "qend", "sstart", "send", 
                         "qseq", "sseq", "evalue", "bitscore")
lungfish <- lungfish[c(lungfish$sstart == 1 | lungfish$send == 1), ]   # 和下面一行一起组成：能匹配到ucr全长
lungfish <- lungfish[c(lungfish$length == (abs(lungfish$send - lungfish$sstart)) + 1), ]
lungfish <- lungfish[c(lungfish$pident == 100), ]
# 0 UCR

# symmary -----------------------------------------------------------------

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
