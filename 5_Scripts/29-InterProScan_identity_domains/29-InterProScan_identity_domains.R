# 识别UCR所在的外显子发挥什么功能

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(ensembldb)

# 使用ensDb()函数连接到ENSEMBL数据库，并选择相应的物种
ensembl_db <- EnsDb("human")  # 替换"species"为你所研究的物种，如"human"或"mouse"














