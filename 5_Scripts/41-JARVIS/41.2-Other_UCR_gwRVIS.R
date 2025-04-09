# �鿴other UCR��other RF��gwRVIS��û����������

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(ggview)
library(ggpubr)
library(ggsci)
library(ggbeeswarm)
library(Hmisc)
library(car)

# other UCR gwRVIS --------------------------------------------------------

other_UCR_gwRVIS <- read.table(file = "02-analysis/41-JARVIS/other_UCR_gwRVIS.bed", sep = "\t")
other_UCR_gwRVIS <- other_UCR_gwRVIS %>% 
  mutate(start_distance = V2 - V6, 
         end_distance = V7 - V3) %>% 
  filter(start_distance > 0 & end_distance > 0) %>% 
  select(V4, V8) %>% 
  mutate(type = "other_UCR")

# other rf gwRVIS ---------------------------------------------------------

other_RF_gwRVIS <- read.table(file = "02-analysis/41-JARVIS/other_RF_gwRVIS.bed", sep = "\t")
other_RF_gwRVIS <- other_RF_gwRVIS %>% 
  mutate(start_distance = V2 - V6, 
         end_distance = V7 - V3) %>% 
  filter(start_distance > 0 & end_distance > 0) %>% 
  select(V4, V8) %>% 
  mutate(type = "other_RF")

wilcox.test(other_UCR_gwRVIS$V8, other_RF_gwRVIS$V8)
t.test(other_UCR_gwRVIS$V8, other_RF_gwRVIS$V8)

other_UCR_RF_gwRVIS <- rbind(other_UCR_gwRVIS, 
                             other_RF_gwRVIS)

# ��֤�����Ƿ�����̬�ֲ�
qqPlot(lm(V8~type, data = other_UCR_RF_gwRVIS), simulate = TRUE, main = "QQ plot", labels = FALSE)

# Shapiro-Wilk ���飬���ҽ������� p ֵ������ 0.05 ʱ�������ݷ�����̬�ֲ�
shapiro <- tapply(other_UCR_RF_gwRVIS$V8, other_UCR_RF_gwRVIS$type, shapiro.test)
shapiro$other_UCR$p.value
shapiro$other_RF$p.value
# ������̬�ֲ�
















