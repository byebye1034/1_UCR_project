# 分析一下划痕实验和transwell实验的结果

setwd(dir = "D:/R_project/UCR_project/")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(tidyverse)
library(readxl)
library(ggplot2)
library(reshape2)
library(ggprism)
library(ggview)
library(ggsignif)

# caki1 -------------------------------------------------------------------

caki1_ctrl_0h <- read_table("D:/A_projects/UCR/Wound_assay/caki1/ctrl/0h/Results.xls")
caki1_ctrl_0h <- caki1_ctrl_0h %>% 
  mutate(group = "ctrl", time = "0h")

caki1_ctrl_12h <- read_table("D:/A_projects/UCR/Wound_assay/caki1/ctrl/12h/Results.xls")
caki1_ctrl_12h <- caki1_ctrl_12h %>% 
  mutate(group = "ctrl", time = "12h")

caki1_ctrl_24h <- read_table("D:/A_projects/UCR/Wound_assay/caki1/ctrl/24h/Results.xls")
caki1_ctrl_24h <- caki1_ctrl_24h %>% 
  mutate(group = "ctrl", time = "24h")

caki1_oe_0h <- read_table("D:/A_projects/UCR/Wound_assay/caki1/oe/0h/Results.xls")
caki1_oe_0h <- caki1_oe_0h %>% 
  mutate(group = "oe", time = "0h")

caki1_oe_12h <- read_table("D:/A_projects/UCR/Wound_assay/caki1/oe/12h/Results.xls")
caki1_oe_12h <- caki1_oe_12h %>% 
  mutate(group = "oe", time = "12h")

caki1_oe_24h <- read_table("D:/A_projects/UCR/Wound_assay/caki1/oe/24h/Results.xls")
caki1_oe_24h <- caki1_oe_24h %>% 
  mutate(group = "oe", time = "24h")

merge_results <- rbind(caki1_ctrl_0h, 
                       caki1_ctrl_12h, 
                       caki1_ctrl_24h, 
                       caki1_oe_0h, 
                       caki1_oe_12h, 
                       caki1_oe_24h)

merge_results <- arrange(merge_results, 
                         group, Area, time)

merge_results$Area <- factor(merge_results$Area, 
                             levels = c("1", "2", "3", "4"))
merge_results$time <- factor(merge_results$time, 
                             levels = c("0h", "12h", "24h"))

wound_results <- ggplot(data = merge_results, 
                        mapping = aes(x = group, y = Mean)) +
  
  geom_point(mapping = aes(color = time)) +
  
  facet_wrap(~ Area, ncol = 4) +
  
  theme_classic()

wound_results

closure_rate <- merge_results %>% 
  group_by(group, Area) %>% 
  mutate(closure_rate = 1 - (Mean / max(Mean))) %>% 
  ungroup()

h12_h24_closure_rate <- closure_rate %>% 
  filter(time %in% c("12h", "24h")) %>% 
  arrange(time, group)

ggplot(data = h12_h24_closure_rate) +
  
  geom_bar(mapping = aes(x = time, y = closure_rate, 
                               fill = group, color = group))

# transwell ---------------------------------------------------------------

lwd_pt <- .pt*72.27/96

# caki1
caki1_transwell <- tibble(
  ctrl = c(335, 362, 340), 
  oe = c(260, 259, 338)
)

caki1_transwell_l <- melt(data = caki1_transwell, 
                        variable.name = "group", 
                        value.name = "count")

t.test(caki1_transwell$ctrl, 
       caki1_transwell$oe)

p_transwell_caki1 <- ggplot(data = caki1_transwell_l) +
  
  geom_col(mapping = aes(x = group, y = count, fill = NA), 
           position = position_dodge(1), width = 0.7) +
  geom_jitter(mapping = aes(x = group, y = count), size = 0.25) +
  
  scale_y_continuous(limits = c(0, 400), expand = c(0, 0)) +
  
  theme_prism(base_size = 7, base_fontface = "plain", 
              base_line_size = 0.5/lwd_pt, base_rect_size = 0.5/lwd_pt) +
  
  theme(
    legend.position = "none"
  )

p_transwell_caki1

p_transwell_caki1 +
  
  geom_signif(y_position = c(350), xmin = c(0.5), xmax = c(1), 
              annotation = c("ns"), tip_length = 0, size = 0.5/lwd_pt, textsize = 5)

# make a p-value table
df_p_val <- data.frame(
  group1 = "ctrl",
  group2 = "oe",
  p.adj = 0.139,
  y.position = 380, 
  label.size = 7, 
  bracket.size = 0.5/lwd_pt
)

p_transwell_caki1 + add_pvalue(df_p_val)

ggview(ggplot2::last_plot(), width = 4.8, height = 6, units = "cm", dpi = 1200)

# 786-O
ctrl <- c(519, 512, 603)
oe <- c(366, 433, 358)
t.test(ctrl, oe)
















