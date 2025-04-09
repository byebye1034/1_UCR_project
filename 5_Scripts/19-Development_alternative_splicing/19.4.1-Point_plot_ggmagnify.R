# 使用ggmagnify绘制散点图（除此之外，类似的方法还有ggforce包）

install.packages("ggmagnify", repos = c("https://hughjonesd.r-universe.dev", 
                                        "https://cloud.r-project.org"))

setwd(dir = "D:/R_project/UCR_project/")
rm(list = ls())

library(tidyverse)
library(ggmagnify)
library(ggview)
library(cowplot)

load("D:/R_project/UCR_project/02-analysis/19-Development_alternative_splicing/devCE_hexamers_summary_merged.Rdata")

lwd_pt <- .pt*72.27/96

range(hexamers_summary_merged$frequency.up)
# [1] 0.000000000 0.001567612
range(hexamers_summary_merged$frequency.down)
# [1] 0.000000000 0.002387757

# 绘制原始的点图
p1 <- ggplot(data = hexamers_summary_merged) +
  geom_point(mapping = aes(
    x = frequency.down, y = frequency.up
  ), 
  size = 0.15, position = "jitter") +
  
  # 添加从左下到右上的虚线
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  
  labs(
    x = "human brain down devCE", 
    y = "human brain up devCE", 
    title = ""
  ) +
  
  theme(
    panel.grid = element_blank(), 
    panel.background = element_blank(),
    panel.border = element_rect(linewidth = 0.5/lwd_pt, 
                                color = "black", 
                                fill = NA), 
    
    axis.text = element_text(size = 5, color = "black"), 
    axis.line.x = element_line(linewidth = 0.5/lwd_pt), 
    axis.line.y = element_line(linewidth = 0.5/lwd_pt), 
    axis.ticks = element_line(linewidth = 0.5/lwd_pt), 
    axis.title = element_text(size = 5, color = "black"), 
    
    aspect.ratio = 1:1, 
    plot.margin = unit(c(0,0,0,0), "cm")
  ) +
  
  scale_x_continuous(limits = c(0, 0.00250), 
                     expand = c(0.000125, 0), 
                     labels = scales::scientific) +
  scale_y_continuous(limits = c(0, 0.00250), 
                     expand = c(0.000125, 0), 
                     labels = scales::scientific)
p1

p2 <- p1 +
  geom_segment(mapping = aes(x = 0.0005, xend = 0.0012, y = 0.0005, yend = 0.0005), linewidth = 0.5/lwd_pt) +
  geom_segment(mapping = aes(x = 0.0005, xend = 0.0012, y = 0.0012, yend = 0.0012), linewidth = 0.5/lwd_pt) +
  geom_segment(mapping = aes(x = 0.0005, xend = 0.0005, y = 0.0005, yend = 0.0012), linewidth = 0.5/lwd_pt) +
  geom_segment(mapping = aes(x = 0.0012, xend = 0.0012, y = 0.0005, yend = 0.0012), linewidth = 0.5/lwd_pt)
p2

ggview(p2, width = 4.8, height = 4.8, units = "cm", dpi = 1200)

save(p2, file = "03-results/19-Development_alternative_splicing/p2.Rdata")
save(p2, file = "03-results/figures/fig3/the_motif_enrichment_of_srsf1_in_human_down_devce_p2.Rdata")

# 绘制局部散点图 -----------------------------------------------------------------

# 按照down进行排序
sorted_data <- hexamers_summary_merged %>%
  arrange(desc(frequency.down))

# 取出down排名前10的数据
down_top_10 <- head(sorted_data, 10)

# 按照up进行排序
sorted_data <- hexamers_summary_merged %>%
  arrange(desc(frequency.up))

# 取出up排名前10的数据
up_top_10 <- head(sorted_data, 10)

part_hexamers_summary_merged <- hexamers_summary_merged %>% 
  filter(frequency.up > 0.0005 & frequency.up < 0.00120) %>% 
  filter(frequency.down > 0.0005 & frequency.down < 0.00120)

part_hexamers_summary_merged_notop <- anti_join(part_hexamers_summary_merged, up_top_10)
part_hexamers_summary_merged_notop <- anti_join(part_hexamers_summary_merged_notop, down_top_10)

hexamers_point_part_p <- ggplot(data = part_hexamers_summary_merged_notop) +
  geom_point(mapping = aes(x = frequency.down, 
                           y = frequency.up), 
             size = 0.15, position = "jitter") +
  
  # 添加从左下到右上的虚线
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  
  labs(
    x = "human brain down devCE", 
    y = "human brain up devCE", 
    title = ""
  ) +

  scale_x_continuous(limits = c(0.0005, 0.0012), 
                     expand = c(0, 0), 
                     labels = scales::scientific) +
  scale_y_continuous(limits = c(0.0005, 0.0012), 
                     expand = c(0, 0), 
                     labels = scales::scientific) +
  
  theme(
    panel.grid = element_blank(), 
    panel.background = element_blank(),
    panel.border = element_rect(linewidth = 0.5/lwd_pt, 
                                color = "black", 
                                fill = NA), 
    
    axis.text = element_text(size = 5, color = "black"), 
    axis.line.x = element_line(linewidth = 0.5/lwd_pt), 
    axis.line.y = element_line(linewidth = 0.5/lwd_pt), 
    axis.ticks = element_line(linewidth = 0.5/lwd_pt), 
    axis.title = element_text(size = 5, color = "black"), 
    
    aspect.ratio = 1:1, 
    plot.margin = unit(c(0,0,0,0), "cm")
  )

hexamers_point_part_p

p3 <- hexamers_point_part_p +
  
  geom_point(data = up_top_10, 
             mapping = aes(
               x = frequency.down, 
               y = frequency.up
             ), 
             color = "#4DBBD5FF", 
             size = 0.35) +
  
  geom_text_repel(data = up_top_10, 
                  mapping = aes(
                    x = frequency.down, 
                    y = frequency.up, 
                    label = hexamers
                  ), 
                  hjust = -0.2, vjust = 0.5, 
                  size = 1, color = "#4DBBD5FF", 
                  segment.size = 0.5/lwd_pt, 
                  max.overlaps = 20) +
  
  geom_point(data = down_top_10, 
             mapping = aes(
               x = frequency.down, 
               y = frequency.up
             ), 
             color = "#E64B35FF", 
             size = 0.35) +
  
  geom_text_repel(data = down_top_10, 
                  mapping = aes(
                    x = frequency.down, 
                    y = frequency.up, 
                    label = hexamers
                  ), 
                  hjust = 1, vjust = 0.5, 
                  size = 1, color = "#E64B35FF", 
                  segment.size = 0.5/lwd_pt, 
                  max.overlaps = 20)
p3

ggview(plot = ggplot2::last_plot(), 
       width = 4.8, height = 4.8, units = "cm", dpi = 1200)

save(p3, file = "03-results/19-Development_alternative_splicing/p3.Rdata")
save(p3, file = "03-results/figures/fig3/the_motif_enrichment_of_srsf1_in_human_down_devce_p3.Rdata")

cowplot::plot_grid(p2, p3, 
                   align = c("hv"), axis = c("tb"), 
                   nrow = 1, ncol = 2, rel_widths = 1)

ggview(plot = ggplot2::last_plot(), 
       width = 10, height = 6, units = "cm", dpi = 1200)

# 保存成pdf ------------------------------------------------------------------

pdf(file = "03-results/19-Development_alternative_splicing/2.pdf", width = 1000/254, height = 600/254)

cowplot::plot_grid(p2, p3, 
                   align = c("hv"), axis = c("tb"), 
                   nrow = 1, ncol = 2, rel_widths = 1)

dev.off()




my_theme <- theme(
  panel.grid = element_blank(), 
  panel.background = element_blank(), 
  text = element_text(size = 5, color = "#000000"), 
  line = element_line(linewidth = 0.5/lwd_pt, color = "#000000"), 
  
  axis.title = element_text(size = 5, color = "#000000"), 
  axis.text = element_text(size = 5, color = "#000000"), 
  axis.ticks = element_line(linewidth = 0.5/lwd_pt, color = "#000000"), 
  axis.line = element_line(linewidth = 0.5/lwd_pt, color = "#000000"),
  
  legend.position = "top", 
  legend.key.size = unit(0.25, "cm"), 
  legend.title = element_text(size = 5, color = "#000000"), 
  legend.text = element_text(size = 5, color = "#000000"), 
  
  plot.title = element_text(size = 5, color = "#000000"),
  
  aspect.ratio = 1:1
)





