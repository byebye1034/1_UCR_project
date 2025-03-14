# figure中的配色
custom_palette <- c("#000000", "#E60212", "#083490", "#751384", "#007A34", "#D85F00", "#f1f1f1")
save(custom_palette, file = "D:/R_project/protocol/custom_palette.Rdata")

# "#000000" 转换为 RGB(0, 0, 0)
# "#E60212" 转换为 RGB(230, 2, 18)
# "#083490" 转换为 RGB(8, 52, 144)
# "#751384" 转换为 RGB(117, 19, 132)
# "#007A34" 转换为 RGB(0, 122, 52)
# "#D85F00" 转换为 RGB(216, 95, 0)
# "#f1f1f1" 转换为 RGB(241, 241, 241)

# pdf(path/file=, width=, heigh=)
# width 和 heigh 的单位是inches(英寸)，1 inch = 2.54 centimeter(cm)
# width = 400/254, 代表的是4 cm

# 单栏宽度：8cm, 800/254
# 2/3版图宽度：12-15cm, 1200/254 - 1500/254
# 1.5栏宽度：11cm, 1100/254
# 全版图宽度：17cm, 1700/254

# pdf输出时，图片比例为3：2时，高度设置为480/254，宽度720/254
# pdf输出时，图片比例为1：1时，高度设置为480/254，宽度480/254

# 文字大小：5，在AI当中打开也刚好是5 pt
# 线条粗细：width = 0.5/lwd_pt, (提前设置lwd_pt <- .pt*72.27/96)，在AI中打开为1
# 散点图中的点：7，放在aes的外部才能调节点的大小

# 对坐标轴参数的调整：文字，刻度，线条，标题

# ggplot2中柱形图柱子宽度默认是0.9，我设置为0.5

library(gginnards) # 更改显著性标记的字体大小
library(ggpubr) # stat_compare_means添加P值
library(gghalves) # 绘制云雨图

# PBMC_U937_DEG.R当中的例子
GSE173754 <- read_tsv(file = "02-analysis/08-PBMC_U937_DEG/GSE173754.top.table.tsv")
GSE173754 <- GSE173754 %>% 
  mutate(Type = as.factor(ifelse(log2FoldChange > 1 & padj < 0.05, "up", ifelse(log2FoldChange < -1 & padj < 0.05, "down", "no change"))))

p <- ggplot(data = GSE173754, mapping = aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Type), size = 0.7) +        # 点的大小0.7
  geom_boxplot(mapping = aes(color = Group), 
               linewidth = 0.5/lwd_pt,                                          # 设置箱线图的线条粗细
               outlier.size = 0.25) +                                           # 设置离群值点的大小
  
  scale_color_manual(values = c("#f1f1f1", "#083490", "#E60212")) +
  
  geom_hline(yintercept = c(-log10(0.05)),           # 添加水平方向的线
             linewidth = 0.5/lwd_pt, 
             lty = "dashed") +
  
  geom_vline(xintercept = c(-1, 1),                  # 添加垂直方向的线
             linewidth = 0.5/lwd_pt, 
             lty = "dashed") +
  
  guides(color = guide_legend(nrow = 2, 
                              ncol = 2)) +           # 设置图例为两行两列
  
  theme(
    panel.grid = element_blank(),                    # 背景网格删除
    panel.background = element_blank(),              # 背景删除
    text = element_text(size = 7),                   # 除坐标轴字体大小
    legend.position = "none",                        # 删除标签
    line = element_line(linewidth = 0.5/lwd_pt),     # 线条宽度（坐标轴上的点）
    
    axis.text = element_text(size = 7),              # 坐标轴字体大小
    axis.line = element_line(linewidth = 0.5/lwd_pt),  # 坐标轴线条粗细
	
	  aspect.ratio = 1,                                # 保持绘图区域为1：1
	
	  legend.position = "top",                         # 将图例放置在顶部
	  legned.key.size = unit(0.2, units = "cm")        # 设置图例大小为0.2cm
  )
p

p2$layers[[which_layers(p2, "GeomSignif")]]$aes_params$textsize <- 5/lwd_pt   # 更改显著性标记的字体大小

# 火山图 ---------------------------------------------------------------------

gene_of_interest <- c("PBX3", "BCL11A", "LCOR", "MECOM", "ZFHX3")
p + geom_text_repel(data = filter(GSE173754, Symbol %in% gene_of_interest), 
                    aes(label = Symbol), size = 4)

pdf(file = "03-results/figures/PBMC_U937_DEG.pdf", height = 480/254, width = 480/254)
p + geom_point(data = gene_of_interest, 
               mapping = aes(x = log2FoldChange, y = -log10(padj), fill = Symbol), 
               size = 0.7) +
  geom_text_repel(data = filter(GSE173754, Symbol %in% gene_of_interest$Symbol), 
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 50),
                  aes(label = Symbol), size = 3)
dev.off()

# 柱形圖 ---------------------------------------------------------------------

custom_palette <- c("#000000", "#E60212", "#083490", "#751384", "#007A34", "#D85F00")

lwd_pt <- .pt*72.27/96

ggplot(data = UCR_class, mapping = aes(x = UCR_type, y = UCR_number, color = UCR_type)) +
  geom_col(width = 0.8, fill = NA, size = 0.5/lwd_pt) +             # width:柱子宽度 size:边框的粗细
  
  scale_color_manual(values = custom_palette) +
  scale_y_continuous(limits = c(0, 250), expand = c(0, 0)) +
  
  theme_prism(palette = "black_and_white", 
              base_size = 7, 
              base_family = "serif", 
              base_fontface = "bold", 
              base_line_size = 0.5/lwd_pt, 
              base_rect_size = 0.5/lwd_pt, 
              axis_text_angle = 45, 
              border = FALSE) +
  
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(), 
    legend.key.size = unit(0.2, "cm"),                              # 调整图例大小：图例的key.size为0.2 cm
    legend.spacing.y = unit(1, "cm"),                               # 调整图例之间的垂直距离
    legend.position = c(1, 1),                                      # 调整图例位置到右上角
    legend.justification = c(1, 1)                                  # 设置图例的锚点在右上角
  )

# 常用参数
lwd_pt <- .pt*72.27/96

my_theme <- theme(
  panel.grid = element_blank(), 
  panel.background = element_blank(), 
  text = element_text(size = 7, color = "#000000"), 
  line = element_line(linewidth = 0.5/lwd_pt, color = "#000000"), 
  
  axis.title = element_text(size = 7, color = "#000000"), 
  axis.text = element_text(size = 7, color = "#000000"), 
  axis.ticks = element_line(linewidth = 0.5/lwd_pt, color = "#000000"), 
  axis.line = element_line(linewidth = 0.5/lwd_pt, color = "#000000"),
  
  legend.position = "top", 
  legend.key.size = unit(0.25, "cm"), 
  legend.title = element_text(size = 7, color = "#000000"), 
  legend.text = element_text(size = 7, color = "#000000"), 
  
  plot.title = element_text(size = 7, color = "#000000"),
  
  aspect.ratio = 1:1
)

library(sysfonts)
library(showtext)

font_add(family = "Arial", regular = "C:/Windows/Fonts/arial.ttf")
font_families()

pdf(file = "path_to_file/filename.pdf", width = 450/254, height = 600/254, dpi = 1200)
showtext_begin()

prism_theme <- theme_prism(base_size = 7, 
              base_family = "Arial", 
              base_fontface = "plain", 
              base_line_size = 0.5/lwd_pt, 
              base_rect_size = 0.5/lwd_pt, 
              border = TRUE) +
  
  theme(
    legend.position = "top", 
    legend.key.size = unit(0.25, "cm"), 
	aspect.ratio = 1
  )

showtext_end()
dev.off()

# ggsci的配色 ----------------------------------------------------------------

library(ggplot2)
library(ggsci)

# 获取 npg 配色方案的颜色
npg_palette <- ggsci::pal_npg()(10)

# 打印颜色的十六进制表示
print(npg_palette)

npg_palette <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF",
                 "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF")

# 获取jco（Journal of Clinical Oncology）配色方案的颜色
jco_palette <- ggsci::pal_jco()(10)

# 打印颜色的十六进制表示
print(jco_palette)

jco_palette <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#7AA6DCFF", 
                 "#003C67FF", "#8F7700FF", "#3B3B3BFF", "#A73030FF", "#4A6990FF")

# 基于ggplot2绘制图形的次刻度添加ggh4x ------------------------------------------------

# 基于ggplot2绘制一张常规图形
library(ggplot2)
p1<- ggplot(iris, aes(Species, Sepal.Length))+
  stat_boxplot(geom = "errorbar", width=0.1,linewidth=0.8)+
  geom_boxplot(aes(fill=Species), outlier.color = NA )+
  geom_jitter(aes(fill=Species), shape=21,size=2.5,alpha=0.6,width=0.3)+
  theme_bw()
p1

# 通过设置标签方式实现（缺陷：主刻度与次刻度长度一致）
p1 + scale_y_continuous(breaks = c(4,4.5,5,5.5,6,6.5,7,7.5,8),
                      labels = c(4,"","","",6,"","","",8),
                      limits = c(4,8))

# 基于ggh4x设置次刻度
library(ggh4x)
p1+scale_y_continuous(guide = "axis_minor",
                      minor_breaks = seq(4, 8, by = 0.2),#需要设置次刻度的范围及间隔
                      limits = c(4,8))+#坐标轴整体范围
  theme(axis.ticks.length = unit(1.5, "mm"),#主刻度长度设置
        ggh4x.axis.ticks.length.minor = rel(0.5)) #主刻度与次刻度比例

# 字体 ----------------------------------------------------------------------

#调用 windows 的字体格式，本示例分别以 Arial 和 Times New Roman 为例进行演示
windowsFonts(Arial = windowsFont('Arial'), TNR = windowsFont('Times New Roman'))

#将坐标轴和图例中的字体改为 Arial
p + theme(
  axis.title = element_text(family = 'Arial'), 
  axis.text = element_text(family = 'Arial'), 
  legend.title = element_text(family = 'Arial'), 
  legend.text = element_text(family = 'Arial'))

#将坐标轴和图例中的字体改为 Times New Roman
p + theme(
  axis.title = element_text(family = 'TNR'), 
  axis.text = element_text(family = 'TNR'), 
  legend.title = element_text(family = 'TNR'), 
  legend.text = element_text(family = 'TNR'))

# 使用sysfonts和showtext调用字体 -------------------------------------------------
# 如果只使用ggsave（）保存文件，可以使用前面调用字体的方式
# 如果需要输出为pdf文件，需要使用sysfonts和showtext包，方法如下

# 第一步，加载包
library(sysfonts)
library(showtext)

# 第二步，添加字体
font_add(family = "Arial", regular = "C:/Windows/Fonts/arial.ttf")
font_families() # 查看当前可用字体

# 第三步，保存pdf文件
# 打开图形设备
pdf(file = "03-results/37-PLEK/coding_pre_p.pdf", width = 480/254, height = 600/254)

# 使用showtext渲染字体
showtext_begin()

# 绘图
ggplot()

# 关闭showtext渲染字体
showtext_end()

# 关闭图形设备
dev.off()

# 2024-04-19 学习投稿级结果图的绘制 --------------------------------------------------

library(ggplot2)

# heat可以获得颜色比较接近的颜色
# RColorBrewer这个包里有一套颜色，可以使用display.brewer.all()展示

# geom_hisgram(col = "XX", fill = "XX") col:描边颜色 fill:填充颜色
# 如果要将变量group值映射为颜色，就要放在aes()当中：geom_hisgram(mapping = aes(col = factor(group))

# 图例的标题:Group
# labs(col = "Group")

# 调色函数的使用 RColorBrewer
library(RColorBrewer)

# 利用brewer.pal获得8种颜色
cols <- c(brewer.pal(4, "Set1"), brewer.pal(4, "Paired"))
cols_paired <- c(brewer.pal(8, "Paired"))
# 用8种颜色绘图
p1 <- ggplot(data = diamonds, mapping = aes(x = cut, fill = clarity)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = cols_paired)
p1

# 调整x轴刻度
# scale_x_continuous(breaks = seq(from = 3000, to = 6000, by = 500))

# 调整标题
# hjust调整水平方向[0, 1]，vjust调整垂直方向[0, 1]，方向是文字的水平和垂直
# y轴标题如果想上下调节，应该使用水平调节参数hjust
# title = element_text(family = “Arial", size = 7, face = "plain", hjust = 0.5) # 所有标题（x轴，y轴和图片上方的标题）
# plot.title = element_text() # 图片上方主标题
# axis.title = element_text() # 坐标轴标题
# axis.title.x = element_text() # x轴标题

# axis.text = element_text() # 刻度字体大小
# axis.ticks = element_line() # 刻度
# axis.line = element_line() # 线条

# 修改图例
# 移除图例
# 去除特定的图例：guides()，scale_fill_discrete()
# 去除全部的图例：theme(legend.position='none')
p2 <- p1 +
  theme(legend.position = "none") # 移除图例
p2
# 修改图例标题
p3 <- ggplot(data = diamonds, mapping = aes(x = cut, fill = clarity)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = cols_paired, 
                    name = "Clarity") # 因为在这里图例的标题就是指明颜色代表的是什么
p3
# 修改图例顺序
# 如果图例比较少的话，可以直接使用：breaks = c("normal", "cancer")
p4 <- ggplot(data = diamonds, mapping = aes(x = cut, fill = clarity)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = cols_paired, 
                    name = "Clarity", 
                    breaks = diamonds$clarity) # fill就改scale_fill_mannual,color就改scale_color_mannual,shape就改scale_shape_mannual
p4
# 重新定义标签名称，比如数据里是1和2，其实代表的是male和female(当然最好是在数据当中修改)
label <- diamonds$clarity
label <- gsub("VVS1", "\\[VVS1\\]", label) # 添加中括号【】，也可以用于去除【】
p5 <- ggplot(data = diamonds, mapping = aes(x = cut, fill = clarity)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = cols_paired, 
                    name = "Clarity", 
                    breaks = diamonds$clarity, 
                    labels = label)
p5

# 设置点的形状(直接指定使用参数shape，使用映射则放在aes（）里)
p6 <- ggplot(data = mtcars, mapping = aes(x = wt, y = mpg)) +
  geom_point(mapping = aes(col = as.factor(cyl), shape = as.factor(cyl))) +
  scale_color_manual(values = c("red", "green", "blue"), name = "Cyl") +
  scale_fill_manual(values = c("red", "green", "blue"))
p6
ggsave("03-results/figures/p6.pdf", width = 9.6, height = 9.6, units = "cm", dpi = 1200)

# 线的类型(水平线：geom_hline() 垂直线：geom_vline())
df <- data.frame(x = 1:25, y = 1:25)
p7 <- ggplot(data = df, mapping = aes(x = x, y = y)) +
  geom_abline(linetype = 2)
p7
# 坐标的线（theme()当中修改的是和数据无关的内容，如果和数据存在映射关系就放在aes()里）
# theme(
#   axis.line = element_line(size = 0.5/lwd_pt)
# )

# 添加标注信息
# geom_text(aes(label = data$label), size = xx, vjust = xx)：添加和数据有关的标注
# 如果只是单独标记某个点的话，就使用annotate()：“text” “rect”
p8 <- ggplot(data = df, mapping = aes(x = x, y = y)) +
  geom_point() +
  annotate("text", x = 10, y = 10, label = "point",  # 添加文本"text"
           hjust = 2, color = "red")
p8

p9 <- ggplot(data = df, mapping = aes(x = x, y = y)) +
  geom_point() +
  annotate("rect", xmin = 8, xmax = 10, ymin = 8, ymax = 10, # 添加图形"rect"
          color = "red", fill = NA)
p9

# 图片布局
# 分面：一行多列
p1 <- ggplot(data = mtcars, mapping = aes(x = mpg, y = wt)) +
  stat_smooth(method = lm) +
  geom_point() +
  facet_wrap(~ cyl) # 将cyl作为一个factor进行分面
p1
# 分面：一列多行
p2<- ggplot(data = mtcars, mapping = aes(x = mpg, y = wt)) +
  stat_smooth(method = lm) +
  geom_point() +
  facet_wrap(~ cyl, nrow = 3) # 将cyl作为一个factor进行分面
p2

# 坐标的调节
p3 <- ggplot(data = mtcars, mapping = aes(x = mpg, y = wt)) +
  stat_smooth(method = lm) +
  geom_point() +
  facet_wrap(~ cyl, scales = "free") # 将cyl作为一个factor进行分面
p3

p4 <- ggplot(data = mtcars, mapping = aes(x = mpg, y = wt)) +
  stat_smooth(method = lm) +
  geom_point() +
  facet_wrap(~ cyl, scales = "free_y") # x轴保持一致
p4

p5 <- ggplot(data = mtcars, mapping = aes(x = mpg, y = wt)) +
  stat_smooth(method = lm) +
  geom_point() +
  facet_wrap(~ cyl, scales = "free_x") # y轴保持一致
p5

# 一页多图：ggpubr

# 来自公众号的配色 ----------------------------------------------------------------

# 双色
two_color_set1 <- c("#E3738B", "#713948", 
                    "#F9D5DD", "#D7C6CA")

two_color_set2 <- c("#8CA5EA", "#495373", 
                    "#DCE4FA", "#C9CBD6")

two_color_set3 <- c("#FCB462", "#795B34", 
                    "#FFE8CE", "#DACEC2")

two_color_set4 <- c("#7BC4C5", "#3F6561", 
                    "#D9EEEE", "#C5D1D2")

two_color_set5 <- c("#C17F9E", "#8BACD1")

# 三色
three_color_set1 <- c("#A3A5A6", "#F8B496", "#B0B7DB")

three_color_set2 <- c("#F7D8B7", "#E7C1D7", "#AFB0D0", 
                      "#F0B47C", "#BE7A9A", "#7985B3")

# 多色
four_color_set1 <- c("#989A9C", "#F7D08D", "#BF83A5", "#8684B0")

five_color_set1 <- c("#9DB4CE", "#F9C08A", "#EDA1A4", "#B3D8D5", "#A4CB9E")

# ggside：边际组合图的简易语法 -------------------------------------------------------
# https://mp.weixin.qq.com/s/FDHNUGHw4-o5THl2ZsGJDg

# significance test
NS <- c(0.8567453, 1.3971413, 0.7461134)
LPS <- c(0.7182861, 0.5770795, 0.9473288, 0.3749752)
x <- c(NS, LPS)
group <- c(rep("NS", 3), rep("LPS", 4))

shapiro.test(NS)
shapiro.test(LPS)

bartlett.test(x~group) # 方差齐性检验，p接近1说明方差齐

t.test(NS, LPS, paired = FALSE, var.equal = T)

# pie plot
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



