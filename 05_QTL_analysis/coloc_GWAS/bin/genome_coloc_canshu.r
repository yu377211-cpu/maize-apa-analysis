library(ggplot2)
library(gridExtra)
library(dplyr)
library(RColorBrewer)

# 阻止自动生成Rplots.pdf
if(!interactive()) pdf(NULL)

# 解析命令行参数
args <- commandArgs(trailingOnly = TRUE)

# 检查参数数量
if(length(args) != 3) {
  stop("Usage: Rscript genome_coloc_meihua.r <match_results_file> <chr_length_file> <output_prefix>", call. = FALSE)
}

# 获取参数
match_results_file <- args[1]  #500kb_match_results.txt
chr_length_file <- args[2]  #chrNameLength_sub.txt
output_prefix <- args[3]

# 读取数据
data <- read.table(match_results_file, header=TRUE, sep="\t")
chr_lengths <- read.table(chr_length_file, header=TRUE)

# 预处理数据
data <- data %>%
  mutate(
    chr = factor(chr.x, levels = paste0("Chr", 1:10)),
    logpval.x = -log10(pvalue.x),
    logpval.y = -log10(pvalue.y)
  )

# 计算染色体累计位置（用于x轴）
chr_lengths <- chr_lengths %>% 
  mutate(
    chr = factor(Chr, levels = paste0("Chr",1:10)),
    cumlen = cumsum(Length) - Length,
    midpos = cumlen + Length/2
  ) %>%
  arrange(chr)

# 合并染色体位置信息
plot_data <- data %>%
  left_join(chr_lengths, by = "chr") %>%
  mutate(
    global_pos.x = pos.x + cumlen,
    global_pos.y = pos.y + cumlen
  )

# 1. 数据准备 ----------------------------------------------------------------
# 假设plot_data已经准备好
plot_data <- plot_data %>%
  mutate(
    chr_num = as.numeric(gsub("Chr", "", Chr)),  # 提取染色体编号用于颜色分组
    color_group = ifelse(chr_num %% 2 == 1, "odd", "even")  # 奇偶染色体分组
  )

# 2. 颜色设置 ----------------------------------------------------------------
palette <- brewer.pal(8, "Set1")  # 使用RColorBrewer调色板
color_mapping <- c(
  "odd" = alpha(palette[1], 0.6),  # 奇数染色体颜色（带透明度）
  "even" = alpha(palette[2], 0.6)   # 偶数染色体颜色
)

# 3. 绘制.x位点图（cis-3aQTL）----------------------------------------------
p1 <- ggplot(plot_data, aes(x = global_pos.x, y = logpval.x)) +
  geom_point(aes(color = color_group), size = 1.5, alpha = 0.7) +
  scale_color_manual(values = color_mapping) +
  geom_hline(
    yintercept = -log10(2.17e-14), #修改
    linetype = "dashed", 
    color = "grey30",
    linewidth = 0.5
  ) +
  scale_x_continuous(
    name = "Chromosome",
    breaks = chr_lengths$midpos,
    labels = 1:10,
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    name = expression(-log[10](italic(p))), 
    #limits = c(0, max(plot_data$logpval.x, plot_data$logpval.y) * 1.05),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(title = "Seedling trans-3aQTL") + #修改
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  # 添加染色体分隔线
  geom_vline(
    xintercept = c(0, chr_lengths$cumlen[-1]), 
    color = "grey80", 
    linewidth = 0.3
  )

# 4. 绘制.y位点图（Zeamap GWAS）--------------------------------------------
p2 <- ggplot(plot_data, aes(x = global_pos.y, y = logpval.y)) +
  geom_point(aes(color = color_group), size = 1.5, alpha = 0.7) +
  scale_color_manual(values = color_mapping) +
  geom_hline(
    yintercept = -log10(6.52e-7), #修改 
    linetype = "dashed", 
    color = "grey30",
    linewidth = 0.5
  ) +
  scale_x_continuous(
    name = "Chromosome",
    breaks = chr_lengths$midpos,
    labels = 1:10,
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    name = expression(-log[10](italic(p))), 
    #limits = c(0, max(plot_data$logpval.x, plot_data$logpval.y) * 1.05), #强制y轴范围
    expand = expansion(mult = c(0, 0.05)) # 自动适配y轴范围
  ) +
  labs(title = "GWAS LeafAngle") + #修改
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  # 添加染色体分隔线
  geom_vline(
    xintercept = c(0, chr_lengths$cumlen[-1]), 
    color = "grey80", 
    linewidth = 0.3
  )

# 5. 合并图形并保存 --------------------------------------------------------
combined_plot <- grid.arrange(
  p1, p2, 
  ncol = 1,
  heights = c(1, 1.1)  # 下边稍大以容纳x轴标签
)

# 生成输出文件名
png_file <- paste0(output_prefix, ".png")
pdf_file <- paste0(output_prefix, ".pdf")

ggsave(
  png_file,
  combined_plot,
  width = 10, 
  height = 7,
  dpi = 300,
  bg = "white"
)
ggsave(
  pdf_file,
  combined_plot,
  width = 10,
  height = 7,
  device = cairo_pdf  # 获得更好的PDF渲染质量
)

# 彻底关闭所有图形设备
invisible(dev.off())
