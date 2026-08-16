#!/usr/bin/env Rscript

# ================================================================
# 相关性分布直方图（堆叠，火山图风格分类，颜色修正）
# 用法：Rscript plot_correlation_histogram_stacked.R -i corr.txt -o output
# ================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
  library(ggplot2)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "输入文件（含 Correlation 和 FDR 列）"),
  make_option(c("-o", "--output"), type = "character", default = "Correlation_histogram",
              help = "输出前缀"),
  make_option(c("--fdr"), type = "numeric", default = 0.05,
              help = "FDR 阈值 [默认 0.05]"),
  make_option(c("--strength"), type = "numeric", default = 0.5,
              help = "相关性强度阈值 [默认 0.5]"),
  make_option(c("--width"), type = "numeric", default = 7,
              help = "图片宽度（英寸）"),
  make_option(c("--height"), type = "numeric", default = 5,
              help = "图片高度（英寸）")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input)) stop("必须提供 -i 参数")

# 读取数据
df <- read.delim(opt$input, stringsAsFactors = FALSE, check.names = FALSE)
df <- df[!is.na(df$Correlation), ]
df$FDR[is.na(df$FDR)] <- 1

# 分类（顺序固定，便于颜色映射）
df <- df %>%
  mutate(
    Category = factor(
      case_when(
        FDR < opt$fdr & Correlation > opt$strength  ~ "Significant Positive",
        FDR < opt$fdr & Correlation < -opt$strength ~ "Significant Negative",
        FDR < opt$fdr & abs(Correlation) <= opt$strength ~ "Significant Weak",
        TRUE ~ "Not significant"
      ),
      levels = c("Significant Positive", "Significant Negative", 
                 "Significant Weak", "Not significant")  # 固定顺序
    )
  )

# 统计
cat_counts <- table(df$Category)
cat("各类别数量：\n")
print(cat_counts)

# 颜色映射（按因子顺序）
color_map <- c(
  "Significant Positive" = "#E41A1C",   # 红
  "Significant Negative" = "#377EB8",   # 蓝
  "Significant Weak" = "#FDB462",       # 橙/黄
  "Not significant" = "grey70"          # 灰
)

# 绘图（堆叠直方图）
p <- ggplot(df, aes(x = Correlation, fill = Category)) +
  geom_histogram(binwidth = 0.02, color = "white", size = 0.15, position = "stack") +
  scale_fill_manual(
    values = color_map,
    labels = paste0(names(cat_counts), " (n = ", cat_counts, ")")
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", size = 0.5) +
  labs(
    title = "Spearman Correlation Distribution",
    x = "Correlation (rho)",
    y = "Number of transcripts",
    fill = "Category"
  ) +
  theme_classic(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.8, 0.8),
    legend.background = element_rect(fill = "white", color = "black", size = 0.3),
    legend.title = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(paste0(opt$output, ".pdf"), p, width = opt$width, height = opt$height)
cat("堆叠直方图已保存至: ", paste0(opt$output, ".pdf"), "\n")
