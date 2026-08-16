#!/usr/bin/env Rscript

# ================================================================
# 绘制 Spearman 相关系数密度图，标注高/低阈值计数
# 用法：Rscript plot_correlation_density.R -i correlation.txt -o output_prefix
# ================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
  library(ggplot2)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "输入文件（含 Correlation 列）"),
  make_option(c("-o", "--output"), type = "character", default = "Correlation_density",
              help = "输出文件前缀"),
  make_option(c("--upper"), type = "numeric", default = 0.5,
              help = "正阈值 [默认 0.5]"),
  make_option(c("--lower"), type = "numeric", default = -0.5,
              help = "负阈值 [默认 -0.5]"),
  make_option(c("--width"), type = "numeric", default = 7,
              help = "图片宽度（英寸）"),
  make_option(c("--height"), type = "numeric", default = 5,
              help = "图片高度（英寸）")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input)) stop("必须提供 -i 参数")

# -------- 读取数据 --------
df <- read.delim(opt$input, stringsAsFactors = FALSE, check.names = FALSE)

# 检查是否有 Correlation 列
if (!"Correlation" %in% colnames(df)) {
  stop("输入文件必须包含 'Correlation' 列")
}

# 提取相关系数列
cor_vec <- df$Correlation
cor_vec <- cor_vec[!is.na(cor_vec)]  # 去除 NA

cat("总转录本数（非 NA 相关性）: ", length(cor_vec), "\n")

# 计数高于/低于阈值
upper_count <- sum(cor_vec > opt$upper, na.rm = TRUE)
lower_count <- sum(cor_vec < opt$lower, na.rm = TRUE)
cat("rho > ", opt$upper, ": ", upper_count, "\n", sep = "")
cat("rho < ", opt$lower, ": ", lower_count, "\n", sep = "")

# -------- 绘图 --------
p <- ggplot(data.frame(x = cor_vec), aes(x = x)) +
  geom_density(fill = "steelblue", alpha = 0.6, color = "black", size = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", size = 0.6) +
  geom_vline(xintercept = opt$upper, linetype = "dotted", color = "red", size = 0.8) +
  geom_vline(xintercept = opt$lower, linetype = "dotted", color = "blue", size = 0.8) +
  annotate("text", x = opt$upper + 0.05, y = Inf,
           label = paste0("rho > ", opt$upper, "\n", upper_count),
           hjust = 0, vjust = 1.5, color = "red", size = 4) +
  annotate("text", x = opt$lower - 0.05, y = Inf,
           label = paste0("rho < ", opt$lower, "\n", lower_count),
           hjust = 1, vjust = 1.5, color = "blue", size = 4) +
  labs(
    title = "Distribution of Spearman Correlation Coefficients",
    x = "Correlation (rho)",
    y = "Density"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

# 保存 PDF
out_pdf <- paste0(opt$output, ".pdf")
ggsave(out_pdf, p, width = opt$width, height = opt$height)
cat("密度图已保存至: ", out_pdf, "\n")

# 同时保存一个 PNG 副本（可选）
# ggsave(paste0(opt$output, ".png"), p, width = opt$width, height = opt$height, dpi = 300)

cat("完成！\n")
