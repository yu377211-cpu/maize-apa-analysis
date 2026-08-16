#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("用法: Rscript plot_na_distribution.R <输入文件> <输出图片> [图片标题]")
}

input_file  <- args[1]
output_file <- args[2]
plot_title  <- ifelse(length(args) >= 3, args[3], "转录本缺失率分布")

# ---------- 读取数据 ----------
# 假设文件为两列，无表头，用制表符或空格分隔
df <- read.table(input_file, header = FALSE, col.names = c("NA_count", "transcript_count"))

# 检查数据有效性
if (nrow(df) == 0) stop("输入文件为空")
max_na <- max(df$NA_count)
if (max_na == 0) stop("最大 NA 数为 0，无法划分百分比区间")

# ---------- 按每 10% 划分区间 ----------
breaks <- seq(0, max_na, length.out = 11)  # 产生 0%, 10%, ..., 100% 的切点
labels <- paste0(seq(0, 90, by = 10), "-", seq(10, 100, by = 10), "%")

df$bin <- cut(df$NA_count,
              breaks = breaks,
              include.lowest = TRUE,
              labels = labels)

# 汇总每个区间内的转录本总数
result <- aggregate(transcript_count ~ bin, data = df, sum)

# ---------- 绘制条形图 ----------
p <- ggplot(result, aes(x = bin, y = transcript_count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "缺失率区间", y = "转录本数量", title = plot_title) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---------- 保存图片 ----------
# ggsave 会根据文件扩展名自动选择格式（png, pdf, jpeg 等）
ggsave(output_file, p, width = 10, height = 6, dpi = 300)

# 保存统计结果到文本文件
stat_file <- sub("\\.[^.]+$", ".txt", output_file)   # 替换扩展名为 .txt
write.table(result, file = stat_file, row.names = FALSE, quote = FALSE, sep = "\t")

cat("✅ 图片已保存至:", output_file, "\n")
