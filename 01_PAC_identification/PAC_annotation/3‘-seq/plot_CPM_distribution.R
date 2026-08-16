#!/usr/bin/env Rscript

# 用法: Rscript plot_CPM_distribution.R TMM_CPM_matrix.txt [output_prefix]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1) {
  stop("Usage: Rscript plot_CPM_distribution.R TMM_CPM_matrix.txt [output_prefix]")
}

input_file <- args[1]
output_prefix <- ifelse(length(args) > 1, args[2], "CPM_distribution")

# 读取矩阵（行=PAC，列=样本）
cat("Reading matrix:", input_file, "\n")
mat <- fread(input_file, header = TRUE)
rownames(mat) <- mat[[1]]
mat <- as.matrix(mat[, -1])  # 去掉第一列（PAC ID）

# 将所有值展平为一维向量
all_values <- as.vector(mat)
cat("Total values:", length(all_values), "\n")
cat("Zero values:", sum(all_values == 0), " (", 
    round(100 * sum(all_values == 0) / length(all_values), 1), "%)\n")
cat("Positive values:", sum(all_values > 0), "\n")

# 去掉 NA（如果有）
all_values <- all_values[!is.na(all_values)]

# ================================
# 图1：原始 CPM 直方图（含 0 值，使用 log10 x 轴）
# ================================
df <- data.frame(CPM = all_values)

p1 <- ggplot(df, aes(x = CPM)) +
  geom_histogram(bins = 100, fill = "#2b83ba", color = "white", alpha = 0.8) +
  scale_x_log10(
    breaks = c(0, 0.1, 1, 10, 100, 1000),
    labels = c("0", "0.1", "1", "10", "100", "1000")
  ) +
  labs(
    title = "TMM CPM distribution (log10 x-axis)",
    x = "CPM (log10 scale)",
    y = "Count (PAC × sample)"
  ) +
  theme_bw()

ggsave(paste0(output_prefix, "_raw.png"), p1, width = 8, height = 5)
ggsave(paste0(output_prefix, "_raw.pdf"), p1, width = 8, height = 5)

# ================================
# 图2：log2(CPM+1) 直方图（更接近正态，便于观察）
# ================================
df$log2CPM <- log2(df$CPM + 1)

p2 <- ggplot(df, aes(x = log2CPM)) +
  geom_histogram(bins = 100, fill = "#d95f02", color = "white", alpha = 0.8) +
  labs(
    title = "log2(CPM+1) distribution",
    x = "log2(CPM + 1)",
    y = "Count (PAC × sample)"
  ) +
  theme_bw()

ggsave(paste0(output_prefix, "_log2.png"), p2, width = 8, height = 5)
ggsave(paste0(output_prefix, "_log2.pdf"), p2, width = 8, height = 5)

# ================================
# 图3：非零 CPM 的分布（排除 0 值，只看表达部分）
# ================================
df_pos <- df[df$CPM > 0, ]

if(nrow(df_pos) > 0) {
  p3 <- ggplot(df_pos, aes(x = log10(CPM))) +
    geom_histogram(bins = 80, fill = "#238b45", color = "white", alpha = 0.8) +
    labs(
      title = "Positive CPM distribution (log10, excluding zeros)",
      x = "log10(CPM)",
      y = "Count (PAC × sample)"
    ) +
    theme_bw()
  ggsave(paste0(output_prefix, "_positive.png"), p3, width = 8, height = 5)
  ggsave(paste0(output_prefix, "_positive.pdf"), p3, width = 8, height = 5)
}

# 输出统计摘要
cat("\n=== Summary statistics ===\n")
cat("Min (non-zero):", min(df$CPM[df$CPM > 0]), "\n")
cat("1st quartile (non-zero):", quantile(df$CPM[df$CPM > 0], 0.25), "\n")
cat("Median (all):", median(df$CPM), "\n")
cat("Mean (all):", mean(df$CPM), "\n")
cat("Max:", max(df$CPM), "\n")

cat("\nDone. Output files:\n")
cat(paste0(output_prefix, "_raw.{png,pdf}\n"))
cat(paste0(output_prefix, "_log2.{png,pdf}\n"))
if(nrow(df_pos) > 0) {
  cat(paste0(output_prefix, "_positive.{png,pdf}\n"))
}
