#!/usr/bin/env Rscript

# ==========================================================
# 对 PDUI 矩阵进行逐转录本（PAC）统计，并输出整体分布图
# 用法: Rscript PAC_statistics_summary.R -i <PDUI_matrix> -o <output_prefix> [--tol <tolerance>]
# ==========================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(reshape2)  # 用于画多个分布图
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "PDUI matrix file (tab-separated, rows=PACs, columns=samples, first column is PAC ID)"),
  make_option(c("-o", "--out"), type = "character", default = "PAC_stats", help = "Output file prefix [default %default]"),
  make_option(c("--tol"), type = "double", default = 1e-9, help = "Tolerance for considering a value as 0 or 1 (|x - 0| <= tol -> 0, |x - 1| <= tol -> 1) [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ==========================================================
# 1. 读取数据
# ==========================================================
cat("Reading PDUI matrix:", opt$input, "\n")
mat <- fread(opt$input, data.table = FALSE)
rownames(mat) <- mat[, 1]
mat <- mat[, -1]  # 去掉第一列
mat <- as.matrix(mat)
mode(mat) <- "numeric"
cat("Matrix dimensions:", nrow(mat), "PACs ×", ncol(mat), "samples\n")

# ==========================================================
# 2. 辅助函数：统计每行的指标
# ==========================================================
calc_row_stats <- function(x, tol = 1e-9) {
  # x 是一个数值向量（一个 PAC 在全部样本中的值）
  # 返回命名向量
  
  # 缺失统计
  na_count <- sum(is.na(x))
  na_ratio <- na_count / length(x)
  
  # 非缺失值
  x_non_na <- x[!is.na(x)]
  n_non_na <- length(x_non_na)
  
  if (n_non_na == 0) {
    # 如果全部缺失，则方差、IQR 等设为 NA
    variance <- NA
    iqr <- NA
    count_0 <- NA
    count_1 <- NA
    non_extreme_ratio <- NA
  } else {
    variance <- var(x_non_na)
    iqr <- IQR(x_non_na)
    
    # 基于容差判断 0 和 1
    count_0 <- sum(abs(x_non_na - 0) <= tol)
    count_1 <- sum(abs(x_non_na - 1) <= tol)
    
    # 非极端值比例：0 < x < 1（严格，不考虑容差，或者可以用 > tol & < 1-tol）
    # 为了与之前保持一致，使用严格 >0 & <1，但考虑到浮点误差，这里也用容差
    non_extreme <- sum(x_non_na > tol & x_non_na < 1 - tol)
    non_extreme_ratio <- non_extreme / n_non_na
  }
  
  c(NA_count = na_count, NA_ratio = na_ratio,
    Variance = variance, IQR = iqr,
    Count_0 = count_0, Count_1 = count_1,
    NonExtreme_ratio = non_extreme_ratio,
    NonNA_count = n_non_na)  # 同时保留非缺失样本数，便于参考
}

# ==========================================================
# 3. 逐行计算（可能耗时，但可接受）
# ==========================================================
cat("Calculating per-PAC statistics...\n")
stats_list <- apply(mat, 1, calc_row_stats, tol = opt$tol)
# 转置为数据框
stats_df <- as.data.frame(t(stats_list))
# 添加 PAC ID 列
stats_df$PAC <- rownames(mat)
# 重新排列列顺序
stats_df <- stats_df[, c("PAC", "NonNA_count", "NA_count", "NA_ratio", 
                         "Variance", "IQR", "Count_0", "Count_1", "NonExtreme_ratio")]

# 转换为数值（apply 可能产生字符，再转换）
stats_df[, -1] <- lapply(stats_df[, -1], as.numeric)

cat("Finished calculation.\n")

# ==========================================================
# 4. 输出每个 PAC 的统计表格
# ==========================================================
out_table <- paste0(opt$out, "_per_PAC.txt")
write.table(stats_df, file = out_table, row.names = FALSE, quote = FALSE, sep = "\t")
cat("Per-PAC statistics saved to:", out_table, "\n")

# ==========================================================
# 5. 汇总统计（各指标的总体分布）
# ==========================================================
# 对关键指标生成汇总表（如中位数、均值、分位数等）
summary_cols <- c("Variance", "IQR", "NonExtreme_ratio", "NA_ratio", "Count_0", "Count_1")
summary_stats <- data.frame(
  Metric = character(),
  Min = numeric(),
  Q1 = numeric(),
  Median = numeric(),
  Mean = numeric(),
  Q3 = numeric(),
  Max = numeric(),
  stringsAsFactors = FALSE
)

for (col in summary_cols) {
  vals <- stats_df[[col]]
  # 去掉 NA 值（例如全缺失的行会有 NA）
  vals <- vals[!is.na(vals)]
  if (length(vals) > 0) {
    q <- quantile(vals, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
    summary_stats <- rbind(summary_stats,
                           data.frame(Metric = col,
                                      Min = q[1], Q1 = q[2], Median = q[3],
                                      Mean = mean(vals, na.rm = TRUE),
                                      Q3 = q[4], Max = q[5]))
  }
}
cat("\n=== Summary statistics of per-PAC metrics ===\n")
print(summary_stats, digits = 3)

# 保存汇总表
write.table(summary_stats, file = paste0(opt$out, "_summary_stats.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

# ==========================================================
# 6. 绘制这些指标的分布图（整体分布）
# ==========================================================
cat("Generating distribution plots...\n")

# 我们想画几个关键指标的直方图：Variance, IQR, NonExtreme_ratio, NA_ratio
# 注意：很多值可能是0或极小，可以加上对数变换，但此处先用原始值，并允许用户后续调整

# 定义绘图函数
plot_dist <- function(data, x, binwidth = NULL, fill = "steelblue", title = NULL, xlab = NULL, log_scale = FALSE) {
  p <- ggplot(data, aes_string(x = x)) +
    geom_histogram(bins = 100, fill = fill, color = "black") +
    theme_bw() +
    labs(title = title, x = xlab, y = "Number of PACs")
  if (log_scale) {
    p <- p + scale_x_log10()
  }
  return(p)
}

# 准备数据：只包含有非缺失值的行（去除全缺失的）
plot_data <- stats_df[!is.na(stats_df$Variance), ]

# 1. 方差分布（加 log10 尺度，因为方差通常右偏）
p_var <- plot_dist(plot_data, "Variance", title = "Variance distribution per PAC",
                   xlab = "Variance of PDUI (log10 scale)", log_scale = TRUE)
ggsave(paste0(opt$out, "_variance_dist.png"), p_var, width = 6, height = 5)

# 2. IQR 分布
p_iqr <- plot_dist(plot_data, "IQR", title = "IQR distribution per PAC",
                   xlab = "Interquartile range (IQR)")
ggsave(paste0(opt$out, "_IQR_dist.png"), p_iqr, width = 6, height = 5)

# 3. 非极端值比例分布（0~1之间）
p_nonext <- plot_dist(plot_data, "NonExtreme_ratio", title = "Non-extreme ratio distribution (0<PDUI<1)",
                      xlab = "Proportion of non-extreme values per PAC")
ggsave(paste0(opt$out, "_nonExtreme_ratio_dist.png"), p_nonext, width = 6, height = 5)

# 4. NA 比例分布
p_na <- plot_dist(plot_data, "NA_ratio", title = "NA ratio distribution per PAC",
                  xlab = "Proportion of missing values")
ggsave(paste0(opt$out, "_NA_ratio_dist.png"), p_na, width = 6, height = 5)

# 5. 0 的数目分布（因为很多 PAC 可能有大量 0，可显示）
p0 <- plot_dist(plot_data, "Count_0", title = "Count of zeros per PAC",
                xlab = "Number of samples with PDUI = 0 (within tolerance)")
ggsave(paste0(opt$out, "_count0_dist.png"), p0, width = 6, height = 5)

# 6. 1 的数目分布
p1 <- plot_dist(plot_data, "Count_1", title = "Count of ones per PAC",
                xlab = "Number of samples with PDUI = 1 (within tolerance)")
ggsave(paste0(opt$out, "_count1_dist.png"), p1, width = 6, height = 5)

cat("All distribution plots saved with prefix:", opt$out, "\n")

# 可选：如果只想看某些指标，可自行注释。

# ==========================================================
# 7. 可选的散点图（Variance vs IQR）等，探索关系
# ==========================================================
p_scatter <- ggplot(plot_data, aes(x = IQR, y = Variance)) +
  geom_point(alpha = 0.1, size = 0.5) +
  scale_x_log10() + scale_y_log10() +
  theme_bw() +
  labs(title = "IQR vs Variance (log-log)", x = "IQR", y = "Variance")
ggsave(paste0(opt$out, "_IQR_vs_Variance.png"), p_scatter, width = 6, height = 5)

cat("Done!\n")
