#!/usr/bin/env Rscript
# ================================================================
# 比较两种分组方式对 APA 特异性的贡献（如组织 vs 亚群）
# 输出合并表格、二维散点图（Hex bin）、ΔTau 直方图及统计结果
# 用法：Rscript compare_tissue_vs_population.R -i PDUI.txt -g group1.txt -s group2.txt -o out
# ================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(ggplot2)
  library(reshape2)
})

# ------------------------- 参数解析 -------------------------
option_list <- list(
  make_option(c("-i", "--input"), type="character", help="PDUI 矩阵文件（tab分隔，第一列为转录本ID）"),
  make_option(c("-g", "--group1"), type="character", help="第一分组文件（有标题，第一列样本ID，第二列分组名）"),
  make_option(c("-s", "--group2"), type="character", help="第二分组文件（有标题，第一列样本ID，第二列分组名）"),
  make_option(c("-o", "--output"), type="character", default="Comparison", help="输出前缀")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$group1) || is.null(opt$group2)) {
  stop("必须提供 -i, -g, -s 三个参数", call. = FALSE)
}

# ------------------------- 1. 读取数据 -------------------------
cat("读取 PDUI 矩阵...\n")
pdui <- fread(opt$input, header=TRUE, sep="\t", data.table=FALSE)
rownames(pdui) <- pdui[,1]
pdui <- as.matrix(pdui[,-1])

# 读取分组文件（自动提取第二列名作为分组类型）
cat("读取第一分组文件...\n")
g1 <- fread(opt$group1, header=TRUE, sep="\t", data.table=FALSE)
group1_name <- colnames(g1)[2]          # 取第二列标题
g1_vec <- setNames(g1[,2], g1[,1])

cat("读取第二分组文件...\n")
g2 <- fread(opt$group2, header=TRUE, sep="\t", data.table=FALSE)
group2_name <- colnames(g2)[2]
g2_vec <- setNames(g2[,2], g2[,1])

# 取交集样本
common <- intersect(colnames(pdui), intersect(names(g1_vec), names(g2_vec)))
if (length(common) == 0) stop("样本ID完全不匹配")
pdui <- pdui[, common, drop=FALSE]
g1_vec <- g1_vec[common]
g2_vec <- g2_vec[common]

g1_factor <- factor(g1_vec)
g2_factor <- factor(g2_vec)
cat(sprintf("样本数: %d, %s数: %d, %s数: %d\n", 
            ncol(pdui), group1_name, nlevels(g1_factor), group2_name, nlevels(g2_factor)))

# ------------------------- 2. 定义辅助函数 -------------------------
calc_tau <- function(x) {
  # x 为数值向量，含 NA
  if (sum(!is.na(x)) < 2) return(NA)
  x_max <- max(x, na.rm=TRUE)
  if (is.na(x_max) || x_max == 0) return(0)
  x_norm <- x / x_max
  sum(1 - x_norm, na.rm=TRUE) / (length(x) - 1)
}

get_dominant <- function(x) {
  # 返回最大值对应的分组名（若并列，用分号连接）
  if (all(is.na(x))) return(NA)
  m <- max(x, na.rm=TRUE)
  idx <- which(x == m)
  paste(names(x)[idx], collapse=";")
}

# ------------------------- 3. 计算两组均值矩阵 -------------------------
cat(sprintf("计算 %s 特异性 ...\n", group1_name))
levels1 <- levels(g1_factor)
agg1 <- sapply(levels1, function(lev) {
  cols <- which(g1_factor == lev)
  if (length(cols) == 1) pdui[, cols] else rowMeans(pdui[, cols, drop=FALSE], na.rm=TRUE)
})
# 保持行名一致
agg1 <- agg1[!apply(agg1, 1, function(x) all(is.na(x))), , drop=FALSE]  # 全NA行在计算Tau时会返回NA，但保留
tau1 <- apply(agg1, 1, calc_tau)
dom1 <- apply(agg1, 1, get_dominant)
mean1 <- rowMeans(agg1, na.rm=TRUE)   # 平均PDUI（所有组）
sd1 <- apply(agg1, 1, sd, na.rm=TRUE)

cat(sprintf("计算 %s 特异性 ...\n", group2_name))
levels2 <- levels(g2_factor)
agg2 <- sapply(levels2, function(lev) {
  cols <- which(g2_factor == lev)
  if (length(cols) == 1) pdui[, cols] else rowMeans(pdui[, cols, drop=FALSE], na.rm=TRUE)
})
agg2 <- agg2[!apply(agg2, 1, function(x) all(is.na(x))), , drop=FALSE]
tau2 <- apply(agg2, 1, calc_tau)
dom2 <- apply(agg2, 1, get_dominant)
mean2 <- rowMeans(agg2, na.rm=TRUE)
sd2 <- apply(agg2, 1, sd, na.rm=TRUE)

# ------------------------- 4. 合并结果 -------------------------
all_genes <- union(rownames(agg1), rownames(agg2))
result_df <- data.frame(
  Transcript = all_genes,
  Mean_PDUI_1 = mean1[match(all_genes, rownames(agg1))],
  SD_1 = sd1[match(all_genes, rownames(agg1))],
  Tau_1 = tau1[match(all_genes, rownames(agg1))],
  Dominant_1 = dom1[match(all_genes, rownames(agg1))],
  Mean_PDUI_2 = mean2[match(all_genes, rownames(agg2))],
  SD_2 = sd2[match(all_genes, rownames(agg2))],
  Tau_2 = tau2[match(all_genes, rownames(agg2))],
  Dominant_2 = dom2[match(all_genes, rownames(agg2))],
  stringsAsFactors = FALSE
)
# 计算 DeltaTau (Tau_1 - Tau_2)
result_df$DeltaTau <- result_df$Tau_1 - result_df$Tau_2

# 重命名列（使用实际分组名称）
colnames(result_df)[grep("_1$", colnames(result_df))] <- 
  gsub("_1$", paste0("_", group1_name), colnames(result_df)[grep("_1$", colnames(result_df))])
colnames(result_df)[grep("_2$", colnames(result_df))] <- 
  gsub("_2$", paste0("_", group2_name), colnames(result_df)[grep("_2$", colnames(result_df))])

# 保存表格
output_table <- paste0(opt$output, "_Tau_Comparison.txt")
write.table(result_df, file=output_table, sep="\t", row.names=FALSE, quote=FALSE)
cat("合并结果已保存至:", output_table, "\n")

# ------------------------- 5. 统计和绘图准备 -------------------------
# 去除两列Tau均为NA的行（无法绘图）
plot_df <- result_df[!is.na(result_df[[paste0("Tau_", group1_name)]]) & 
                     !is.na(result_df[[paste0("Tau_", group2_name)]]), ]

epsilon <- 1e-6
above <- sum(plot_df[[paste0("Tau_", group1_name)]] > plot_df[[paste0("Tau_", group2_name)]] + epsilon)
below <- sum(plot_df[[paste0("Tau_", group2_name)]] > plot_df[[paste0("Tau_", group1_name)]] + epsilon)
equal <- nrow(plot_df) - above - below

cat(sprintf("%s 效应更强 (Tau_%s > Tau_%s): %d (%.1f%%)\n", 
            group1_name, group1_name, group2_name, above, above/nrow(plot_df)*100))
cat(sprintf("%s 效应更强 (Tau_%s > Tau_%s): %d (%.1f%%)\n", 
            group2_name, group2_name, group1_name, below, below/nrow(plot_df)*100))
cat(sprintf("近似相等 (|diff| <= 1e-6): %d (%.1f%%)\n", equal, equal/nrow(plot_df)*100))

# 配对 Wilcoxon 检验
wt <- wilcox.test(plot_df[[paste0("Tau_", group1_name)]], 
                  plot_df[[paste0("Tau_", group2_name)]], paired=TRUE)
cat(sprintf("配对 Wilcoxon 检验 P值: %g\n", wt$p.value))

# 四象限统计 (以0.5为界)
tau1_col <- paste0("Tau_", group1_name)
tau2_col <- paste0("Tau_", group2_name)
plot_df$Quadrant <- with(plot_df, 
  ifelse(get(tau1_col) >= 0.5 & get(tau2_col) < 0.5, "I",
  ifelse(get(tau1_col) >= 0.5 & get(tau2_col) >= 0.5, "II",
  ifelse(get(tau1_col) < 0.5 & get(tau2_col) < 0.5, "III", "IV")))
)
quad_counts <- table(plot_df$Quadrant)
cat("四象限统计 (阈值0.5):\n")
print(quad_counts)

# ------------------------- 6. 绘制散点图 (Hex bin + 颜色映射Delta) -------------------------
# 检查是否安装 hexbin
if (!requireNamespace("hexbin", quietly=TRUE)) {
  warning("包 'hexbin' 未安装，使用 geom_point 替代。安装: install.packages('hexbin')")
  geom_func <- geom_point(aes(color = DeltaTau), alpha=0.3, size=1)
} else {
  geom_func <- geom_hex(aes(fill = after_stat(density)), bins=80) + 
                scale_fill_viridis_c(option = "plasma", trans="log1p")
  # 但我们希望颜色映射Delta，但hexbin的fill是密度，不能直接映射Delta。我们可以叠加散点，但更好的方法是使用 geom_point 并设置颜色渐变。
  # 为了突出Delta，建议仍使用geom_point + color gradient，或者用geom_bin2d。这里选择更清晰的 geom_point 但透明度+颜色渐变。
}

# 改用 geom_point + 颜色渐变，同时支持大量点（透明度）
p <- ggplot(plot_df, aes_string(x = tau2_col, y = tau1_col, color = "DeltaTau")) +
  geom_point(alpha = 0.2, size = 0.8) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                        name = expression(Delta * Tau)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 0.8) +
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "grey50", alpha=0.5) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey50", alpha=0.5) +
  coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1)) +
  labs(
    title = paste0(group1_name, " vs ", group2_name, " Specificity (Tau)"),
    subtitle = paste0("Red dash: y = x  (above: ", group1_name, " stronger; below: ", group2_name, " stronger)\n",
                      "Quadrants (0.5 cutoff): I=", quad_counts["I"], ", II=", quad_counts["II"], 
                      ", III=", quad_counts["III"], ", IV=", quad_counts["IV"]),
    x = paste0(group2_name, " Tau"),
    y = paste0(group1_name, " Tau")
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.8),
    plot.subtitle = element_text(color = "grey30", size = 11)
  )

# 保存散点图
out_pdf <- paste0(opt$output, "_ScatterPlot.pdf")
ggsave(out_pdf, p, width = 8, height = 7)
cat("散点图已保存至:", out_pdf, "\n")

# ------------------------- 7. 绘制 DeltaTau 直方图 -------------------------
h <- ggplot(plot_df, aes(x = DeltaTau)) +
  geom_histogram(binwidth = 0.02, fill = "steelblue", color = "black", alpha=0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.8) +
  geom_vline(xintercept = median(plot_df$DeltaTau, na.rm=TRUE), 
             linetype = "dotted", color = "darkgreen", size = 0.8) +
  annotate("text", x = median(plot_df$DeltaTau, na.rm=TRUE), y = Inf, 
           label = paste("Median =", round(median(plot_df$DeltaTau, na.rm=TRUE), 3)), 
           hjust = -0.1, vjust = 1.5, color = "darkgreen", size = 4) +
  labs(
    title = paste0("Distribution of DeltaTau (", group1_name, " - ", group2_name, ")"),
    x = expression(Delta * Tau),
    y = "Number of transcripts"
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.8))

hist_pdf <- paste0(opt$output, "_DeltaTau_Histogram.pdf")
ggsave(hist_pdf, h, width = 7, height = 5)
cat("DeltaTau 直方图已保存至:", hist_pdf, "\n")

cat("\n===== 全部完成 =====\n")
