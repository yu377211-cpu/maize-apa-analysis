#!/usr/bin/env Rscript
# ================================================================
# 组织特异性 APA 分析（PDUI 矩阵）
# 用法：Rscript specific_APA.R -i PDUI_matrix.txt -g sample_groups.txt -o output_prefix
# ================================================================

# 加载必要包（请提前安装）
suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(pheatmap)
  library(BiocParallel)
  library(ggplot2)  # 用于调色
})

# ------------------------- 解析命令行参数 -------------------------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "PDUI 矩阵文件（tab分隔，第一列为转录本ID，其余列为样本）", metavar = "file"),
  make_option(c("-g", "--group"), type = "character", default = NULL,
              help = "分组文件（tab分隔，第一列样本ID，第二列组织名，有标题行）", metavar = "file"),
  make_option(c("-o", "--output"), type = "character", default = "APA_result",
              help = "输出文件前缀 [默认 %default]", metavar = "prefix"),
  make_option(c("-t", "--threshold"), type = "double", default = 0.05,
              help = "FDR 显著性阈值 [默认 %default]", metavar = "num"),
  make_option(c("-c", "--cores"), type = "integer", default = 4,
              help = "并行计算核心数 [默认 %default]", metavar = "num"),
  make_option(c("-m", "--method"), type = "character", default = "kruskal",
              help = "检验方法: kruskal 或 anova [默认 %default]", metavar = "method"),
  make_option(c("--agg"), action = "store_true", default = FALSE,
              help = "是否使用组织聚合矩阵进行聚类（推荐） [默认启用，若加此参数则使用原始样本矩阵]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$group)) {
  print_help(opt_parser)
  stop("必须提供 -i 和 -g 参数", call. = FALSE)
}

# ------------------------- 1. 读入数据 -------------------------
cat("正在读取 PDUI 矩阵...\n")
pdui <- fread(opt$input, header = TRUE, sep = "\t", data.table = FALSE)
rownames(pdui) <- pdui[, 1]
pdui <- as.matrix(pdui[, -1])

cat("正在读取分组文件...\n")
groups <- fread(opt$group, header = TRUE, sep = "\t", data.table = FALSE)
# 取前两列，第一列样本ID，第二列组织
sample_id <- groups[, 1]
tissue <- groups[, 2]
names(tissue) <- sample_id

# 确保矩阵列名与分组一致
common <- intersect(colnames(pdui), names(tissue))
if (length(common) == 0) stop("样本ID完全不匹配！")
pdui <- pdui[, common, drop = FALSE]
tissue <- tissue[common]
tissue_factor <- factor(tissue)

n_samples <- ncol(pdui)
n_tissues <- nlevels(tissue_factor)
cat(sprintf("共有 %d 个样本，%d 个组织类型\n", n_samples, n_tissues))

# ------------------------- 2. 差异检验（并行） -------------------------
# 自定义检验函数
if (opt$method == "kruskal") {
  test_func <- function(row_vals) {
    valid <- !is.na(row_vals)
    if (sum(valid) < 3) return(NA)  # 至少3个有效值
    x <- row_vals[valid]
    g <- tissue_factor[valid]
    if (any(table(g) < 2)) return(NA)
    tryCatch(kruskal.test(x, g)$p.value, error = function(e) NA)
  }
} else if (opt$method == "anova") {
  test_func <- function(row_vals) {
    valid <- !is.na(row_vals)
    if (sum(valid) < 3) return(NA)
    x <- row_vals[valid]
    g <- tissue_factor[valid]
    if (any(table(g) < 2)) return(NA)
    tryCatch(summary(aov(x ~ g))[[1]]$`Pr(>F)`[1], error = function(e) NA)
  }
} else {
  stop("method 必须是 'kruskal' 或 'anova'")
}

cat("正在进行", opt$method, "检验（并行", opt$cores, "核）...\n")
BPPARAM <- MulticoreParam(workers = opt$cores, progressbar = TRUE)
p_values <- bplapply(seq_len(nrow(pdui)), function(i) {
  test_func(pdui[i, ])
}, BPPARAM = BPPARAM)
p_values <- unlist(p_values)

cat("有效 p 值数量:", sum(!is.na(p_values)), "\n")

# 校正
padj <- rep(NA, length(p_values))
p_adj <- p.adjust(p_values[!is.na(p_values)], method = "fdr")
padj[!is.na(p_values)] <- p_adj

sig_idx <- which(padj < opt$threshold & !is.na(padj))
cat(sprintf("显著转录本数量 (FDR < %g): %d\n", opt$threshold, length(sig_idx)))

if (length(sig_idx) == 0) {
  cat("没有显著转录本，程序退出。\n")
  quit(save = "no")
}

# 保存显著结果
sig_data <- pdui[sig_idx, , drop = FALSE]
sig_results <- data.frame(
  transcript = rownames(sig_data),
  p_value = p_values[sig_idx],
  padj = padj[sig_idx]
)
out_file <- paste0(opt$output, "_sig_APA.txt")
write.table(sig_results, file = out_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("显著转录本列表已保存至:", out_file, "\n")

# ------------------------- 3. 聚类热图（组织聚合模式，推荐） -------------------------
# 若用户加了 --agg 参数，则使用原始样本矩阵（列数多，可能需要处理NA）
if (!opt$agg) {
  cat("使用组织聚合矩阵进行聚类（按组织均值，忽略NA）...\n")
  # 按组织计算均值（宽格式转长格式再聚合）
  # 方法：将矩阵转为数据框，split by tissue，计算行均值
  # 高效做法：使用 rowsum + 分组
  # 先构建一个组织向量与列对应
  tissue_vec <- tissue_factor
  # 获取每个组织的列索引
  tissue_levels <- levels(tissue_vec)
  agg_list <- lapply(tissue_levels, function(lev) {
    cols <- which(tissue_vec == lev)
    if (length(cols) == 1) {
      return(sig_data[, cols, drop = FALSE])
    } else {
      # 按行计算均值，na.rm = TRUE
      rowMeans(sig_data[, cols, drop = FALSE], na.rm = TRUE)
    }
  })
  # 合并成矩阵（列名=组织）
  agg_mat <- do.call(cbind, agg_list)
  colnames(agg_mat) <- tissue_levels
  # 检查是否有整行全为NA（理论上显著转录本至少有些值）
  na_rows <- apply(agg_mat, 1, function(x) all(is.na(x)))
  if (any(na_rows)) {
    cat("剔除", sum(na_rows), "个在聚合后全为NA的行\n")
    agg_mat <- agg_mat[!na_rows, , drop = FALSE]
  }
  # 标准化（按行）
  scaled <- t(scale(t(agg_mat)))
  # 处理极端值
  scaled[scaled > 3] <- 3
  scaled[scaled < -3] <- -3
  # 如果显著转录本太多，按方差筛选
  if (nrow(scaled) > 2000) {
    cat("显著转录本过多 (", nrow(scaled), ")，将按行方差选取前 2000 个显示\n")
    row_var <- apply(scaled, 1, var, na.rm = TRUE)
    keep <- order(row_var, decreasing = TRUE)[1:2000]
    scaled <- scaled[keep, , drop = FALSE]
  }
  # 列注释（组织）
  annotation_col <- data.frame(Tissue = colnames(scaled))
  rownames(annotation_col) <- colnames(scaled)
  # 颜色
  tissue_cols <- rainbow(n_tissues)
  names(tissue_cols) <- tissue_levels
  ann_colors <- list(Tissue = tissue_cols)
  
  # 绘制热图（列数为组织数，简洁）
  pdf_file <- paste0(opt$output, "_heatmap_agg.pdf")
  pdf(pdf_file, width = 8 + 0.3 * n_tissues, height = 10)
  pheatmap(scaled,
           main = "组织特异性 APA (PDUI) - 组织均值",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           show_rownames = FALSE,
           show_colnames = TRUE,
           fontsize_col = 10,
           color = colorRampPalette(c("navy", "white", "red"))(100),
           scale = "none")
  dev.off()
  cat("热图已保存至:", pdf_file, "\n")
} else {
  # 原始样本级热图（不推荐，列数过多，且需处理NA）
  cat("使用原始样本矩阵进行聚类（列数多，注意内存）...\n")
  # 需要对每行进行标准化，但保留NA，然后聚类时使用可处理NA的距离
  # 此处简化：用0填充NA（仅用于聚类），但会改变生物学解释，不建议。
  # 更稳健：使用相关性距离（use="pairwise.complete.obs"），但pheatmap不直接支持。
  # 我们可以自定义距离函数，但比较复杂。这里给出替代方案：
  # 对于样本级聚类，建议先对每个转录本进行组织内均值填充（或中位数），再标准化。
  # 为了快速解决，我们直接使用聚合方式（忽略 --agg 参数）。
  warning("原始样本级聚类可能因NA值失败，建议使用默认的组织聚合方式。")
  # 简单处理：用0填充NA（只用于热图显示）
  mat_fill <- sig_data
  mat_fill[is.na(mat_fill)] <- 0
  scaled <- t(scale(t(mat_fill)))
  scaled[scaled > 3] <- 3
  scaled[scaled < -3] <- -3
  if (nrow(scaled) > 2000) {
    row_var <- apply(scaled, 1, var, na.rm = TRUE)
    keep <- order(row_var, decreasing = TRUE)[1:2000]
    scaled <- scaled[keep, , drop = FALSE]
  }
  annotation_col <- data.frame(Tissue = tissue_factor)
  rownames(annotation_col) <- colnames(scaled)
  tissue_cols <- rainbow(n_tissues)
  names(tissue_cols) <- levels(tissue_factor)
  ann_colors <- list(Tissue = tissue_cols)
  
  pdf_file <- paste0(opt$output, "_heatmap_raw.pdf")
  pdf(pdf_file, width = 25, height = 12)
  pheatmap(scaled,
           main = "组织特异性 APA (PDUI) - 原始样本（NA填0）",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           show_rownames = FALSE,
           show_colnames = FALSE,
           fontsize_row = 4,
           color = colorRampPalette(c("navy", "white", "red"))(100),
           scale = "none")
  dev.off()
  cat("热图已保存至:", pdf_file, "\n")
}

cat("分析完成！\n")
