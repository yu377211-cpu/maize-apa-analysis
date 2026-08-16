#!/usr/bin/env Rscript

# ============================================================================
# 功能: DESeq2 差异表达分析，支持自定义设计公式（含协变量），raw count 输入。
# 用法:
#   Rscript DESeq2_covariate.R <count_matrix> <sample_info> <output_png> [design_formula]
# 参数:
#   count_matrix  : TSV，行为基因，列为样本，第一列为基因名（行名）
#   sample_info   : TSV，第一列为样本名（与 count_matrix 列名一致），其余列为样本属性
#                   文件名建议采用 "AvsB.condition" 格式，用于自动确定对比方向
#   output_png    : 热图输出路径（png）
#   design_formula: (可选) R 公式字符串，如 "~ genotype + tissue"
# 示例:
#   Rscript DESeq2_covariate.R counts.tsv GRootvsGShoot.condition heatmap.png "~ genotype + tissue"
# ============================================================================

library(data.table)
library(dplyr)
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(Cairo)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: ... <count_matrix> <sample_info> <output_png> [design_formula]")

# ---------- 1. 稳健读取表达矩阵 (raw count) ----------
cat("Reading count matrix...\n")
dt <- fread(args[1], header = TRUE, sep = "\t")
if (anyDuplicated(dt[[1]])) {
  warning("Gene names in first column are not unique. Adding suffix to duplicates.")
  dt[[1]] <- make.unique(dt[[1]])
}
count_matrix <- as.matrix(dt[, -1])
rownames(count_matrix) <- dt[[1]]
count_matrix <- round(count_matrix)
cat("Count matrix dimensions:", dim(count_matrix), "\n")

# ---------- 2. 读取样本信息 ----------
cat("Reading sample info...\n")
sample_info <- read.table(args[2], header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
cat("Sample info rows:", nrow(sample_info), "\n")

# ---------- 3. 取交集样本 ----------
common_samples <- intersect(colnames(count_matrix), rownames(sample_info))
if (length(common_samples) == 0) {
  stop("No common samples between count matrix columns and sample info rownames.")
}
cat("Number of common samples:", length(common_samples), "\n")
count_matrix <- count_matrix[, common_samples, drop = FALSE]
sample_info <- sample_info[common_samples, , drop = FALSE]
count_matrix <- count_matrix[, rownames(sample_info)]

# ---------- 4. 从样本信息文件名解析对比组 A vs B ----------
sample_file <- basename(args[2])
# 去除扩展名（支持 .condition, .txt, .tsv 等）
sample_name <- sub("\\.[^.]+$", "", sample_file)
A <- B <- NULL
if (grepl("vs", sample_name, ignore.case = TRUE)) {
  parts <- strsplit(sample_name, "vs", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    A <- trimws(parts[1])
    B <- trimws(parts[2])
    cat("Parsed comparison from filename:", A, "vs", B, "\n")
  } else {
    warning("Filename contains 'vs' but could not split into two parts. Falling back to auto-detection.")
  }
}

# ---------- 5. 确定分组列（group column） ----------
group_col <- NULL
if (!is.null(A) && !is.null(B)) {
  # 在 sample_info 中寻找同时包含 A 和 B 的列
  candidates <- c()
  for (col in colnames(sample_info)) {
    vals <- as.character(sample_info[[col]])
    if (A %in% vals && B %in% vals) {
      candidates <- c(candidates, col)
    }
  }
  if (length(candidates) == 1) {
    group_col <- candidates[1]
    cat("Detected group column (containing both", A, "and", B, "):", group_col, "\n")
  } else if (length(candidates) > 1) {
    warning("Multiple columns contain both A and B. Using first: ", candidates[1])
    group_col <- candidates[1]
  } else {
    warning("Could not find a column containing both", A, "and", B, ". Falling back to auto-detection.")
  }
}

# 若未通过文件名解析到 group_col，则使用原有的自动检测逻辑
if (is.null(group_col)) {
  col_names <- colnames(sample_info)
  group_col_candidates <- grep("tissue|group|condition", col_names, ignore.case = TRUE, value = TRUE)
  if (length(group_col_candidates) == 0) {
    group_col <- col_names[1]
    warning("No 'tissue/group/condition' column found, using first column (", group_col, ") as grouping.")
  } else {
    group_col <- group_col_candidates[1]
  }
  cat("Auto-detected group column (fallback):", group_col, "\n")
  # 此时 A/B 可能未定义，后续使用因子水平顺序
}

# ---------- 6. 设计公式 ----------
if (length(args) >= 4) {
  design_formula <- as.formula(args[4])
  cat("Using user-specified design:", deparse(design_formula), "\n")
  # 检查 group_col 是否在公式中
  if (!(group_col %in% all.vars(design_formula))) {
    stop("The detected group column '", group_col, "' is not present in the user-provided design formula. Please include it or adjust the formula.")
  }
} else {
  # 自动生成设计公式：协变量（除了 group_col 之外的列）+ group_col（放在最后）
  covariate_cols <- setdiff(colnames(sample_info), group_col)
  if (length(covariate_cols) > 0) {
    formula_str <- paste("~", paste(c(covariate_cols, group_col), collapse = " + "))
  } else {
    formula_str <- paste("~", group_col)
  }
  design_formula <- as.formula(formula_str)
  cat("Auto-generated design:", deparse(design_formula), "\n")
}

# ---------- 7. 构建 DESeq2 对象并运行 ----------
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_info,
                              design = design_formula)
dds <- DESeq(dds)

# ---------- 8. 提取差异结果 ----------
# 如果成功解析出 A 和 B，则直接使用它们构造 contrast
if (!is.null(A) && !is.null(B) && group_col %in% colnames(sample_info)) {
  # 确保 A 和 B 确实存在于该列中
  if (!(A %in% sample_info[[group_col]] && B %in% sample_info[[group_col]])) {
    stop("The parsed group values (", A, ", ", B, ") are not found in column '", group_col, "'.")
  }
  contrast <- c(group_col, B, A)
  cat("Using explicit contrast:", contrast[1], contrast[2], "vs", contrast[3], "\n")
} else {
  # 回退：使用因子水平第二个 vs 第一个
  group_levels <- as.character(unique(sample_info[[group_col]]))
  if (length(group_levels) < 2) stop(paste(group_col, "has fewer than 2 levels."))
  contrast <- c(group_col, group_levels[2], group_levels[1])
  cat("Using default contrast (levels order):", contrast[1], contrast[2], "vs", contrast[3], "\n")
}

res <- results(dds, contrast = contrast)
res_df <- as.data.frame(res, stringsAsFactors = FALSE)
res_df <- res_df[order(res_df$pvalue, -abs(res_df$log2FoldChange)), ]

# 筛选显著基因 (|log2FC| >= 1 & padj < 0.05；若 padj 缺失则用 pvalue)
pval_col <- if ("padj" %in% colnames(res_df)) "padj" else "pvalue"
res_sig <- res_df[which(abs(res_df$log2FoldChange) >= 1 & res_df[[pval_col]] < 0.05 & !is.na(res_df[[pval_col]])), ]

# 输出差异基因
outfile <- paste0(args[2], ".diff.txt")
write.table(res_sig, file = outfile, sep = "\t", col.names = NA, quote = FALSE)
cat("Differential results saved to", outfile, "\n")

# ---------- 9. 热图（使用 vst） ----------
if (nrow(res_sig) == 0) {
  cat("No significant genes. Heatmap skipped.\n")
  quit()
}
vsd <- vst(dds, blind = FALSE)
heatmap_data <- assay(vsd)[rownames(res_sig), , drop = FALSE]

show_rows <- if (nrow(heatmap_data) > 50) FALSE else TRUE
show_cols <- if (ncol(heatmap_data) > 50) FALSE else TRUE

annotation_col <- sample_info
CairoPNG(file = args[3], width = 1200, height = 800)
pheatmap(heatmap_data,
         annotation_col = annotation_col,
         show_rownames = show_rows,
         show_colnames = show_cols,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_method = "ward.D2",
         fontsize_row = 6,
         fontsize_col = 6,
         main = paste("DEPACs (|log2FC|>=1, adj.p<0.05)", contrast[2], "vs", contrast[3]))
dev.off()
cat("Heatmap saved to", args[3], "\n")
