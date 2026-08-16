#!/usr/bin/env Rscript

# ==========================================================
# 统计不同过滤条件下剩余 PAC 数量
# 用法: Rscript count_PAC_by_conditions.R -i <PDUI_matrix>
# ==========================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "PDUI matrix file (rows=PACs, columns=samples, first column is PAC ID)")
)

opt <- parse_args(OptionParser(option_list = option_list))

# 读取数据
cat("Reading:", opt$input, "\n")
mat <- fread(opt$input, data.table = FALSE)
rownames(mat) <- mat[, 1]
mat <- mat[, -1]
mat <- as.matrix(mat)
mode(mat) <- "numeric"

n_total <- nrow(mat)
n_samples <- ncol(mat)
cat("Total PACs:", n_total, "\n")
cat("Total samples:", n_samples, "\n")

# ==========================================================
# 逐行计算 4 个核心指标
# ==========================================================
na_ratio <- rowMeans(is.na(mat))  # 缺失率

# 方差（处理全缺失行：返回 NA，后续过滤时排除）
var_val <- apply(mat, 1, var, na.rm = TRUE)

# IQR（处理全缺失行：返回 NA）
iqr_val <- apply(mat, 1, IQR, na.rm = TRUE)

# 非极端值比例: (0 < PDUI < 1) 的样本数 / 非缺失样本数
non_extreme_count <- rowSums(mat > 1e-9 & mat < 1 - 1e-9, na.rm = TRUE)
non_na_count <- rowSums(!is.na(mat))
non_extreme_ratio <- non_extreme_count / non_na_count
# 如果全部缺失，比例设为 0（而不是 NaN）
non_extreme_ratio[is.na(non_extreme_ratio)] <- 0

# ==========================================================
# 4 个递进条件的布尔向量
# ==========================================================
cond1 <- na_ratio < 0.5
cond2 <- cond1 & var_val > 0
cond3 <- cond2 & iqr_val > 0
cond4 <- cond3 & non_extreme_ratio > 0.05

# ==========================================================
# 结果汇总
# ==========================================================
results <- data.frame(
  Condition = c(
    "NA < 0.5",
    "NA < 0.5 + Var > 0",
    "NA < 0.5 + Var > 0 + IQR > 0",
    "NA < 0.5 + Var > 0 + IQR > 0 + NonExtreme > 0.05"
  ),
  Remaining_PAC = c(
    sum(cond1, na.rm = TRUE),
    sum(cond2, na.rm = TRUE),
    sum(cond3, na.rm = TRUE),
    sum(cond4, na.rm = TRUE)
  ),
  Percentage = c(
    round(sum(cond1, na.rm = TRUE) / n_total * 100, 2),
    round(sum(cond2, na.rm = TRUE) / n_total * 100, 2),
    round(sum(cond3, na.rm = TRUE) / n_total * 100, 2),
    round(sum(cond4, na.rm = TRUE) / n_total * 100, 2)
  )
)

cat("\n========================================\n")
print(results)
cat("========================================\n")

# 可选：将结果保存到文件
write.table(results, file = "PAC_filter_counts.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
cat("\nResults saved to: PAC_filter_counts.txt\n")
