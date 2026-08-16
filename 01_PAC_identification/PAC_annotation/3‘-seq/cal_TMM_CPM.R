#!/usr/bin/env Rscript

# ============================================================
# 脚本功能：从原始 PAC count 矩阵计算 TMM 标准化的 CPM
# 输入格式：制表符分隔，第一列为 PACid，后续列为各样本的 raw counts
# 输出：两个文件
#   1) TMM_CPM_matrix.txt  : TMM 标准化的 CPM 值（未 log 变换）
#   2) TMM_log2CPM_matrix.txt : log2(CPM + 1) 变换后的值（用于 PCA/热图）
# 使用方法：Rscript this_script.R input_matrix.txt
# ============================================================

# 加载必要的包
library(edgeR)

# 从命令行获取输入文件名
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript cal_TMM_CPM.R <input_count_matrix.txt>")
}
input_file <- args[1]

# 读取数据
cat("Reading count matrix from", input_file, "...\n")
count_data <- read.table(input_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# 检查数据是否为整数矩阵
if (!all(count_data == round(count_data))) {
  warning("Count matrix contains non-integer values. edgeR requires integer counts.")
}

# 将数据框转换为矩阵（edgeR 要求）
count_matrix <- as.matrix(count_data)

# 创建 DGEList 对象
dge <- DGEList(counts = count_matrix)

# 计算 TMM 标准化因子
cat("Calculating TMM normalization factors...\n")
dge <- calcNormFactors(dge, method = "TMM")

# 获取 TMM 标准化的 CPM（未 log 变换）
cat("Computing TMM-normalized CPM...\n")
tmm_cpm <- cpm(dge, log = FALSE)

# 获取 log2(CPM + 1) 变换后的值（用于可视化）
# prior.count = 1 避免 log2(0) 的问题
log2_tmm_cpm <- cpm(dge, log = TRUE, prior.count = 1)

# 输出结果
output_base <- gsub("\\.txt$|\\.tsv$", "", input_file)  # 去掉扩展名

cpm_file <- paste0(output_base, "_TMM_CPM.txt")
log2_file <- paste0(output_base, "_TMM_log2CPM.txt")

# 写入制表符分隔的文件，保留行名
cat("Writing TMM CPM to", cpm_file, "...\n")
write.table(tmm_cpm, file = cpm_file, sep = "\t", quote = FALSE, col.names = NA)

cat("Writing log2(TMM CPM + 1) to", log2_file, "...\n")
write.table(log2_tmm_cpm, file = log2_file, sep = "\t", quote = FALSE, col.names = NA)

cat("Done!\n")
