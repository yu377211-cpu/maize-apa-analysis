# 读取协变量文件（第一列为协变量名，第一行为样本名）
cov <- read.table("./Matrix_eQTL/Covariate_matrix.txt", header=TRUE, check.names=FALSE)
rownames(cov) <- cov[,1]
cov <- cov[,-1]  # 此时行为协变量，列为样本

# 转换为数值矩阵
mat <- as.matrix(cov)
mode(mat) <- "numeric"

# 1. 计算方差，保留大于阈值（如 1e-8）的行
row_var <- apply(mat, 1, var, na.rm = TRUE)
keep <- row_var > 1e-8
cat("保留前：", nrow(mat), "行，保留后：", sum(keep), "行\n")

mat_filtered <- mat[keep, ]

# 2. 剔除高度共线的行（相关系数 > 0.9999）
if (nrow(mat_filtered) > 1) {
  cor_mat <- cor(t(mat_filtered), use = "pairwise.complete.obs")
  high_cor <- which(abs(cor_mat) > 0.9999 & upper.tri(cor_mat), arr.ind = TRUE)
  if (nrow(high_cor) > 0) {
    remove_idx <- unique(high_cor[, 2])
    mat_filtered <- mat_filtered[-remove_idx, ]
    cat("剔除共线行：", length(remove_idx), "个\n")
  }
}

# 3. 写回新文件（行名作为第一列）
write.table(mat_filtered, file = "./Matrix_eQTL/Covariate_matrix_filtered.txt", 
            quote = FALSE, sep = "\t", col.names = NA)
