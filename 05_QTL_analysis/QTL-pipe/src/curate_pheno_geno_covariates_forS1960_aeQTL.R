#!/opt/app/languages/R-3.6.3/bin/Rscript
# curate_pheno_geno_covariates_forS1960_aeQTL.R
# 适用于 3'aeQTL 分析，输入为 TMM_CPM 矩阵（含0值）
# 过滤在原始 CPM 层面进行，过滤后进行 log2(CPM+1) 变换
# 修正了 PLINK .eigenvec 文件解析（FID + IID）

library(optparse)
library(dplyr)
library(peer)

option_list <- list(
  make_option(c("-p", "--pheno_data"), type="character", 
              default="TMM_CPM_matrix.txt", 
              action="store", help="TMM_CPM matrix (raw CPM, with zeros)"),
  make_option(c("-g","--geno_pca"), type="character", 
              default="./Matrix_eQTL/genotype_pca.eigenvec", 
              action="store", help="PLINK format .eigenvec (FID + IID + PCs)"),
  make_option(c("-c","--known_covs"), type="character", default="NA", 
              action="store", help="known covariates file (first column: sample ID)"),
  make_option(c("-n","--top_N_pca"), type="integer", default="5", 
              action="store", help="top N PCA components to use"),
  make_option(c("--cpm_threshold"), type="double", default=5, 
              action="store", help="CPM threshold for low-expression filter"),
  make_option(c("--sample_pct"), type="double", default=0.3, 
              action="store", help="minimum fraction of samples with CPM >= threshold")
)
opt <- parse_args(OptionParser(option_list=option_list, usage="usage: %prog [options]"))

cat('Arguments:\n',
    '--pheno_data', opt$pheno_data, '\n',
    '--geno_pca', opt$geno_pca, '\n',
    '--known_covs', opt$known_covs, '\n',
    '--top_N_pca', opt$top_N_pca, '\n',
    '--cpm_threshold', opt$cpm_threshold, '\n',
    '--sample_pct', opt$sample_pct, '\n')

# ============================================================
# 1. 构建协变量（PC + known covariates）
# ============================================================
# ---- 读取 PCA（PLINK 格式：FID + IID + PCs） ----
gt_pca_raw <- read.table(opt$geno_pca, header=F, sep=" ", stringsAsFactors=F)

N <- as.integer(opt$top_N_pca)

# 检查列数是否足够
if(ncol(gt_pca_raw) < (N + 2)) {
  stop("ERROR: genotype_pca.eigenvec has fewer than ", N+2, " columns. Check -n parameter.")
}

gt_pca <- data.frame(
  subject_id = gt_pca_raw[, 2],  # IID（第2列）
  gt_pca_raw[, 3:(N+2)]          # PC1 ~ PCn
)
colnames(gt_pca) <- c("subject_id", paste0("PC_", 1:N))

cat("Genotype PCA samples:", nrow(gt_pca), "\n")
cat("Using top", N, "PCs\n")

# ---- 添加已知协变量 ----
if(opt$known_covs != "NA"){
  known_cov <- read.table(opt$known_covs, header=T, sep="\t", stringsAsFactors=F)
  
  # 确保第一列名为 subject_id
  colnames(known_cov)[1] <- "subject_id"
  
  # 将非数值列转为数值因子
  for(i in 2:ncol(known_cov)){
    if(class(known_cov[,i]) == "character"){
      known_cov[,i] <- as.numeric(as.factor(known_cov[,i]))
    }
  }
  
  cat("Known covariates samples:", nrow(known_cov), "\n")
  
  # 合并
  gt_pca <- merge(gt_pca, known_cov, by="subject_id")
  
  cat("Matched samples after merging:", nrow(gt_pca), "\n")
}

if(nrow(gt_pca) == 0){
  stop("ERROR: Covariate matrix is empty. Check sample IDs in -g and -c files.")
}

# 转换为矩阵
rownames(gt_pca) <- gt_pca$subject_id
gt_pca <- as.matrix(gt_pca[, -1])
cat('Final covariate dimensions:', dim(gt_pca), '\n')

# ============================================================
# 2. 读取表型数据（TMM_CPM，含0值）
# ============================================================
cat("Reading phenotype matrix (TMM_CPM):", opt$pheno_data, "\n")
pheno_mat <- read.table(opt$pheno_data, header=T, sep="\t", check.names=F, stringsAsFactors=F)

# 提取与协变量样本匹配的列
sample_ids <- rownames(gt_pca)
pheno_sel <- pheno_mat %>% dplyr::select(all_of(sample_ids))
pheno_sel <- as.matrix(pheno_sel)
rownames(pheno_sel) <- pheno_mat[, 1]

cat("Original phenotype dimensions:", dim(pheno_sel), "\n")

# ---- 检查输入是否为 log 值（陷阱检测） ----
if(any(pheno_sel > 20, na.rm=TRUE)) {
  cat("WARNING: Detected values > 20 in input matrix.\n")
  cat("         Are you sure this is raw TMM_CPM (not log-transformed)?\n")
  cat("         If this is log2(CPM+1), please use raw CPM as input.\n")
}

# ---- 检查 NA ----
if(any(is.na(pheno_sel))){
  na_genes <- rownames(pheno_sel)[rowSums(is.na(pheno_sel)) > 0]
  stop("ERROR: Found NA values in phenotype matrix.\n",
       "  Affected genes (first 10): ", paste(head(na_genes, 10), collapse=", "), "\n",
       "  Please check your input data.")
}

# ============================================================
# 3. 低表达过滤（直接在原始 CPM 上进行）
# ============================================================
cpm_threshold <- opt$cpm_threshold
sample_pct <- opt$sample_pct
n_samples <- ncol(pheno_sel)
min_samples <- max(3, ceiling(n_samples * sample_pct))

cat("Filter: CPM >= ", cpm_threshold, " in >= ", min_samples, " samples (", 
    sample_pct*100, "% of ", n_samples, " samples)\n")

keep_genes <- rowSums(pheno_sel >= cpm_threshold) >= min_samples
pheno_sel_raw <- pheno_sel[keep_genes, ]

cat("Genes retained (CPM >= ", cpm_threshold, "): ", nrow(pheno_sel_raw), 
    " out of ", length(keep_genes), "\n")

if(nrow(pheno_sel_raw) == 0){
  stop("ERROR: No genes passed the expression filter. Please relax threshold.")
}

# ============================================================
# 4. 转换为 log2(CPM+1)（用于 PEER）
# ============================================================
pheno_log <- log2(pheno_sel_raw + 1)
cat("Transformed to log2(CPM+1) for PEER. Shape:", dim(pheno_log), "\n")

# ---- 输出表型矩阵（log2(CPM+1)） ----
id_order <- colnames(pheno_log)
pheno_out <- cbind(rownames(pheno_log), pheno_log)
colnames(pheno_out) <- c("Gene", id_order)
write.table(pheno_out, file="./Matrix_eQTL/Phenotype_matrix.txt", 
            row.names=F, col.names=T, quote=F, sep="\t")
cat("Output phenotype matrix: ./Matrix_eQTL/Phenotype_matrix.txt\n")

# ============================================================
# 5. PEER 因子估计（使用 log2(CPM+1)）
# ============================================================
cat("Running PEER on log2(CPM+1) matrix...\n")
model <- PEER()

# 设置协变量
PEER_setCovariates(model, gt_pca)

# 设置表型（转置：样本为列，基因为行）
PEER_setPhenoMean(model, t(pheno_log))

# 确定 PEER 因子数量
if (n_samples < 150) {
  numcov <- 15
} else if (n_samples < 250) {
  numcov <- 30
} else {
  numcov <- 35
}
cat("Number of PEER factors:", numcov, "\n")

PEER_setNk(model, numcov)
PEER_update(model)

# ---- 诊断图 ----
pdf('peer.diag.pdf', width=6, height=8)
PEER_plotModel(model)
dev.off()
cat("PEER diagnostic plot: peer.diag.pdf\n")

# ---- 提取因子 ----
factors <- t(PEER_getX(model))
cat(class(factors), "\n")
cat(dim(factors), "\n")
rownames(factors) <- c(colnames(gt_pca), paste0("PEER_", 1:numcov))
colnames(factors) <- colnames(pheno_log)

# ---- 输出残差（可选） ----
residuals <- t(PEER_getResiduals(model))
rownames(residuals) <- rownames(pheno_log)
colnames(residuals) <- colnames(pheno_log)
residual.df <- data.frame(id = rownames(residuals), residuals, check.names = FALSE)
write.table(residual.df, 
            file="pdui.peer.residuals.txt", 
            row.names=F, col.names=c("id", colnames(residuals)), 
            quote=F, sep='\t')
cat("PEER residuals: pdui.peer.residuals.txt\n")

# ============================================================
# 6. 基因型和协变量矩阵输出
# ============================================================
cat("Loading genotype_matrix.bed...\n")
gt_mat <- read.table("./Matrix_eQTL/genotype_matrix.bed", header=T, sep="\t", check.names=F)
gt_mat.reorder <- gt_mat %>% dplyr::select("id", all_of(id_order))
write.table(gt_mat.reorder, file="./Matrix_eQTL/Genotype_matrix.txt", 
            quote=F, sep="\t", row.names=F, col.names=T)
cat("Output genotype matrix: ./Matrix_eQTL/Genotype_matrix.txt\n")

# ---- 协变量矩阵 ----
factors.df <- data.frame(id = rownames(factors), factors, check.names = FALSE)
colnames(factors.df) <- c("id", colnames(factors))
factors.reorder <- factors.df %>% dplyr::select("id", all_of(id_order))
write.table(factors.reorder, file="./Matrix_eQTL/Covariate_matrix.txt", 
            row.names=F, quote=F, sep='\t', col.names=T)
cat("Output covariates matrix: ./Matrix_eQTL/Covariate_matrix.txt\n")

# ============================================================
# 7. 总结信息
# ============================================================
cat("\n=== Done ===\n")
cat("Input PACs:", nrow(pheno_sel), "\n")
cat("PACs retained after filter (CPM>=", cpm_threshold, " in >=", min_samples, " samples):", 
    nrow(pheno_log), "\n")
cat("Samples:", n_samples, "\n")
cat("PEER factors:", numcov, "\n")
cat("\nOutput files:\n")
cat("  ./Matrix_eQTL/Phenotype_matrix.txt (log2(CPM+1), n=", nrow(pheno_log), ")\n")
cat("  ./Matrix_eQTL/Genotype_matrix.txt\n")
cat("  ./Matrix_eQTL/Covariate_matrix.txt\n")
cat("  pdui.peer.residuals.txt\n")
cat("  peer.diag.pdf\n")
