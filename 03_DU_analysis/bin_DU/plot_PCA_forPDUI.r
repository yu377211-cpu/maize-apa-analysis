#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

# --------------------------
# 1. 参数
# --------------------------
option_list <- list(
  make_option(c("-e","--expr"), type="character", help="PDUI matrix (rows=transcripts, cols=samples)"),
  make_option(c("-m","--meta"), type="character", help="Sample metadata (sample, tissue, subpop)"),
  make_option(c("-o","--out"), type="character", default="PDUI_PCA", help="Output prefix"),
  make_option(c("--min_prop"), type="numeric", default=0.2, help="Minimum proportion of non-NA per transcript"),
  make_option(c("--width"), type="numeric", default=6, help="Figure width"),
  make_option(c("--height"), type="numeric", default=5, help="Figure height")
)
opt <- parse_args(OptionParser(option_list=option_list))

# --------------------------
# 2. 读入数据
# --------------------------
pdui <- fread(opt$expr, data.table=FALSE)
rownames(pdui) <- pdui[,1]
pdui <- pdui[,-1, drop=FALSE]

meta <- fread(opt$meta, data.table=FALSE)

# --------------------------
# 3. 样本匹配
# --------------------------
common.samples <- intersect(colnames(pdui), meta$sample)
pdui <- pdui[, common.samples, drop=FALSE]
meta <- meta[match(common.samples, meta$sample), ]
cat("Samples retained:", ncol(pdui), "\n")

# --------------------------
# 4. 删除全 NA 转录本
# --------------------------
keep.nonNA <- rowSums(!is.na(pdui)) > 0
pdui <- pdui[keep.nonNA, ]
cat("Transcripts after removing all-NA:", nrow(pdui), "\n")

# --------------------------
# 5. coverage filtering
# --------------------------
keep.coverage <- rowMeans(!is.na(pdui)) >= opt$min_prop
pdui <- pdui[keep.coverage, ]
cat("Transcripts after coverage filtering:", nrow(pdui), "\n")

# --------------------------
# 6. QC histogram: PDUI distribution
# --------------------------
pdui.values <- numeric(0)
for (i in 1:ncol(pdui)) {
  col_vals <- as.numeric(pdui[[i]])
  col_vals <- col_vals[!is.na(col_vals)]
  pdui.values <- c(pdui.values, col_vals)
}

if (length(pdui.values) > 0) {
  pdui.long <- data.frame(value = pdui.values)
  
  p1 <- ggplot(pdui.long, aes(x=value)) +
    geom_histogram(bins=50, fill="darkorange", color="black") +
    theme_bw() +
    labs(title="Global PDUI distribution", x="PDUI", y="Count")
  
  ggsave(paste0(opt$out, "_PDUI_dist.png"), p1, width=opt$width, height=opt$height)
  ggsave(paste0(opt$out, "_PDUI_dist.pdf"), p1, width=opt$width, height=opt$height)
  cat("PDUI distribution plot saved\n")
}

# --------------------------
# 7. 数据预处理和NA填充
# --------------------------
pdui.mat <- as.matrix(pdui)
mode(pdui.mat) <- "numeric"

# NA填充（每行 median）
pdui.impute <- matrix(NA, nrow=nrow(pdui.mat), ncol=ncol(pdui.mat))
rownames(pdui.impute) <- rownames(pdui.mat)
colnames(pdui.impute) <- colnames(pdui.mat)

for (i in 1:nrow(pdui.mat)) {
  row_vals <- pdui.mat[i, ]
  med_val <- median(row_vals, na.rm=TRUE)
  
  if (is.na(med_val)) {
    pdui.impute[i, ] <- 0.5
  } else {
    row_vals[is.na(row_vals)] <- med_val
    pdui.impute[i, ] <- row_vals
  }
}

# --------------------------
# 8. QC histogram: variance
# --------------------------
var.values <- apply(pdui.impute, 1, var)
valid.var <- var.values[!is.na(var.values) & var.values > 0]

if (length(valid.var) > 0) {
  var.df <- data.frame(variance = valid.var)
  
  p2 <- ggplot(var.df, aes(x=variance)) +
    geom_histogram(bins=100, fill="steelblue", color="black") +
    theme_bw() +
    labs(title="PDUI variance distribution", x="Variance of PDUI", y="Transcript count")
  
  if (max(valid.var) / min(valid.var) > 100) {
    p2 <- p2 + scale_x_log10()
  }
  
  ggsave(paste0(opt$out, "_PDUI_variance_dist.png"), p2, width=opt$width, height=opt$height)
  ggsave(paste0(opt$out, "_PDUI_variance_dist.pdf"), p2, width=opt$width, height=opt$height)
  cat("PDUI variance distribution plot saved\n")
}

# --------------------------
# 9. 去掉零方差
# --------------------------
pdui.var <- apply(pdui.impute, 1, var)
keep.var <- !is.na(pdui.var) & pdui.var > 0
pdui.filtered <- pdui.impute[keep.var, , drop=FALSE]
cat("Transcripts after variance filtering:", nrow(pdui.filtered), "\n")

if (nrow(pdui.filtered) < 2 || ncol(pdui.filtered) < 2) {
  stop("Not enough data after filtering.")
}

# --------------------------
# 10. PCA
# --------------------------
cat("Final data dimensions:", dim(pdui.filtered), "\n")

tryCatch({
  pca <- prcomp(t(pdui.filtered), center=TRUE, scale.=FALSE)
  pc1.var <- round(summary(pca)$importance[2,1]*100, 2)
  pc2.var <- round(summary(pca)$importance[2,2]*100, 2)
  
  pca.df <- data.frame(
    sample = rownames(pca$x),
    PC1 = pca$x[,1],
    PC2 = pca$x[,2],
    tissue = meta$tissue,
    subpop = meta$subpop
  )
  
  # --------------------------
  # 11. PCA plot - 与参考格式保持一致
  # --------------------------
  p <- ggplot(pca.df, aes(PC1, PC2, color = tissue, shape = subpop)) +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_shape_manual(values = c(16, 17, 15, 3, 4, 7, 8)) +
    theme_bw() +
    labs(
      x = paste0("PC1 (", pc1.var, "%)"),
      y = paste0("PC2 (", pc2.var, "%)"),
      title = "Usage based PCA(3'seq)"  # 按照您的要求修改标题
    )
  
  ggsave(paste0(opt$out, ".png"), p, width = 6, height = 5)
  ggsave(paste0(opt$out, ".pdf"), p, width = 6, height = 5)
  
  # 保存PCA坐标
  write.table(pca.df, paste0(opt$out, "_coordinates.txt"), 
              sep="\t", quote=FALSE, row.names=FALSE)
  
  cat("\n✅ Analysis completed!\n")
  cat("PCA plot saved:", paste0(opt$out, ".png"), "\n")
  
}, error = function(e) {
  cat("\n❌ PCA failed:", e$message, "\n")
})
