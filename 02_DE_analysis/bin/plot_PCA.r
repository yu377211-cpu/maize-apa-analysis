#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(data.table)
})

option_list <- list(
  make_option(c("-e", "--expr"), type="character", help="PAC TPM matrix (row=PAC, col=sample)"),
  make_option(c("-m", "--meta"), type="character", help="Sample metadata"),
  make_option(c("-o", "--out"),  type="character", default="PCA")
)

opt <- parse_args(OptionParser(option_list=option_list))

# ===============================
# 1. Read data
# ===============================
expr <- fread(opt$expr, data.table = FALSE)
rownames(expr) <- expr[,1]
expr <- expr[,-1]

meta <- fread(opt$meta, data.table = FALSE)

# ===============================
# 2. Match samples
# ===============================
common.samples <- intersect(colnames(expr), meta$sample)
expr <- expr[, common.samples, drop=FALSE]
meta <- meta[match(common.samples, meta$sample), ]

cat("Samples retained:", ncol(expr), "\n")

# ===============================
# 3. Remove all-zero PACs
# ===============================
expr <- expr[rowSums(expr) > 0, ]
cat("PACs after removing all-zero:", nrow(expr), "\n")

# ===============================
# 4. log2 transform
# ===============================
expr.log <- log2(expr + 1)

# ===============================
# 5. Variance distribution (QC)
# ===============================
pac.var.raw <- apply(expr.log, 1, var)

var.df <- data.frame(variance = pac.var.raw)

p.var <- ggplot(var.df, aes(x = variance)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "black") +
  scale_x_log10() +
  theme_bw() +
  labs(
    title = "PAC variance distribution (log10)",
    x = "Variance of log2(TPM + 1)",
    y = "PAC count"
  )

ggsave(paste0(opt$out, "_variance_dist.png"), p.var, width = 6, height = 5)
ggsave(paste0(opt$out, "_variance_dist.pdf"), p.var, width = 6, height = 5)

# ===============================
# 6. Variance filtering
# ===============================
expr.log <- expr.log[pac.var.raw > 0, ]
cat("PACs after variance filtering:", nrow(expr.log), "\n")

# ===============================
# 7. PCA
# ===============================
pca <- prcomp(t(expr.log), center = TRUE, scale. = TRUE)

pca.df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  tissue = meta$tissue,
  subpop = meta$subpop
)

pc1.var <- round(summary(pca)$importance[2,1] * 100, 2)
pc2.var <- round(summary(pca)$importance[2,2] * 100, 2)

# ===============================
# 8. Plot PCA
# ===============================
p <- ggplot(pca.df, aes(PC1, PC2, color = tissue, shape = subpop)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_shape_manual(values = c(16, 17, 15, 3, 4, 7, 8)) +
  theme_bw() +
  labs(
    x = paste0("PC1 (", pc1.var, "%)"),
    y = paste0("PC2 (", pc2.var, "%)"),
    title = "Expression-based PCA (3′ seq)"
  )

ggsave(paste0(opt$out, ".png"), p, width = 6, height = 5)
ggsave(paste0(opt$out, ".pdf"), p, width = 6, height = 5)

