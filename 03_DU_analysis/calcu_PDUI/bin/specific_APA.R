#!/usr/bin/env Rscript
# ================================================================
# 组织特异性 APA 分析（含 Tau 指数，支持亚群分组）
# 输出基因列表供后续 GO 分析，不做富集
# 用法：Rscript specific_APA_with_subgroup.R -i PDUI.txt -g tissue_groups.txt -s subgroup_groups.txt -o output
# ================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(pheatmap)
  library(BiocParallel)
  library(ggplot2)
  library(reshape2)
})

# ------------------------- 参数解析 -------------------------
option_list <- list(
  make_option(c("-i", "--input"), type="character", help="PDUI 矩阵文件（tab分隔，第一列转录本ID）"),
  make_option(c("-g", "--group"), type="character", help="组织分组文件（有标题，第一列样本ID，第二列组织）"),
  make_option(c("-s", "--subgroup"), type="character", default=NULL, 
              help="亚群分组文件（可选，有标题，第一列样本ID，第二列亚群）"),
  make_option(c("-o", "--output"), type="character", default="APA_result", help="输出前缀"),
  make_option(c("-t", "--fdr"), type="double", default=0.05, help="FDR 阈值"),
  make_option(c("--tau"), type="double", default=0.5, help="Tau 特异性阈值（0~1）"),
  make_option(c("-c", "--cores"), type="integer", default=4, help="并行核心数"),
  make_option(c("--top_genes"), type="integer", default=2000, help="热图显示最大基因数")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$group)) stop("必须提供 -i 和 -g")

# ------------------------- 1. 读取数据 -------------------------
cat("读取 PDUI 矩阵...\n")
pdui <- fread(opt$input, header=TRUE, sep="\t", data.table=FALSE)
rownames(pdui) <- pdui[,1]
pdui <- as.matrix(pdui[,-1])

cat("读取组织分组...\n")
groups <- fread(opt$group, header=TRUE, sep="\t", data.table=FALSE)
sample_id <- groups[,1]
tissue <- groups[,2]
names(tissue) <- sample_id

# 若有亚群文件
if (!is.null(opt$subgroup)) {
  cat("读取亚群分组...\n")
  subg <- fread(opt$subgroup, header=TRUE, sep="\t", data.table=FALSE)
  sub_sample <- subg[,1]
  subgroup <- subg[,2]
  names(subgroup) <- sub_sample
  # 确保顺序与矩阵一致
  common <- intersect(colnames(pdui), intersect(names(tissue), names(subgroup)))
  if (length(common) == 0) stop("样本ID在三个文件中完全不匹配")
  pdui <- pdui[, common, drop=FALSE]
  tissue <- tissue[common]
  subgroup <- subgroup[common]
  tissue_factor <- factor(tissue)
  subgroup_factor <- factor(subgroup)
  cat(sprintf("样本: %d, 组织: %d, 亚群: %d\n", ncol(pdui), nlevels(tissue_factor), nlevels(subgroup_factor)))
} else {
  common <- intersect(colnames(pdui), names(tissue))
  pdui <- pdui[, common, drop=FALSE]
  tissue <- tissue[common]
  tissue_factor <- factor(tissue)
  subgroup_factor <- NULL
  cat(sprintf("样本: %d, 组织: %d\n", ncol(pdui), nlevels(tissue_factor)))
}

# ------------------------- 2. Kruskal 检验 (并行) -------------------------
cat("进行 Kruskal-Wallis 检验...\n")
BPPARAM <- MulticoreParam(workers=opt$cores, progressbar=TRUE)
p_values <- bplapply(seq_len(nrow(pdui)), function(i) {
  row_vals <- pdui[i, ]
  valid <- !is.na(row_vals)
  if (sum(valid) < 3) return(NA)
  x <- row_vals[valid]
  g <- tissue_factor[valid]
  if (any(table(g) < 2)) return(NA)
  tryCatch(kruskal.test(x, g)$p.value, error=function(e) NA)
}, BPPARAM=BPPARAM)
p_values <- unlist(p_values)

padj <- rep(NA, length(p_values))
padj[!is.na(p_values)] <- p.adjust(p_values[!is.na(p_values)], method="fdr")
sig_idx <- which(padj < opt$fdr & !is.na(padj))
cat(sprintf("显著转录本: %d\n", length(sig_idx)))
if (length(sig_idx) == 0) quit(save="no")

sig_data <- pdui[sig_idx, , drop=FALSE]

# ------------------------- 3. 计算组织均值矩阵 (忽略亚群) -------------------------
cat("计算组织均值矩阵...\n")
tissue_levels <- levels(tissue_factor)
agg_list <- lapply(tissue_levels, function(lev) {
  cols <- which(tissue_factor == lev)
  if (length(cols) == 1) sig_data[, cols, drop=FALSE] else rowMeans(sig_data[, cols], na.rm=TRUE)
})
agg_mat <- do.call(cbind, agg_list)
colnames(agg_mat) <- tissue_levels
agg_mat <- agg_mat[!apply(agg_mat, 1, function(x) all(is.na(x))), , drop=FALSE]

# ------------------------- 4. 计算总体 Tau 和 Dominant -------------------------
calc_tau <- function(x) {
  if (all(is.na(x))) return(NA)
  x_max <- max(x, na.rm=TRUE)
  if (x_max == 0) return(0)
  x_norm <- x / x_max
  sum(1 - x_norm, na.rm=TRUE) / (length(x) - 1)
}

get_dominant <- function(x) {
  x_sorted <- sort(x, decreasing=TRUE, na.last=NA)
  if (length(x_sorted) < 2) return(data.frame(dom=NA, delta=NA, max_val=NA))
  max_val <- x_sorted[1]
  second_val <- x_sorted[2]
  dom_tissue <- names(x_sorted)[1]
  return(data.frame(dom=dom_tissue, delta=max_val - second_val, max_val=max_val))
}

tau_vec <- apply(agg_mat, 1, calc_tau)
dom_info <- do.call(rbind, apply(agg_mat, 1, get_dominant))
result_table <- data.frame(
  Transcript = rownames(agg_mat),
  Dominant_Tissue = dom_info$dom,
  Max_PDUI = dom_info$max_val,
  Delta = dom_info$delta,
  Tau = tau_vec,
  P_value = p_values[match(rownames(agg_mat), rownames(pdui))],
  FDR = padj[match(rownames(agg_mat), rownames(pdui))]
)
result_table <- result_table[order(result_table$Tau, decreasing=TRUE), ]
write.table(result_table, file=paste0(opt$output, "_Tau_results.txt"), 
            sep="\t", row.names=FALSE, quote=FALSE)

# ------------------------- 5. 输出基因列表（忽略亚群） -------------------------
high_tau_idx <- which(result_table$Tau >= opt$tau & !is.na(result_table$Tau))
high_tau_genes <- result_table[high_tau_idx, ]
cat(sprintf("Tau >= %.2f 的高特异性基因 (忽略亚群): %d\n", opt$tau, nrow(high_tau_genes)))

# 按组织输出列表
for (tiss in unique(high_tau_genes$Dominant_Tissue)) {
  sub <- high_tau_genes[high_tau_genes$Dominant_Tissue == tiss, "Transcript"]
  write.table(sub, file=paste0(opt$output, "_", tiss, "_specific_genes.txt"), 
              row.names=FALSE, col.names=FALSE, quote=FALSE)
}
cat("组织特异性基因列表已输出（忽略亚群）\n")

# 柱状图
tissue_counts <- table(high_tau_genes$Dominant_Tissue)
df_bar <- data.frame(Tissue = names(tissue_counts), Count = as.numeric(tissue_counts))
p_bar <- ggplot(df_bar, aes(x=reorder(Tissue, -Count), y=Count, fill=Tissue)) +
  geom_bar(stat="identity") + geom_text(aes(label=Count), vjust=-0.3) +
  labs(title=paste0("Tissue-specific genes (Tau >= ", opt$tau, ")"), x="Tissue", y="Count") +
  theme_minimal() + theme(legend.position="none")
ggsave(paste0(opt$output, "_barplot.pdf"), p_bar, width=8, height=6)

# ------------------------- 6. 若提供亚群，计算亚群内组织特异性 -------------------------
if (!is.null(opt$subgroup)) {
  cat("计算各亚群内的组织特异性...\n")
  sub_levels <- levels(subgroup_factor)
  # 存储每个亚群的结果
  all_sub_results <- list()
  conserved_list <- list()  # 用于记录在所有亚群中dominant一致的基因
  
  for (sub in sub_levels) {
    sub_cols <- which(subgroup_factor == sub)
    sub_tissue <- tissue_factor[sub_cols]
    sub_data <- sig_data[, sub_cols, drop=FALSE]
    # 计算该亚群内的组织均值
    sub_agg <- lapply(tissue_levels, function(lev) {
      cols <- which(sub_tissue == lev)
      if (length(cols) == 0) return(rep(NA, nrow(sub_data)))
      if (length(cols) == 1) sub_data[, cols, drop=FALSE] else rowMeans(sub_data[, cols], na.rm=TRUE)
    })
    sub_agg <- do.call(cbind, sub_agg)
    colnames(sub_agg) <- tissue_levels
    # 移除全NA行
    sub_agg <- sub_agg[!apply(sub_agg, 1, function(x) all(is.na(x))), , drop=FALSE]
    # 计算Tau和Dominant
    tau_sub <- apply(sub_agg, 1, calc_tau)
    dom_sub <- do.call(rbind, apply(sub_agg, 1, get_dominant))
    res_sub <- data.frame(
      Transcript = rownames(sub_agg),
      Dominant_Tissue = dom_sub$dom,
      Tau = tau_sub,
      stringsAsFactors=FALSE
    )
    all_sub_results[[sub]] <- res_sub
    # 筛选该亚群内高特异性基因
    high_sub <- res_sub[which(res_sub$Tau >= opt$tau & !is.na(res_sub$Tau)), ]
    if (nrow(high_sub) > 0) {
      for (tiss in unique(high_sub$Dominant_Tissue)) {
        genes <- high_sub[high_sub$Dominant_Tissue == tiss, "Transcript"]
        fname <- paste0(opt$output, "_", tiss, "_", sub, "_specific_genes.txt")
        write.table(genes, file=fname, row.names=FALSE, col.names=FALSE, quote=FALSE)
      }
    }
    # 记录该亚群每个基因的dominant组织
    conserved_list[[sub]] <- setNames(res_sub$Dominant_Tissue, res_sub$Transcript)
  }
  
  # 找出在所有亚群中 dominant 组织一致的基因
  all_genes <- unique(unlist(lapply(conserved_list, names)))
  conserved_genes <- c()
  for (gene in all_genes) {
    doms <- sapply(conserved_list, function(x) x[gene])
    if (!any(is.na(doms)) && length(unique(doms)) == 1) {
      # 且该基因在总体中也满足Tau阈值（可选）
      if (gene %in% result_table$Transcript && result_table[result_table$Transcript == gene, "Tau"] >= opt$tau) {
        conserved_genes <- c(conserved_genes, gene)
      }
    }
  }
  if (length(conserved_genes) > 0) {
    write.table(conserved_genes, file=paste0(opt$output, "_conserved_across_subgroups.txt"),
                row.names=FALSE, col.names=FALSE, quote=FALSE)
    cat(sprintf("保守组织特异性基因（所有亚群一致）: %d\n", length(conserved_genes)))
  } else {
    cat("没有发现跨亚群保守的组织特异性基因\n")
  }
}

# ------------------------- 7. 热图（组织均值） -------------------------
cat("绘制热图（基于组织均值）...\n")
mat_heat <- agg_mat[rownames(agg_mat) %in% high_tau_genes$Transcript, , drop=FALSE]
if (nrow(mat_heat) > opt$top_genes) {
  keep <- high_tau_genes$Transcript[1:opt$top_genes]
  mat_heat <- mat_heat[keep, , drop=FALSE]
}
mat_scaled <- t(scale(t(mat_heat)))
mat_scaled[mat_scaled > 3] <- 3
mat_scaled[mat_scaled < -3] <- -3

annotation_col <- data.frame(Tissue = colnames(mat_scaled))
rownames(annotation_col) <- colnames(mat_scaled)
tissue_colors <- rainbow(nlevels(tissue_factor))
names(tissue_colors) <- tissue_levels

pdf(paste0(opt$output, "_heatmap.pdf"), width=10, height=12)
pheatmap(mat_scaled,
         main = paste0("Tissue-specific APA (Tau >= ", opt$tau, ")"),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         annotation_colors = list(Tissue = tissue_colors),
         show_rownames = FALSE,
         show_colnames = TRUE,
         color = colorRampPalette(c("navy", "white", "red"))(100),
         scale = "none")
dev.off()
cat("热图已保存\n")

cat("\n===== 全部完成 =====\n")
cat("输出文件:\n")
cat("  - ", opt$output, "_Tau_results.txt (完整结果表)\n")
cat("  - ", opt$output, "_barplot.pdf (组织特异基因数量柱状图)\n")
cat("  - ", opt$output, "_heatmap.pdf\n")
cat("  - *_specific_genes.txt (按组织/组织-亚群的特异基因列表)\n")
if (!is.null(opt$subgroup)) {
  cat("  - ", opt$output, "_conserved_across_subgroups.txt (所有亚群一致的组织特异基因)\n")
}
cat("这些基因列表可直接用于后续 GO 富集分析。\n")
