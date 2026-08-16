#!/usr/bin/env Rscript
# ================================================================
# 比较组织 vs 亚群特异性（Tau指数）
# 输出：密度分布图（多曲线）、四色散点图、Top ΔTau热图
# 用法：Rscript compare_tissue_vs_population_final.R \
#         -i PDUI_matrix.txt \
#         -p subgroup_groups.txt \
#         -t tissue_groups.txt \
#         -o output_prefix \
#         --topN 30
# ================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(ggplot2)
  library(reshape2)
  library(pheatmap)
  library(RColorBrewer)
  library(BiocParallel)   # 用于并行KW检验
})

# ------------------------- 参数解析 -------------------------
option_list <- list(
  make_option(c("-i", "--input"), type="character", help="PDUI矩阵（tab分隔，第一列为转录本ID）"),
  make_option(c("-p", "--pop"), type="character", help="亚群分组文件（有标题，第一列样本ID，第二列亚群）"),
  make_option(c("-t", "--tiss"), type="character", help="组织分组文件（有标题，第一列样本ID，第二列组织）"),
  make_option(c("-o", "--output"), type="character", default="Tissue_vs_Population", help="输出前缀"),
  make_option(c("--topN"), type="integer", default=30, help="热图展示的top ΔTau转录本数"),
  make_option(c("--tau_cut"), type="double", default=0.5, help="四象限阈值"),
  make_option(c("--fdr"), type="double", default=0.05, help="KW检验FDR阈值"),
  make_option(c("--cores"), type="integer", default=4, help="并行核心数")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$pop) || is.null(opt$tiss)) {
  stop("必须提供 -i, -p, -t 三个参数", call. = FALSE)
}

POP_NAME <- "Population"
TISS_NAME <- "Tissue"

# ------------------------- 1. 读取数据 -------------------------
cat("读取PDUI矩阵...\n")
pdui <- fread(opt$input, header=TRUE, sep="\t", data.table=FALSE)
rownames(pdui) <- pdui[,1]
pdui <- as.matrix(pdui[,-1])

cat("读取亚群分组...\n")
pop_df <- fread(opt$pop, header=TRUE, sep="\t", data.table=FALSE)
pop_vec <- setNames(pop_df[,2], pop_df[,1])

cat("读取组织分组...\n")
tiss_df <- fread(opt$tiss, header=TRUE, sep="\t", data.table=FALSE)
tiss_vec <- setNames(tiss_df[,2], tiss_df[,1])

common <- intersect(colnames(pdui), intersect(names(pop_vec), names(tiss_vec)))
if (length(common) == 0) stop("样本ID在三文件中完全不匹配")
pdui <- pdui[, common, drop=FALSE]
pop_vec <- pop_vec[common]
tiss_vec <- tiss_vec[common]
pop_factor <- factor(pop_vec)
tiss_factor <- factor(tiss_vec)
cat(sprintf("样本数: %d, 亚群数: %d, 组织数: %d\n", 
            ncol(pdui), nlevels(pop_factor), nlevels(tiss_factor)))

# ------------------------- 2. Kruskal-Wallis 检验（并行） -------------------------
# 分别对组织和亚群进行检验，得到FDR
cat("进行组织Kruskal-Wallis检验...\n")
BPPARAM <- MulticoreParam(workers = opt$cores, progressbar = TRUE)

kw_test <- function(row_vals, group_factor) {
  valid <- !is.na(row_vals)
  if (sum(valid) < 3) return(NA)
  x <- row_vals[valid]
  g <- group_factor[valid]
  if (any(table(g) < 2)) return(NA)
  tryCatch(kruskal.test(x, g)$p.value, error = function(e) NA)
}

p_tiss <- bplapply(seq_len(nrow(pdui)), function(i) kw_test(pdui[i, ], tiss_factor), BPPARAM = BPPARAM)
p_tiss <- unlist(p_tiss)
padj_tiss <- rep(NA, length(p_tiss))
padj_tiss[!is.na(p_tiss)] <- p.adjust(p_tiss[!is.na(p_tiss)], method = "fdr")
sig_tiss <- padj_tiss < opt$fdr & !is.na(padj_tiss)

cat("进行亚群Kruskal-Wallis检验...\n")
p_pop <- bplapply(seq_len(nrow(pdui)), function(i) kw_test(pdui[i, ], pop_factor), BPPARAM = BPPARAM)
p_pop <- unlist(p_pop)
padj_pop <- rep(NA, length(p_pop))
padj_pop[!is.na(p_pop)] <- p.adjust(p_pop[!is.na(p_pop)], method = "fdr")
sig_pop <- padj_pop < opt$fdr & !is.na(padj_pop)

cat(sprintf("组织显著转录本: %d, 亚群显著转录本: %d\n", sum(sig_tiss, na.rm=TRUE), sum(sig_pop, na.rm=TRUE)))

# ------------------------- 3. 计算组织/亚群均值矩阵及Tau -------------------------
calc_tau <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(NA)
  x_max <- max(x)
  if (x_max == 0) return(0)
  sum(1 - x / x_max) / (length(x) - 1)
}
get_dominant <- function(x) {
  if (all(is.na(x))) return(NA)
  m <- max(x, na.rm=TRUE)
  paste(names(x)[which(x == m)], collapse=";")
}

# 组织均值
tiss_levels <- levels(tiss_factor)
tiss_agg <- sapply(tiss_levels, function(lev) {
  cols <- which(tiss_factor == lev)
  if (length(cols) == 1) pdui[, cols] else rowMeans(pdui[, cols, drop=FALSE], na.rm=TRUE)
})
tiss_agg <- tiss_agg[!apply(tiss_agg, 1, function(x) all(is.na(x))), , drop=FALSE]
tiss_tau <- apply(tiss_agg, 1, calc_tau)
tiss_dom <- apply(tiss_agg, 1, get_dominant)

# 亚群均值
pop_levels <- levels(pop_factor)
pop_agg <- sapply(pop_levels, function(lev) {
  cols <- which(pop_factor == lev)
  if (length(cols) == 1) pdui[, cols] else rowMeans(pdui[, cols, drop=FALSE], na.rm=TRUE)
})
pop_agg <- pop_agg[!apply(pop_agg, 1, function(x) all(is.na(x))), , drop=FALSE]
pop_tau <- apply(pop_agg, 1, calc_tau)
pop_dom <- apply(pop_agg, 1, get_dominant)

# ------------------------- 4. 合并结果表 -------------------------
all_genes <- union(rownames(tiss_agg), rownames(pop_agg))
result_df <- data.frame(
  Transcript = all_genes,
  TissueTau = tiss_tau[match(all_genes, rownames(tiss_agg))],
  TissueDominant = tiss_dom[match(all_genes, rownames(tiss_agg))],
  PopulationTau = pop_tau[match(all_genes, rownames(pop_agg))],
  PopulationDominant = pop_dom[match(all_genes, rownames(pop_agg))],
  stringsAsFactors = FALSE
)
result_df$DeltaTau <- result_df$TissueTau - result_df$PopulationTau

# 加入KW显著标记
result_df$TissueSig <- sig_tiss[match(all_genes, rownames(pdui))]
result_df$PopulationSig <- sig_pop[match(all_genes, rownames(pdui))]

# 保存表格
write.table(result_df, file = paste0(opt$output, "_Tau_Comparison.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

# 用于绘图的数据（去除Tau为NA）
plot_df <- result_df[!is.na(result_df$TissueTau) & !is.na(result_df$PopulationTau), ]

# ------------------------- 5. 统计计算 -------------------------
epsilon <- 1e-6
above <- sum(plot_df$TissueTau > plot_df$PopulationTau + epsilon)
below <- sum(plot_df$PopulationTau > plot_df$TissueTau + epsilon)
equal <- nrow(plot_df) - above - below
cat(sprintf("组织效应更强: %d (%.1f%%), 亚群效应更强: %d (%.1f%%), 近似相等: %d (%.1f%%)\n",
            above, 100*above/nrow(plot_df), below, 100*below/nrow(plot_df), equal, 100*equal/nrow(plot_df)))

wt <- wilcox.test(plot_df$TissueTau, plot_df$PopulationTau, paired = TRUE)
cat(sprintf("配对Wilcoxon检验 P值: %g\n", wt$p.value))

# 四象限分类
tau_cut <- opt$tau_cut
plot_df$Quadrant <- with(plot_df,
  ifelse(TissueTau >= tau_cut & PopulationTau < tau_cut, "Tissue-specific",
  ifelse(TissueTau >= tau_cut & PopulationTau >= tau_cut, "Both",
  ifelse(TissueTau < tau_cut & PopulationTau < tau_cut, "None", "Population-specific")))
)
quad_counts <- table(plot_df$Quadrant)
cat("四象限计数 (阈值=", tau_cut, "):\n", sep="")
print(quad_counts)

# ------------------------- 6. 密度分布图（多曲线） -------------------------
cat("绘制密度分布图...\n")

# 定义子集
plot_df$TissueSig <- as.logical(plot_df$TissueSig)
plot_df$PopulationSig <- as.logical(plot_df$PopulationSig)

# 创建长格式数据框用于ggplot：每一行是一个转录本的TissueTau或PopulationTau，并带有分类标签
# 我们为每条曲线准备一个数据框，然后rbind
curves <- list()

# 全部组织Tau
curves$all_tiss <- data.frame(
  Tau = plot_df$TissueTau,
  Group = "All transcripts (Tissue)",
  Type = "All",
  Category = "All"
)
# 全部亚群Tau
curves$all_pop <- data.frame(
  Tau = plot_df$PopulationTau,
  Group = "All transcripts (Population)",
  Type = "All",
  Category = "All"
)

# 组织显著的组织Tau
idx <- which(plot_df$TissueSig == TRUE)
if (length(idx) > 0) {
  curves$sig_tiss <- data.frame(
    Tau = plot_df$TissueTau[idx],
    Group = "Tissue significant (FDR<0.05)",
    Type = "TissueSig",
    Category = "TissueSig"
  )
}
# 组织显著且TissueTau>=0.5
idx2 <- which(plot_df$TissueSig == TRUE & plot_df$TissueTau >= tau_cut)
if (length(idx2) > 0) {
  curves$sig_tiss_high <- data.frame(
    Tau = plot_df$TissueTau[idx2],
    Group = "Tissue specific (Tau≥0.5 & FDR<0.05)",
    Type = "TissueHigh",
    Category = "TissueHigh"
  )
}
# 亚群显著的亚群Tau
idx3 <- which(plot_df$PopulationSig == TRUE)
if (length(idx3) > 0) {
  curves$sig_pop <- data.frame(
    Tau = plot_df$PopulationTau[idx3],
    Group = "Population significant (FDR<0.05)",
    Type = "PopulationSig",
    Category = "PopulationSig"
  )
}
# 亚群显著且PopulationTau>=0.5
idx4 <- which(plot_df$PopulationSig == TRUE & plot_df$PopulationTau >= tau_cut)
if (length(idx4) > 0) {
  curves$sig_pop_high <- data.frame(
    Tau = plot_df$PopulationTau[idx4],
    Group = "Population specific (Tau≥0.5 & FDR<0.05)",
    Type = "PopulationHigh",
    Category = "PopulationHigh"
  )
}

# 合并所有曲线
curve_df <- do.call(rbind, curves)
# 定义线型、颜色映射
line_types <- c("All" = "solid", "TissueSig" = "solid", "TissueHigh" = "solid",
                "PopulationSig" = "dashed", "PopulationHigh" = "dashed")
line_colors <- c("All" = "grey30", "TissueSig" = "blue", "TissueHigh" = "red",
                 "PopulationSig" = "#66CC66", "PopulationHigh" = "#006600")
# 注意：All有两个（组织全部和亚群全部），我们通过Group区分，但使用相同的颜色和线型？为了区分，我们给All组织用实线灰色，All亚群用虚线灰色。
# 但上面我们将All组织Type="All"且Category="All"，All亚群也是Type="All"，它们会共用线型，我们需要单独处理。
# 更稳健：为每个曲线单独指定linetype和color。
curve_df$LineType <- ifelse(curve_df$Group == "All transcripts (Population)", "dashed", "solid")
curve_df$LineColor <- ifelse(curve_df$Group == "All transcripts (Population)", "grey50", 
                             ifelse(curve_df$Group == "All transcripts (Tissue)", "grey30",
                                    ifelse(curve_df$Type == "TissueSig", "blue",
                                           ifelse(curve_df$Type == "TissueHigh", "red",
                                                  ifelse(curve_df$Type == "PopulationSig", "#66CC66",
                                                         ifelse(curve_df$Type == "PopulationHigh", "#006600", "black"))))))

# 绘制
p_density <- ggplot(curve_df, aes(x = Tau, color = Group, linetype = Group)) +
  geom_density(size = 1.2, adjust = 1.5) +
  scale_color_manual(values = c("All transcripts (Tissue)" = "grey30",
                                "All transcripts (Population)" = "grey50",
                                "Tissue significant (FDR<0.05)" = "blue",
                                "Tissue specific (Tau≥0.5 & FDR<0.05)" = "red",
                                "Population significant (FDR<0.05)" = "#66CC66",
                                "Population specific (Tau≥0.5 & FDR<0.05)" = "#006600")) +
  scale_linetype_manual(values = c("All transcripts (Tissue)" = "solid",
                                   "All transcripts (Population)" = "dashed",
                                   "Tissue significant (FDR<0.05)" = "solid",
                                   "Tissue specific (Tau≥0.5 & FDR<0.05)" = "solid",
                                   "Population significant (FDR<0.05)" = "dashed",
                                   "Population specific (Tau≥0.5 & FDR<0.05)" = "dashed")) +
  labs(
    title = "Tissue vs Population APA Specificity (Tau)",
    x = "Tau",
    y = "Density"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.8),
    legend.position = c(0.8, 0.8),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", color = "black", size = 0.5)
  ) +
  xlim(0, 1) +
  # 添加Wilcoxon p值
  annotate("text", x = 0.95, y = Inf, 
           label = paste0("Paired Wilcoxon\nP = ", format(wt$p.value, scientific = TRUE, digits = 3)),
           hjust = 1, vjust = 1.5, size = 4, color = "black")

# 保存
ggsave(paste0(opt$output, "_Tau_Density.pdf"), p_density, width = 8, height = 6)
cat("密度图已保存\n")

# ------------------------- 7. 四色散点图 -------------------------
cat("绘制四色散点图...\n")
# 定义象限颜色
quad_colors <- c("Tissue-specific" = "red", 
                 "Both" = "purple", 
                 "None" = "grey70", 
                 "Population-specific" = "blue")
p_scatter <- ggplot(plot_df, aes(x = PopulationTau, y = TissueTau, color = Quadrant)) +
  geom_point(alpha = 0.4, size = 0.8) +
  scale_color_manual(values = quad_colors) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = tau_cut, linetype = "dotted", color = "grey30", size = 0.8) +
  geom_hline(yintercept = tau_cut, linetype = "dotted", color = "grey30", size = 0.8) +
  coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1)) +
  labs(
    title = "Quadrant Classification",
    x = "Population Tau",
    y = "Tissue Tau",
    color = "Category"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.8),
    legend.position = "right"
  ) +
  # 添加象限计数文本（在角落显示）
  annotate("text", x = 0.1, y = 0.9, label = paste("Tissue-specific:", quad_counts["Tissue-specific"]), color = "red", hjust = 0, size = 4) +
  annotate("text", x = 0.1, y = 0.8, label = paste("Both:", quad_counts["Both"]), color = "purple", hjust = 0, size = 4) +
  annotate("text", x = 0.1, y = 0.7, label = paste("None:", quad_counts["None"]), color = "grey30", hjust = 0, size = 4) +
  annotate("text", x = 0.1, y = 0.6, label = paste("Population-specific:", quad_counts["Population-specific"]), color = "blue", hjust = 0, size = 4)

ggsave(paste0(opt$output, "_Quadrant_Scatter.pdf"), p_scatter, width = 8, height = 7)
cat("四色散点图已保存\n")

# ------------------------- 8. Top ΔTau 热图 -------------------------
cat("绘制Top ΔTau热图...\n")
topN <- min(opt$topN, nrow(plot_df))
top_genes <- plot_df[order(plot_df$DeltaTau, decreasing = TRUE), "Transcript"][1:topN]
sub_pdui <- pdui[rownames(pdui) %in% top_genes, , drop = FALSE]

if (nrow(sub_pdui) > 0) {
  # 计算亚群均值
  pop_sub <- sapply(pop_levels, function(lev) {
    cols <- which(pop_factor == lev)
    if (length(cols) == 1) sub_pdui[, cols] else rowMeans(sub_pdui[, cols, drop = FALSE], na.rm = TRUE)
  })
  rownames(pop_sub) <- rownames(sub_pdui)
  # 组织均值
  tiss_sub <- sapply(tiss_levels, function(lev) {
    cols <- which(tiss_factor == lev)
    if (length(cols) == 1) sub_pdui[, cols] else rowMeans(sub_pdui[, cols, drop = FALSE], na.rm = TRUE)
  })
  rownames(tiss_sub) <- rownames(sub_pdui)
  
  combined <- cbind(pop_sub, tiss_sub)
  colnames(combined) <- c(paste0("Pop_", pop_levels), paste0("Tiss_", tiss_levels))
  
  # Z-score
  mat_scaled <- t(scale(t(combined)))
  mat_scaled[mat_scaled > 3] <- 3
  mat_scaled[mat_scaled < -3] <- -3
  
  ann_col <- data.frame(Group = factor(c(rep("Population", length(pop_levels)), rep("Tissue", length(tiss_levels)))))
  rownames(ann_col) <- colnames(combined)
  ann_colors <- list(Group = c(Population = "#1b9e77", Tissue = "#d95f02"))
  
  pdf(paste0(opt$output, "_Top", topN, "_Heatmap.pdf"), width = 12, height = max(6, topN * 0.3))
  pheatmap(mat_scaled,
           main = paste0("Top ", topN, " ΔTau (Tissue - Population)"),
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           annotation_col = ann_col,
           annotation_colors = ann_colors,
           show_rownames = TRUE,
           show_colnames = TRUE,
           fontsize_row = 6,
           fontsize_col = 8,
           color = colorRampPalette(c("navy", "white", "red"))(100),
           scale = "none")
  dev.off()
  cat("热图已保存\n")
}

cat("\n===== 分析完成 =====\n")
cat("输出文件:\n")
cat("  - ", opt$output, "_Tau_Comparison.txt (完整结果表)\n")
cat("  - ", opt$output, "_Tau_Density.pdf (密度分布图)\n")
cat("  - ", opt$output, "_Quadrant_Scatter.pdf (四色散点图)\n")
cat("  - ", opt$output, "_Top", topN, "_Heatmap.pdf (热图)\n")
