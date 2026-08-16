#!/usr/bin/env Rscript
# ================================================================
# 生成干净的 ECDF、差分ECDF 和箱线图（0~1全范围）
# 用法：Rscript plot_clean.R -i PDUI.txt -p subpop.txt -t tissue.txt -o out
# ================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(ggplot2)
  library(reshape2)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="PDUI矩阵"),
  make_option(c("-p", "--pop"), type="character", help="亚群分组文件"),
  make_option(c("-t", "--tiss"), type="character", help="组织分组文件"),
  make_option(c("-o", "--output"), type="character", default="Tau_Comparison", help="输出前缀")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input) || is.null(opt$pop) || is.null(opt$tiss)) stop("缺少 -i, -p, -t")

# ------------------------- 1. 读取与计算 -------------------------
cat("读取数据...\n")
pdui <- fread(opt$input, header=TRUE, sep="\t", data.table=FALSE)
rownames(pdui) <- pdui[,1]
pdui <- as.matrix(pdui[,-1])

pop_df <- fread(opt$pop, header=TRUE, sep="\t", data.table=FALSE)
pop_vec <- setNames(pop_df[,2], pop_df[,1])
tiss_df <- fread(opt$tiss, header=TRUE, sep="\t", data.table=FALSE)
tiss_vec <- setNames(tiss_df[,2], tiss_df[,1])

common <- intersect(colnames(pdui), intersect(names(pop_vec), names(tiss_vec)))
pdui <- pdui[, common, drop=FALSE]
pop_vec <- pop_vec[common]; tiss_vec <- tiss_vec[common]
pop_factor <- factor(pop_vec); tiss_factor <- factor(tiss_vec)
cat(sprintf("样本: %d, 亚群: %d, 组织: %d\n", ncol(pdui), nlevels(pop_factor), nlevels(tiss_factor)))

calc_tau <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(NA)
  x_max <- max(x)
  if (x_max == 0) return(0)
  sum(1 - x / x_max) / (length(x) - 1)
}

# 组织均值
tiss_agg <- sapply(levels(tiss_factor), function(lev) {
  cols <- which(tiss_factor == lev)
  if (length(cols) == 1) pdui[, cols] else rowMeans(pdui[, cols, drop=FALSE], na.rm=TRUE)
})
tiss_agg <- tiss_agg[!apply(tiss_agg, 1, function(x) all(is.na(x))), , drop=FALSE]
tiss_tau <- apply(tiss_agg, 1, calc_tau)

# 亚群均值
pop_agg <- sapply(levels(pop_factor), function(lev) {
  cols <- which(pop_factor == lev)
  if (length(cols) == 1) pdui[, cols] else rowMeans(pdui[, cols, drop=FALSE], na.rm=TRUE)
})
pop_agg <- pop_agg[!apply(pop_agg, 1, function(x) all(is.na(x))), , drop=FALSE]
pop_tau <- apply(pop_agg, 1, calc_tau)

# 配对数据
all_genes <- union(names(tiss_tau), names(pop_tau))
paired <- data.frame(
  Transcript = all_genes,
  TissueTau = tiss_tau[match(all_genes, names(tiss_tau))],
  PopulationTau = pop_tau[match(all_genes, names(pop_tau))]
)
paired <- paired[!is.na(paired$TissueTau) & !is.na(paired$PopulationTau), ]
med_tiss <- median(paired$TissueTau); med_pop <- median(paired$PopulationTau)
wt <- wilcox.test(paired$TissueTau, paired$PopulationTau, paired=TRUE)
p_val <- wt$p.value
cat(sprintf("组织Tau中位数: %.4f, 亚群Tau中位数: %.4f\n", med_tiss, med_pop))
cat(sprintf("配对Wilcoxon P = %g\n", p_val))

# 长格式
df_long <- rbind(
  data.frame(Tau = paired$TissueTau, Group = "Tissue"),
  data.frame(Tau = paired$PopulationTau, Group = "Population")
)
df_long$Group <- factor(df_long$Group, levels = c("Tissue", "Population"))

# ------------------------- 2. ECDF图（0-1全范围，无网格，白色背景） -------------------------
p_ecdf <- ggplot(df_long, aes(x = Tau, color = Group)) +
  stat_ecdf(size = 1) +
  scale_color_manual(values = c("Tissue" = "#D55E00", "Population" = "#0072B2")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(title = "Cumulative distribution of Tau",
       x = "Tau", y = "ECDF") +
  theme_classic(base_size = 13) +            # 经典白底无网格
  theme(
    panel.grid = element_blank(),            # 确保无网格
    legend.position = c(0.8, 0.2),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", color = "black", size = 0.3)
  ) +
  annotate("text", x = 0.7, y = 0.1,
           label = paste0("Wilcoxon P = ", format(p_val, scientific = TRUE, digits = 2)),
           hjust = 0, size = 3.5)
ggsave(paste0(opt$output, "_ECDF.pdf"), p_ecdf, width = 7, height = 5)

# ------------------------- 3. 差分ECDF图（0-1全范围） -------------------------
ecdf_tiss <- ecdf(paired$TissueTau)
ecdf_pop <- ecdf(paired$PopulationTau)
tau_seq <- seq(0, 1, length.out = 500)
diff_ecdf <- ecdf_tiss(tau_seq) - ecdf_pop(tau_seq)
diff_df <- data.frame(Tau = tau_seq, Diff = diff_ecdf)

p_diff <- ggplot(diff_df, aes(x = Tau, y = Diff)) +
  geom_line(size = 1, color = "purple") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.6) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(title = "Difference in ECDF (Tissue - Population)",
       x = "Tau", y = "ΔECDF") +
  theme_classic(base_size = 13) +
  theme(panel.grid = element_blank()) +
  annotate("text", x = 0.7, y = 0.02,
           label = "Positive: Tissue has more low-Tau transcripts",
           hjust = 0, size = 3.5, color = "grey30") +
  annotate("text", x = 0.7, y = -0.02,
           label = "Negative: Tissue has more high-Tau transcripts",
           hjust = 0, size = 3.5, color = "grey30")
ggsave(paste0(opt$output, "_DiffECDF.pdf"), p_diff, width = 7, height = 5)

# ------------------------- 4. 箱线图（保留全部数据，无网格，白底） -------------------------
p_box <- ggplot(df_long, aes(x = Group, y = Tau, fill = Group)) +
  geom_boxplot(width = 0.5, alpha = 0.7, outlier.size = 0.3, outlier.shape = 16) +
  scale_fill_manual(values = c("Tissue" = "#D55E00", "Population" = "#0072B2")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(title = "Tau distribution", x = "", y = "Tau") +
  theme_classic(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  annotate("text", x = 1.5, y = 1.02,
           label = paste0("P = ", format(p_val, scientific = TRUE, digits = 2)),
           hjust = 0.5, size = 4)
ggsave(paste0(opt$output, "_Boxplot.pdf"), p_box, width = 5, height = 6)

cat("\n===== 完成 =====\n")
cat("输出文件:\n  - ", opt$output, "_ECDF.pdf\n  - ", opt$output, "_DiffECDF.pdf\n  - ", opt$output, "_Boxplot.pdf\n")
