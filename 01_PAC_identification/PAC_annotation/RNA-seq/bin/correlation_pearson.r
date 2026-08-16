# 生物重复一致性评估脚本
# 输入：
#   - 表达矩阵文件（行=基因，列=样本）
#   - 分组信息文件（列：sample, group）
# 输出：
#   - 组内相关性分布图
#   - 组内 vs 组间相关性比较箱线图
#   - 相关性矩阵热图

# ---------------------------
# 1. 加载必要的包
# ---------------------------
library(corrplot)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# ---------------------------
# 2. 数据加载与预处理
# ---------------------------
# 读取表达矩阵（假设为log2归一化后的数据）
expr_matrix <- as.matrix(read.delim("/share2/pub/xingsl/xingsl/data_download/MZAPA/ZY/post_analysis/PAC_CS_RNAseq_feature/SRA055066/filter01/SRA055066.hq.PAC.TPM", 
                                   row.names = 1, 
                                   check.names = FALSE))

# 读取分组信息
group_info <- read.delim("condition.txt", 
                         stringsAsFactors = FALSE)

# 检查样本一致性
if (!all(colnames(expr_matrix) %in% group_info$sample)) {
  stop("错误：表达矩阵中的样本名与分组文件不匹配！")
}

# 确保分组顺序与表达矩阵列名一致
group_info <- group_info[match(colnames(expr_matrix), group_info$sample), ]

# ---------------------------
# 3. 计算相关性矩阵
# ---------------------------
# Pearson相关性矩阵（214x214）
cor_matrix <- cor(expr_matrix, method = "pearson")

# ---------------------------
# 4. 提取组内重复相关性
# ---------------------------
# 假设分组列中的组名为Group1, Group2,...，重复标记为Rep1/Rep2
group_cors <- sapply(unique(group_info$group), function(g) {
  samples <- group_info$sample[group_info$group == g]
  if (length(samples) != 2) {
    warning(sprintf("组 %s 的样本数不是2个，跳过", g))
    return(NA)
  }
  cor_matrix[samples[1], samples[2]]
})

# 移除NA值（无效组）
group_cors <- na.omit(group_cors)

# ---------------------------
# 5. 可视化分析
# ---------------------------
# 5.1 组内相关性分布直方图
pdf("within_group_correlation_hist.pdf", width = 6, height = 5)
hist(group_cors, breaks = 20, col = "lightblue", 
     main = "Distribution of Within-Group Correlations",
     xlab = "Pearson Correlation Coefficient", 
     ylab = "Frequency",
     xlim = c(0.5, 1))
abline(v = median(group_cors), col = "red", lty = 2, lwd = 2)
legend("topright", legend = sprintf("Median = %.3f", median(group_cors)), 
       col = "red", lty = 2, bty = "n")
dev.off()

# 5.2 组内 vs 随机组间相关性比较
set.seed(123)
random_pairs <- combn(sample(colnames(expr_matrix)), 2)
random_cors <- apply(random_pairs, 2, function(p) cor_matrix[p[1], p[2]])

comparison_df <- data.frame(
  Correlation = c(group_cors, random_cors),
  Type = rep(c("Within-Group", "Between-Groups"), 
             times = c(length(group_cors), length(random_cors)))
)

pdf("correlation_comparison_boxplot.pdf", width = 5, height = 6)
ggplot(comparison_df, aes(x = Type, y = Correlation, fill = Type)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.5, alpha = 0.3) +
  scale_fill_manual(values = c("lightgreen", "pink")) +
  labs(title = "Within-Group vs Between-Groups Correlation",
       y = "Pearson R") +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()

# 5.3 相关性矩阵热图（前50个样本，避免过度拥挤）
sample_subset <- head(colnames(expr_matrix), 50)
pdf("correlation_heatmap.pdf", width = 10, height = 9)
pheatmap(cor_matrix[sample_subset, sample_subset],
         color = colorRampPalette(brewer.pal(9, "Blues"))(100),
         clustering_method = "average",
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Sample Correlation Heatmap (Subset)")
dev.off()

# ---------------------------
# 6. 结果统计输出
# ---------------------------
# 生成统计摘要
stats_summary <- data.frame(
  Metric = c("Min", "Median", "Mean", "Max", "SD"),
  Value = c(min(group_cors), 
            median(group_cors),
            mean(group_cors), 
            max(group_cors),
            sd(group_cors))
)

# 保存结果
write.table(stats_summary, "correlation_stats.txt", sep = "\t", quote = FALSE)
write.table(data.frame(Group = names(group_cors), Correlation = group_cors),
            "group_correlations.txt", sep = "\t", row.names = FALSE)

