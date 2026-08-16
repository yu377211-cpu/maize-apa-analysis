# ============================================================================
# APA-Expression 相关性火山图（简化版）
# 专门突出 251 vs 66 的核心发现
# ============================================================================

# 加载必要的包
library(ggplot2)

cat("生成APA相关性火山图...\n")

# 读取全局相关性数据
global_corr <- read.table("PDUI_PAC_global_correlation.txt", 
                          header = TRUE, sep = "\t", 
                          stringsAsFactors = FALSE)

# 数据清理：处理NA值
global_corr$FDR[is.na(global_corr$FDR)] <- 1
global_corr$Correlation[is.na(global_corr$Correlation)] <- 0

# 计算统计
strong_neg <- sum(global_corr$Correlation < -0.5 & global_corr$FDR < 0.05, na.rm = TRUE)
strong_pos <- sum(global_corr$Correlation > 0.5 & global_corr$FDR < 0.05, na.rm = TRUE)
ratio <- round(strong_neg/strong_pos, 1)

cat("统计结果:\n")
cat("强负相关 (r < -0.5):", strong_neg, "\n")
cat("强正相关 (r > 0.5):", strong_pos, "\n")
cat("比例 (负:正):", strong_neg, ":", strong_pos, sprintf("(%.1f:1)\n", ratio))

# 添加分类信息
global_corr$Category <- "Not significant"
global_corr$Category[global_corr$FDR < 0.05 & global_corr$Correlation > 0.5] <- "Strong Positive"
global_corr$Category[global_corr$FDR < 0.05 & global_corr$Correlation < -0.5] <- "Strong Negative"
global_corr$Category[global_corr$FDR < 0.05 & 
                     global_corr$Correlation >= -0.5 & 
                     global_corr$Correlation <= 0.5] <- "Weak Correlation"

# 计算-log10(FDR)
global_corr$neg_log10_FDR <- ifelse(global_corr$FDR == 0, 300, 
                                   -log10(global_corr$FDR + 1e-300))

# 设置颜色
colors <- c("Strong Negative" = "#E41A1C",  # 红色
            "Strong Positive" = "#377EB8",  # 蓝色
            "Weak Correlation" = "#4DAF4A", # 绿色
            "Not significant" = "#999999")  # 灰色

# 创建火山图
p <- ggplot(global_corr, aes(x = Correlation, y = neg_log10_FDR)) +
  geom_point(aes(color = Category), alpha = 0.6, size = 1.5) +
  
  # 阈值线
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray50", size = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50", size = 0.5) +
  
  # 颜色设置
  scale_color_manual(values = colors) +
  
  # 坐标轴
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) +
  
  # 标签
  labs(
    title = "APA-Gene Expression Correlations in Maize",
    subtitle = paste("Strong Negative:", strong_neg, " | Strong Positive:", strong_pos, 
                    " | Ratio:", ratio, ":1"),
    x = "Spearman Correlation Coefficient (r)",
    y = "-log₁₀(FDR)",
    color = "Correlation Type"
  ) +
  
  # 主题
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "red"),
    legend.position = "right"
  ) +
  
  # 添加统计标注
  annotate("text", x = -0.75, y = max(global_corr$neg_log10_FDR, na.rm = TRUE) * 0.9, 
           label = paste("Strong Negative\nn =", strong_neg),
           color = "#E41A1C", size = 4.5, fontface = "bold") +
  
  annotate("text", x = 0.75, y = max(global_corr$neg_log10_FDR, na.rm = TRUE) * 0.9, 
           label = paste("Strong Positive\nn =", strong_pos),
           color = "#377EB8", size = 4.5, fontface = "bold")

# 保存图片
ggsave("APA_Correlation_Volcano_Simple.png", p, 
       width = 9, height = 7, dpi = 300, bg = "white")
ggsave("APA_Correlation_Volcano_Simple.pdf", p, 
       width = 9, height = 7, bg = "white")

cat("\n火山图已保存:\n")
cat("• APA_Correlation_Volcano_Simple.png\n")
cat("• APA_Correlation_Volcano_Simple.pdf\n")
cat("\n完成！\n")
