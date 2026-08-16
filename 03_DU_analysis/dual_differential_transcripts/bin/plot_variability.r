#!/usr/bin/env Rscript
# ==============================================================================
# Script: plot_variability_comparison.R
# Purpose: Create comparison plots for Tissue vs Subgroup variability
# Input: PDUI_PAC_tissue_correlation.txt and PDUI_PAC_subgroup_correlation.txt
# Output: Multiple comparison plots in PNG and PDF formats
# ==============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)

# 设置参数
args <- commandArgs(trailingOnly = TRUE)

# 如果没有提供参数，使用默认文件名
if (length(args) == 0) {
  tissue_file <- "PDUI_PAC_tissue_correlation.txt"
  subgroup_file <- "PDUI_PAC_subgroup_correlation.txt"
  output_dir <- "."
  cat("Using default file names:\n")
  cat("  Tissue file:", tissue_file, "\n")
  cat("  Subgroup file:", subgroup_file, "\n")
  cat("  Output directory:", output_dir, "\n\n")
} else if (length(args) == 2) {
  tissue_file <- args[1]
  subgroup_file <- args[2]
  output_dir <- "."
  cat("Using provided file names:\n")
  cat("  Tissue file:", tissue_file, "\n")
  cat("  Subgroup file:", subgroup_file, "\n")
  cat("  Output directory:", output_dir, "\n\n")
} else if (length(args) == 3) {
  tissue_file <- args[1]
  subgroup_file <- args[2]
  output_dir <- args[3]
  cat("Using provided parameters:\n")
  cat("  Tissue file:", tissue_file, "\n")
  cat("  Subgroup file:", subgroup_file, "\n")
  cat("  Output directory:", output_dir, "\n\n")
} else {
  stop("Usage: Rscript plot_variability_comparison.R [tissue_file] [subgroup_file] [output_dir]\n",
       "  If no arguments provided, uses default file names.")
}

# 检查文件是否存在
if (!file.exists(tissue_file)) {
  stop("Tissue file not found: ", tissue_file)
}
if (!file.exists(subgroup_file)) {
  stop("Subgroup file not found: ", subgroup_file)
}

# 创建输出目录（如果不存在）
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# 1. 读取数据
cat("Reading data...\n")
group_corr_tissue <- read.table(
  tissue_file,
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

group_corr_subgroup <- read.table(
  subgroup_file,
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

# 检查数据格式
cat("\nData summary:\n")
cat("Tissue data - Columns:", paste(colnames(group_corr_tissue), collapse = ", "), "\n")
cat("Tissue data - Rows:", nrow(group_corr_tissue), "\n")
cat("Tissue data - Unique transcripts:", length(unique(group_corr_tissue$Transcript_ID)), "\n")
cat("Tissues analyzed:", paste(sort(unique(group_corr_tissue$Tissue)), collapse = ", "), "\n\n")

cat("Subgroup data - Columns:", paste(colnames(group_corr_subgroup), collapse = ", "), "\n")
cat("Subgroup data - Rows:", nrow(group_corr_subgroup), "\n")
cat("Subgroup data - Unique transcripts:", length(unique(group_corr_subgroup$Transcript_ID)), "\n")
cat("Subgroups analyzed:", paste(sort(unique(group_corr_subgroup$Subgroup)), collapse = ", "), "\n\n")

# 2. 计算变异性（相关系数范围）
cat("Calculating variability...\n")
calculate_variability <- function(corr_data, group_var) {
  variability <- corr_data %>%
    group_by(Transcript_ID) %>%
    summarise(
      group_count = n_distinct(.data[[group_var]]),
      mean_cor = mean(Group_Cor, na.rm = TRUE),
      median_cor = median(Group_Cor, na.rm = TRUE),
      sd_cor = sd(Group_Cor, na.rm = TRUE),
      min_cor = min(Group_Cor, na.rm = TRUE),
      max_cor = max(Group_Cor, na.rm = TRUE),
      range_cor = max_cor - min_cor,
      iqr_cor = IQR(Group_Cor, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(group_count >= 2) %>%
    arrange(desc(range_cor))
  
  return(variability)
}

tissue_variability <- calculate_variability(group_corr_tissue, "Tissue")
subgroup_variability <- calculate_variability(group_corr_subgroup, "Subgroup")

# 设置阈值
high_var_threshold <- 0.5
extreme_var_threshold <- 0.7

# 统计摘要
cat("\nVariability statistics:\n")
cat("Tissue variability:\n")
cat("  Number of transcripts:", nrow(tissue_variability), "\n")
cat("  Mean range:", round(mean(tissue_variability$range_cor, na.rm = TRUE), 4), "\n")
cat("  Median range:", round(median(tissue_variability$range_cor, na.rm = TRUE), 4), "\n")
cat("  Transcripts with range > 0.5:", sum(tissue_variability$range_cor > 0.5, na.rm = TRUE), 
    sprintf("(%.1f%%)", mean(tissue_variability$range_cor > 0.5, na.rm = TRUE) * 100), "\n")
cat("  Transcripts with range > 0.7:", sum(tissue_variability$range_cor > 0.7, na.rm = TRUE), 
    sprintf("(%.1f%%)", mean(tissue_variability$range_cor > 0.7, na.rm = TRUE) * 100), "\n\n")

cat("Subgroup variability:\n")
cat("  Number of transcripts:", nrow(subgroup_variability), "\n")
cat("  Mean range:", round(mean(subgroup_variability$range_cor, na.rm = TRUE), 4), "\n")
cat("  Median range:", round(median(subgroup_variability$range_cor, na.rm = TRUE), 4), "\n")
cat("  Transcripts with range > 0.5:", sum(subgroup_variability$range_cor > 0.5, na.rm = TRUE), 
    sprintf("(%.1f%%)", mean(subgroup_variability$range_cor > 0.5, na.rm = TRUE) * 100), "\n")
cat("  Transcripts with range > 0.7:", sum(subgroup_variability$range_cor > 0.7, na.rm = TRUE), 
    sprintf("(%.1f%%)", mean(subgroup_variability$range_cor > 0.7, na.rm = TRUE) * 100), "\n\n")

# 3. 准备绘图数据
cat("Preparing data for plotting...\n")
tissue_variability_long <- tissue_variability %>%
  mutate(Group_Type = "Tissue") %>%
  select(Transcript_ID, range_cor, Group_Type)

subgroup_variability_long <- subgroup_variability %>%
  mutate(Group_Type = "Subgroup") %>%
  select(Transcript_ID, range_cor, Group_Type)

combined_range_data <- bind_rows(tissue_variability_long, subgroup_variability_long)

# 4. 创建绘图函数
cat("\nCreating plots...\n")

# 函数：保存图形
save_plot <- function(plot_obj, filename, width = 10, height = 6) {
  png_path <- file.path(output_dir, paste0(filename, ".png"))
  pdf_path <- file.path(output_dir, paste0(filename, ".pdf"))
  
  ggsave(png_path, plot = plot_obj, width = width, height = height, dpi = 300)
  ggsave(pdf_path, plot = plot_obj, width = width, height = height)
  
  cat("  Saved:", basename(png_path), "and", basename(pdf_path), "\n")
}

# 绘图1：重叠密度图
p1 <- ggplot(combined_range_data, aes(x = range_cor, fill = Group_Type, color = Group_Type)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  geom_vline(xintercept = high_var_threshold, linetype = "dashed", color = "red", alpha = 0.7, linewidth = 0.8) +
  geom_vline(xintercept = extreme_var_threshold, linetype = "dashed", color = "darkred", alpha = 0.7, linewidth = 0.8) +
  scale_fill_manual(values = c("Tissue" = "steelblue", "Subgroup" = "darkorange")) +
  scale_color_manual(values = c("Tissue" = "steelblue", "Subgroup" = "darkorange")) +
  labs(
    title = "Comparison of Variability Distributions",
    subtitle = sprintf("Tissue (n=%d) vs Subgroup (n=%d)", 
                      nrow(tissue_variability), nrow(subgroup_variability)),
    x = "Range of Correlation Coefficients",
    y = "Density",
    fill = "Group Type",
    color = "Group Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  annotate("text", x = Inf, y = Inf, 
           label = sprintf("Tissue: μ=%.3f, M=%.3f\nSubgroup: μ=%.3f, M=%.3f",
                          mean(tissue_variability$range_cor, na.rm = TRUE),
                          median(tissue_variability$range_cor, na.rm = TRUE),
                          mean(subgroup_variability$range_cor, na.rm = TRUE),
                          median(subgroup_variability$range_cor, na.rm = TRUE)),
           hjust = 1.1, vjust = 1.1, size = 4, color = "black")

save_plot(p1, "Variability_overlay_density")

# 绘图2：直方图+密度曲线
p2 <- ggplot(combined_range_data, aes(x = range_cor, fill = Group_Type)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 40, 
                 aes(y = ..density..), color = "white", size = 0.1) +
  geom_density(alpha = 0.3, aes(color = Group_Type), linewidth = 0.8) +
  geom_vline(xintercept = high_var_threshold, linetype = "dashed", color = "red", alpha = 0.7, linewidth = 0.8) +
  geom_vline(xintercept = extreme_var_threshold, linetype = "dashed", color = "darkred", alpha = 0.7, linewidth = 0.8) +
  scale_fill_manual(values = c("Tissue" = "steelblue", "Subgroup" = "darkorange")) +
  scale_color_manual(values = c("Tissue" = "steelblue", "Subgroup" = "darkorange")) +
  labs(
    title = "Overlay of Variability Distributions",
    subtitle = "Histogram with density curves",
    x = "Range of Correlation Coefficients",
    y = "Density",
    fill = "Group Type",
    color = "Group Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

save_plot(p2, "Variability_overlay_histogram")

# 绘图3：箱线图+小提琴图
# 执行统计检验
wilcox_test <- wilcox.test(range_cor ~ Group_Type, data = combined_range_data)

p3 <- ggplot(combined_range_data, aes(x = Group_Type, y = range_cor, fill = Group_Type)) +
  geom_violin(alpha = 0.5, trim = FALSE, width = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 0.8, color = "gray20") +
  geom_hline(yintercept = high_var_threshold, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_hline(yintercept = extreme_var_threshold, linetype = "dashed", color = "darkred", alpha = 0.7) +
  scale_fill_manual(values = c("Tissue" = "steelblue", "Subgroup" = "darkorange")) +
  labs(
    title = "Distribution Comparison: Tissue vs Subgroup Variability",
    subtitle = sprintf("Tissue median: %.3f, Subgroup median: %.3f",
                      median(tissue_variability$range_cor, na.rm = TRUE),
                      median(subgroup_variability$range_cor, na.rm = TRUE)),
    x = "Group Type",
    y = "Range of Correlation Coefficients",
    fill = "Group Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "none"
  ) +
  annotate("text", x = 1.5, y = max(combined_range_data$range_cor, na.rm = TRUE) * 0.95,
           label = paste("Wilcoxon rank-sum test:\n",
                        sprintf("W = %.0f", wilcox_test$statistic),
                        sprintf("\np = %.2e", wilcox_test$p.value)),
           size = 4, hjust = 0.5, vjust = 1)

save_plot(p3, "Variability_box_violin", width = 8, height = 6)

# 绘图4：累积分布函数图
p4 <- ggplot(combined_range_data, aes(x = range_cor, color = Group_Type)) +
  stat_ecdf(linewidth = 1.2, alpha = 0.8) +
  geom_vline(xintercept = high_var_threshold, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_vline(xintercept = extreme_var_threshold, linetype = "dashed", color = "darkred", alpha = 0.7) +
  scale_color_manual(values = c("Tissue" = "steelblue", "Subgroup" = "darkorange")) +
  labs(
    title = "Cumulative Distribution of Variability",
    subtitle = "Proportion of transcripts below given range value",
    x = "Range of Correlation Coefficients",
    y = "Cumulative Proportion",
    color = "Group Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "bottom"
  ) +
  annotate("text", x = high_var_threshold, y = 0.95, 
           label = sprintf("High variability threshold\n(%.1f)", high_var_threshold),
           hjust = 1.1, color = "red", size = 3.5) +
  annotate("text", x = extreme_var_threshold, y = 0.95, 
           label = sprintf("Extreme variability threshold\n(%.1f)", extreme_var_threshold),
           hjust = -0.1, color = "darkred", size = 3.5)

save_plot(p4, "Variability_ecdf")

# 绘图5：并排直方图（原始版本）
p5a <- ggplot(tissue_variability, aes(x = range_cor)) +
  geom_histogram(bins = 40, fill = "steelblue", alpha = 0.8) +
  geom_vline(xintercept = high_var_threshold, linetype = "dashed", color = "red") +
  geom_vline(xintercept = extreme_var_threshold, linetype = "dashed", color = "darkred") +
  labs(
    title = "Tissue Variability",
    x = "Correlation Range",
    y = "Transcript Count"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p5b <- ggplot(subgroup_variability, aes(x = range_cor)) +
  geom_histogram(bins = 40, fill = "darkorange", alpha = 0.8) +
  geom_vline(xintercept = high_var_threshold, linetype = "dashed", color = "red") +
  geom_vline(xintercept = extreme_var_threshold, linetype = "dashed", color = "darkred") +
  labs(
    title = "Subgroup Variability",
    x = "Correlation Range",
    y = "Transcript Count"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p5 <- p5a + p5b + plot_annotation(title = "Variability Distribution: Tissue vs Subgroup")

save_plot(p5, "Variability_side_by_side", width = 12, height = 5)

# 5. 保存变异性数据
cat("\nSaving variability data...\n")
write.table(
  tissue_variability,
  file = file.path(output_dir, "Tissue_variability_results.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  subgroup_variability,
  file = file.path(output_dir, "Subgroup_variability_results.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  combined_range_data,
  file = file.path(output_dir, "Combined_variability_data.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# 6. 生成总结报告
summary_df <- data.frame(
  Metric = c(
    "Tissue transcripts analyzed",
    "Subgroup transcripts analyzed",
    "Tissue mean range",
    "Tissue median range",
    "Tissue high variability (>0.5)",
    "Tissue extreme variability (>0.7)",
    "Subgroup mean range",
    "Subgroup median range",
    "Subgroup high variability (>0.5)",
    "Subgroup extreme variability (>0.7)",
    "Wilcoxon test statistic (W)",
    "Wilcoxon test p-value"
  ),
  Value = c(
    nrow(tissue_variability),
    nrow(subgroup_variability),
    round(mean(tissue_variability$range_cor, na.rm = TRUE), 4),
    round(median(tissue_variability$range_cor, na.rm = TRUE), 4),
    paste(sum(tissue_variability$range_cor > 0.5, na.rm = TRUE),
          sprintf("(%.1f%%)", mean(tissue_variability$range_cor > 0.5, na.rm = TRUE) * 100)),
    paste(sum(tissue_variability$range_cor > 0.7, na.rm = TRUE),
          sprintf("(%.1f%%)", mean(tissue_variability$range_cor > 0.7, na.rm = TRUE) * 100)),
    round(mean(subgroup_variability$range_cor, na.rm = TRUE), 4),
    round(median(subgroup_variability$range_cor, na.rm = TRUE), 4),
    paste(sum(subgroup_variability$range_cor > 0.5, na.rm = TRUE),
          sprintf("(%.1f%%)", mean(subgroup_variability$range_cor > 0.5, na.rm = TRUE) * 100)),
    paste(sum(subgroup_variability$range_cor > 0.7, na.rm = TRUE),
          sprintf("(%.1f%%)", mean(subgroup_variability$range_cor > 0.7, na.rm = TRUE) * 100)),
    round(wilcox_test$statistic, 2),
    format(wilcox_test$p.value, scientific = TRUE, digits = 3)
  )
)

write.table(
  summary_df,
  file = file.path(output_dir, "Variability_summary.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# 7. 完成
cat("\n========================================\n")
cat("Plotting complete!\n")
cat("Output directory:", output_dir, "\n")
cat("\nGenerated files:\n")
cat("1. Plots (PNG & PDF):\n")
cat("   - Variability_overlay_density.*\n")
cat("   - Variability_overlay_histogram.*\n")
cat("   - Variability_box_violin.*\n")
cat("   - Variability_ecdf.*\n")
cat("   - Variability_side_by_side.*\n")
cat("\n2. Data files:\n")
cat("   - Tissue_variability_results.txt\n")
cat("   - Subgroup_variability_results.txt\n")
cat("   - Combined_variability_data.txt\n")
cat("   - Variability_summary.txt\n")
cat("\nKey finding: Tissue shows ", 
    round(mean(tissue_variability$range_cor, na.rm = TRUE) / 
          mean(subgroup_variability$range_cor, na.rm = TRUE), 2),
    "x higher variability than Subgroup\n", sep = "")
cat("========================================\n")
