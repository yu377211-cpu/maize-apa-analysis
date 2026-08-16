library(tidyverse)
library(ggplot2)
library(patchwork)  # 用于组合图形

# 设置主题 - 白色背景，无网格线
white_theme <- theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA)
  )

# 读取已有结果
global_corr <- read.table(
  "PDUI_PAC_global_correlation.txt",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

group_corr_tissue <- read.table(
  "PDUI_PAC_tissue_correlation.txt",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

group_corr_subgroup <- read.table(
  "PDUI_PAC_subgroup_correlation.txt",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

cat("=== Comprehensive APA Regulation Analysis (Corrected) ===\n")
cat("Analyzing both Tissue and Subgroup specificity\n\n")

# 1. 数据质量检查
cat("1. Data Overview:\n")
cat("Global correlation - Total transcripts:", nrow(global_corr), "\n")
cat("Tissue correlation - Transcripts:", length(unique(group_corr_tissue$Transcript_ID)), "\n")
cat("Subgroup correlation - Transcripts:", length(unique(group_corr_subgroup$Transcript_ID)), "\n")

cat("\nTissues analyzed (", length(unique(group_corr_tissue$Tissue)), "): ",
    paste(sort(unique(group_corr_tissue$Tissue)), collapse = ", "), "\n", sep = "")

cat("Subgroups analyzed (", length(unique(group_corr_subgroup$Subgroup)), "): ",
    paste(sort(unique(group_corr_subgroup$Subgroup)), collapse = ", "), "\n\n", sep = "")

# 2. 全局相关性总结
cat("2. Global Correlation Summary:\n")
significant_global <- global_corr %>% filter(FDR < 0.05)
cat("Significant transcripts (FDR<0.05):", nrow(significant_global),
    sprintf("(%.1f%%)", nrow(significant_global)/nrow(global_corr)*100), "\n")

cor_summary <- data.frame(
  Type = c("Positive", "Negative", "Strong Positive (r>0.5)", "Strong Negative (r<-0.5)"),
  Count = c(
    sum(global_corr$Correlation > 0, na.rm = TRUE),
    sum(global_corr$Correlation < 0, na.rm = TRUE),
    sum(global_corr$Correlation > 0.5, na.rm = TRUE),
    sum(global_corr$Correlation < -0.5, na.rm = TRUE)
  ),
  Percentage = c(
    sum(global_corr$Correlation > 0, na.rm = TRUE)/nrow(global_corr)*100,
    sum(global_corr$Correlation < 0, na.rm = TRUE)/nrow(global_corr)*100,
    sum(global_corr$Correlation > 0.5, na.rm = TRUE)/nrow(global_corr)*100,
    sum(global_corr$Correlation < -0.5, na.rm = TRUE)/nrow(global_corr)*100
  )
)

print(cor_summary)
cat("\n")

# 3. 并行分析Tissue和Subgroup的变异性
analyze_variability <- function(corr_data, group_var, group_name) {
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

cat("3. Variability Analysis:\n")

# 分析Tissue和Subgroup
tissue_variability <- analyze_variability(group_corr_tissue, "Tissue", "Tissue")
subgroup_variability <- analyze_variability(group_corr_subgroup, "Subgroup", "Subgroup")

# 定义阈值
high_var_threshold <- 0.5
extreme_var_threshold <- 0.7

# Tissue结果
tissue_high_var <- tissue_variability %>% filter(range_cor > high_var_threshold)
tissue_extreme_var <- tissue_variability %>% filter(range_cor > extreme_var_threshold)

cat("\nTissue Variability:\n")
cat("  Transcripts analyzed:", nrow(tissue_variability), "\n")
cat("  Average correlation range:", round(mean(tissue_variability$range_cor, na.rm = TRUE), 4), "\n")
cat("  High variability (range>0.5):", nrow(tissue_high_var),
    sprintf("(%.1f%%)", nrow(tissue_high_var)/nrow(tissue_variability)*100), "\n")
cat("  Extreme variability (range>0.7):", nrow(tissue_extreme_var),
    sprintf("(%.1f%%)", nrow(tissue_extreme_var)/nrow(tissue_variability)*100), "\n")

# Subgroup结果
subgroup_high_var <- subgroup_variability %>% filter(range_cor > high_var_threshold)
subgroup_extreme_var <- subgroup_variability %>% filter(range_cor > extreme_var_threshold)

cat("\nSubgroup Variability:\n")
cat("  Transcripts analyzed:", nrow(subgroup_variability), "\n")
cat("  Average correlation range:", round(mean(subgroup_variability$range_cor, na.rm = TRUE), 4), "\n")
cat("  High variability (range>0.5):", nrow(subgroup_high_var),
    sprintf("(%.1f%%)", nrow(subgroup_high_var)/nrow(subgroup_variability)*100), "\n")
cat("  Extreme variability (range>0.7):", nrow(subgroup_extreme_var),
    sprintf("(%.1f%%)", nrow(subgroup_extreme_var)/nrow(subgroup_variability)*100), "\n")

# 4. 比较Tissue和Subgroup特异性 - 修正版
cat("\n4. Comparing Tissue vs Subgroup Specificity (Corrected):\n")

# 合并两个变异性的结果进行比较
combined_variability <- tissue_variability %>%
  select(Transcript_ID, tissue_range = range_cor, tissue_groups = group_count,
         tissue_mean = mean_cor, tissue_min = min_cor, tissue_max = max_cor) %>%
  inner_join(
    subgroup_variability %>%
      select(Transcript_ID, subgroup_range = range_cor, subgroup_groups = group_count,
             subgroup_mean = mean_cor, subgroup_min = min_cor, subgroup_max = max_cor),
    by = "Transcript_ID"
  ) %>%
  mutate(
    range_difference = tissue_range - subgroup_range,
    max_range = pmax(tissue_range, subgroup_range),
    # 修正的分类逻辑：两者都>0.5就是双重特异性
    is_tissue_specific = tissue_range > 0.5,
    is_subgroup_specific = subgroup_range > 0.5,
    specificity_type = case_when(
      is_tissue_specific & is_subgroup_specific ~ "Dual-specific",
      is_tissue_specific & !is_subgroup_specific ~ "Tissue-specific",
      !is_tissue_specific & is_subgroup_specific ~ "Subgroup-specific",
      TRUE ~ "Non-specific"
    )
  )

# 统计特异性类型
specificity_counts <- combined_variability %>%
  group_by(specificity_type) %>%
  summarise(
    count = n(),
    percentage = n()/nrow(combined_variability)*100,
    mean_tissue_range = mean(tissue_range, na.rm = TRUE),
    mean_subgroup_range = mean(subgroup_range, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(count))

cat("\nSpecificity Classification (Corrected):\n")
print(specificity_counts)

# 识别双重特异性转录本
dual_specific <- combined_variability %>%
  filter(specificity_type == "Dual-specific") %>%
  arrange(desc(tissue_range + subgroup_range))

if(nrow(dual_specific) > 0) {
  cat("\nDual-specific transcripts found:", nrow(dual_specific), "\n")
  cat("Top dual-specific transcripts (both Tissue and Subgroup):\n")
  print(head(dual_specific, 10))
} else {
  cat("\nNo dual-specific transcripts found.\n")
}

# 5. 模式分类分析
cat("\n5. Pattern Classification:\n")

categorize_patterns <- function(variability_data, corr_data, group_var, group_name) {
  categorized <- variability_data %>%
    left_join(
      corr_data %>%
        group_by(Transcript_ID) %>%
        summarise(
          group_count = n_distinct(.data[[group_var]]),
          .groups = "drop"
        ),
      by = "Transcript_ID"
    ) %>%
    mutate(
      pattern_type = case_when(
        min_cor > 0.3 & max_cor > 0.5 ~ "Consistently Positive",
        max_cor < -0.3 & min_cor < -0.5 ~ "Consistently Negative",
        min_cor < -0.3 & max_cor > 0.3 & range_cor > 0.6 ~ "Mixed/Context-dependent",
        min_cor < -0.2 & max_cor > 0.2 & range_cor > 0.4 ~ "Variable",
        abs(mean_cor) < 0.2 & range_cor < 0.3 ~ "Weak/No Correlation",
        TRUE ~ "Moderate/Other"
      ),
      group_type = group_name
    )

  return(categorized)
}

# 分类两种分组模式
tissue_categorized <- categorize_patterns(tissue_variability, group_corr_tissue, "Tissue", "Tissue")
subgroup_categorized <- categorize_patterns(subgroup_variability, group_corr_subgroup, "Subgroup", "Subgroup")

# 合并模式结果
combined_patterns <- bind_rows(tissue_categorized, subgroup_categorized)

cat("\nPattern Distribution by Group Type:\n")
pattern_summary <- combined_patterns %>%
  group_by(group_type, pattern_type) %>%
  summarise(
    count = n(),
    .groups = "drop"
  ) %>%
  group_by(group_type) %>%
  mutate(
    percentage = count/sum(count)*100,
    pattern_type = factor(pattern_type, levels = c(
      "Consistently Positive", "Consistently Negative",
      "Mixed/Context-dependent", "Variable",
      "Moderate/Other", "Weak/No Correlation"
    ))
  ) %>%
  arrange(group_type, desc(count))

print(pattern_summary)

# 6. 可视化 (修正版，使用白色主题)
cat("\n6. Generating Visualizations (with white background)...\n")

# 6.1 重叠直方图（替换原来的并排图）
# 准备合并数据
tissue_variability_long <- tissue_variability %>%
  mutate(Group_Type = "Tissue") %>%
  select(Transcript_ID, range_cor, Group_Type)

subgroup_variability_long <- subgroup_variability %>%
  mutate(Group_Type = "Subgroup") %>%
  select(Transcript_ID, range_cor, Group_Type)

combined_range_data <- bind_rows(tissue_variability_long, subgroup_variability_long)

p1 <- ggplot(combined_range_data, aes(x = range_cor, fill = Group_Type)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 40,
                 aes(y = ..density..), color = "white", size = 0.1) +
  geom_density(alpha = 0.3, aes(color = Group_Type), linewidth = 0.8) +
  geom_vline(xintercept = high_var_threshold, linetype = "dashed",
             color = "red", alpha = 0.7, linewidth = 0.8) +
  geom_vline(xintercept = extreme_var_threshold, linetype = "dashed",
             color = "darkred", alpha = 0.7, linewidth = 0.8) +
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
  white_theme +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "bottom",
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

ggsave("Variability_overlay_histogram.png", plot = p1, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("Variability_overlay_histogram.pdf", plot = p1, width = 10, height = 6, bg = "white")

# 6.2 特异性类型分布（修正版）
p2 <- ggplot(specificity_counts, aes(x = reorder(specificity_type, -count), y = count)) +
  geom_bar(stat = "identity", aes(fill = specificity_type), alpha = 0.8, width = 0.7) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", count, percentage)),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c(
    "Dual-specific" = "#E41A1C",
    "Tissue-specific" = "#377EB8",
    "Subgroup-specific" = "#FF7F00",
    "Non-specific" = "#999999"
  )) +
  labs(
    title = "Specificity Classification of Transcripts",
    subtitle = "Based on correlation variability in Tissue and Subgroup",
    x = "Specificity Type",
    y = "Number of Transcripts"
  ) +
  white_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "none",
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave("Specificity_classification.png", plot = p2, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("Specificity_classification.pdf", plot = p2, width = 10, height = 6, bg = "white")

# 6.3 Tissue和Subgroup范围比较散点图（修正版）
p3 <- ggplot(combined_variability, aes(x = tissue_range, y = subgroup_range)) +
  geom_point(aes(color = specificity_type), alpha = 0.7, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", alpha = 0.7, linewidth = 0.8) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red", alpha = 0.7, linewidth = 0.8) +
  scale_color_manual(values = c(
    "Dual-specific" = "#E41A1C",
    "Tissue-specific" = "#377EB8",
    "Subgroup-specific" = "#FF7F00",
    "Non-specific" = "#999999"
  )) +
  labs(
    title = "Comparison of Tissue vs Subgroup Specificity",
    subtitle = sprintf("Each point represents one transcript | Dual-specific: %d",
                      sum(combined_variability$specificity_type == "Dual-specific")),
    x = "Tissue Correlation Range",
    y = "Subgroup Correlation Range",
    color = "Specificity Type"
  ) +
  white_theme +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "right",
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  # 添加象限标注
  annotate("text", x = 1.5, y = 1.5,
           label = paste("Dual-specific\nn =", sum(combined_variability$specificity_type == "Dual-specific")),
           color = "#E41A1C", size = 4.5, fontface = "bold") +
  annotate("text", x = 1.5, y = 0.25,
           label = paste("Tissue-specific\nn =", sum(combined_variability$specificity_type == "Tissue-specific")),
           color = "#377EB8", size = 4.5, fontface = "bold") +
  annotate("text", x = 0.25, y = 1.5,
           label = paste("Subgroup-specific\nn =", sum(combined_variability$specificity_type == "Subgroup-specific")),
           color = "#FF7F00", size = 4.5, fontface = "bold")

ggsave("Tissue_vs_subgroup_scatter.png", plot = p3, width = 12, height = 8, dpi = 300, bg = "white")
ggsave("Tissue_vs_subgroup_scatter.pdf", plot = p3, width = 12, height = 8, bg = "white")

# 6.4 模式分类比较（并排条形图，白色背景）
p4 <- ggplot(pattern_summary, aes(x = pattern_type, y = percentage, fill = group_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.8, width = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", percentage)),
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("Tissue" = "#377EB8", "Subgroup" = "#FF7F00")) +
  labs(
    title = "Correlation Pattern Distribution",
    subtitle = "Comparison between Tissue and Subgroup",
    x = "Pattern Type",
    y = "Percentage of Transcripts",
    fill = "Group Type"
  ) +
  white_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "bottom",
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave("Pattern_comparison.png", plot = p4, width = 12, height = 6, dpi = 300, bg = "white")
ggsave("Pattern_comparison.pdf", plot = p4, width = 12, height = 6, bg = "white")

# 6.5 全局相关性分布图（新增）
p5 <- ggplot(global_corr, aes(x = Correlation)) +
  geom_histogram(bins = 50, fill = "#4DAF4A", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted", color = "gray50", alpha = 0.5) +
  labs(
    title = "Distribution of Global PDUI-PAC Correlations",
    subtitle = sprintf("Total transcripts: %d, Significant (FDR<0.05): %d (%.1f%%)",
                      nrow(global_corr), nrow(significant_global),
                      nrow(significant_global)/nrow(global_corr)*100),
    x = "Spearman Correlation Coefficient",
    y = "Number of Transcripts"
  ) +
  white_theme +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

ggsave("Global_correlation_distribution.png", plot = p5, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("Global_correlation_distribution.pdf", plot = p5, width = 10, height = 6, bg = "white")

# 6.6 双重特异性转录本的热点图
if(nrow(dual_specific) >= 10) {
  set.seed(123)
  sample_dual <- sample(dual_specific$Transcript_ID, min(15, nrow(dual_specific)))

  # 获取Tissue数据
  tissue_sample <- group_corr_tissue %>%
    filter(Transcript_ID %in% sample_dual) %>%
    select(Transcript_ID, Tissue, Group_Cor) %>%
    pivot_wider(names_from = Tissue, values_from = Group_Cor) %>%
    column_to_rownames("Transcript_ID")

  # 获取Subgroup数据
  subgroup_sample <- group_corr_subgroup %>%
    filter(Transcript_ID %in% sample_dual) %>%
    select(Transcript_ID, Subgroup, Group_Cor) %>%
    pivot_wider(names_from = Subgroup, values_from = Group_Cor) %>%
    column_to_rownames("Transcript_ID")

  # 绘制热图（需要pheatmap包）
  if(requireNamespace("pheatmap", quietly = TRUE)) {
    library(pheatmap)

    # 合并数据
    combined_sample <- cbind(
      tissue_sample[rownames(tissue_sample), ],
      subgroup_sample[rownames(tissue_sample), ]
    )

    # 绘制热图（PNG）
    pheatmap(
      as.matrix(combined_sample),
      scale = "row",
      color = colorRampPalette(c("blue", "white", "red"))(100),
      main = "Dual-Specific Transcripts: Tissue and Subgroup Correlation Patterns",
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = 8,
      fontsize_col = 8,
      filename = "Dual_specific_heatmap.png",
      width = 12,
      height = 10
    )

    # 绘制热图（PDF）
    pheatmap(
      as.matrix(combined_sample),
      scale = "row",
      color = colorRampPalette(c("blue", "white", "red"))(100),
      main = "Dual-Specific Transcripts: Tissue and Subgroup Correlation Patterns",
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = 8,
      fontsize_col = 8,
      filename = "Dual_specific_heatmap.pdf",
      width = 12,
      height = 10
    )

    cat("Generated dual-specific heatmap\n")
  }
}

# 7. 保存结果文件
cat("\n7. Saving Results...\n")

# 7.1 变异性分析结果
write.table(
  tissue_variability,
  file = "Tissue_variability_analysis.txt",
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  subgroup_variability,
  file = "Subgroup_variability_analysis.txt",
  sep = "\t", row.names = FALSE, quote = FALSE
)

# 7.2 特异性比较结果（修正版）
write.table(
  combined_variability,
  file = "Specificity_comparison_results_corrected.txt",
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  specificity_counts,
  file = "Specificity_classification_summary_corrected.txt",
  sep = "\t", row.names = FALSE, quote = FALSE
)

# 7.3 模式分类结果
write.table(
  combined_patterns,
  file = "Pattern_classification_results.txt",
  sep = "\t", row.names = FALSE, quote = FALSE
)

# 7.4 双重特异性转录本
if(nrow(dual_specific) > 0) {
  write.table(
    dual_specific,
    file = "Dual_specific_transcripts_corrected.txt",
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}

# 7.5 详细的双重特异性转录本数据
if(nrow(dual_specific) > 0) {
  top_dual <- head(dual_specific$Transcript_ID, 20)

  dual_details <- bind_rows(
    group_corr_tissue %>%
      filter(Transcript_ID %in% top_dual) %>%
      mutate(Group_Type = "Tissue", Group = Tissue) %>%
      select(Transcript_ID, Group_Type, Group, Group_Cor),
    group_corr_subgroup %>%
      filter(Transcript_ID %in% top_dual) %>%
      mutate(Group_Type = "Subgroup", Group = Subgroup) %>%
      select(Transcript_ID, Group_Type, Group, Group_Cor)
  )

  write.table(
    dual_details,
    file = "Top_dual_specific_details.txt",
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}

# 8. 生成总结报告（修正版）
cat("\n8. Generating Corrected Summary Report...\n")

summary_report_corrected <- data.frame(
  Analysis = c(
    "Global Analysis",
    "Global Analysis",
    "Global Analysis",
    "Tissue Variability",
    "Tissue Variability",
    "Tissue Variability",
    "Subgroup Variability",
    "Subgroup Variability",
    "Subgroup Variability",
    "Specificity Comparison (Corrected)",
    "Specificity Comparison (Corrected)",
    "Specificity Comparison (Corrected)",
    "Specificity Comparison (Corrected)",
    "Pattern Analysis",
    "Pattern Analysis"
  ),
  Metric = c(
    "Total transcripts",
    "Significant transcripts (FDR<0.05)",
    "Strong positive correlations (r>0.5)",
    "Transcripts with high variability (range>0.5)",
    "Transcripts with extreme variability (range>0.7)",
    "Average correlation range",
    "Transcripts with high variability (range>0.5)",
    "Transcripts with extreme variability (range>0.7)",
    "Average correlation range",
    "Tissue-specific transcripts",
    "Subgroup-specific transcripts",
    "Dual-specific transcripts (Corrected)",
    "Non-specific transcripts",
    "Most common tissue pattern",
    "Most common subgroup pattern"
  ),
  Value = c(
    nrow(global_corr),
    sprintf("%d (%.1f%%)", nrow(significant_global), nrow(significant_global)/nrow(global_corr)*100),
    sum(global_corr$Correlation > 0.5, na.rm = TRUE),
    sprintf("%d (%.1f%%)", nrow(tissue_high_var), nrow(tissue_high_var)/nrow(tissue_variability)*100),
    sprintf("%d (%.1f%%)", nrow(tissue_extreme_var), nrow(tissue_extreme_var)/nrow(tissue_variability)*100),
    round(mean(tissue_variability$range_cor, na.rm = TRUE), 3),
    sprintf("%d (%.1f%%)", nrow(subgroup_high_var), nrow(subgroup_high_var)/nrow(subgroup_variability)*100),
    sprintf("%d (%.1f%%)", nrow(subgroup_extreme_var), nrow(subgroup_extreme_var)/nrow(subgroup_variability)*100),
    round(mean(subgroup_variability$range_cor, na.rm = TRUE), 3),
    sum(combined_variability$specificity_type == "Tissue-specific", na.rm = TRUE),
    sum(combined_variability$specificity_type == "Subgroup-specific", na.rm = TRUE),
    sum(combined_variability$specificity_type == "Dual-specific", na.rm = TRUE),
    sum(combined_variability$specificity_type == "Non-specific", na.rm = TRUE),
    pattern_summary %>% filter(group_type == "Tissue") %>% slice_max(count) %>% pull(pattern_type),
    pattern_summary %>% filter(group_type == "Subgroup") %>% slice_max(count) %>% pull(pattern_type)
  )
)

write.table(
  summary_report_corrected,
  file = "Comprehensive_analysis_summary_corrected.txt",
  sep = "\t", row.names = FALSE, quote = FALSE
)

cat("\nCorrected Summary Report:\n")
print(summary_report_corrected)

# 9. 最终总结
cat("\n=== ANALYSIS COMPLETE (CORRECTED) ===\n")
cat("\nKey Findings (Corrected):\n")
cat("1. Global correlations:", nrow(significant_global),
    sprintf("significant transcripts (%.1f%%)\n", nrow(significant_global)/nrow(global_corr)*100))
cat("2. Tissue-specific:", sum(combined_variability$specificity_type == "Tissue-specific"),
    "transcripts show tissue-specific patterns\n")
cat("3. Subgroup-specific:", sum(combined_variability$specificity_type == "Subgroup-specific"),
    "transcripts show subgroup-specific patterns\n")
cat("4. Dual-specific (CORRECTED):", sum(combined_variability$specificity_type == "Dual-specific"),
    "transcripts show both tissue and subgroup specificity\n")
cat("5. Non-specific:", sum(combined_variability$specificity_type == "Non-specific"),
    "transcripts\n")

cat("\nFiles Generated:\n")
cat("- Visualizations (with white background):\n")
cat("  1. Variability_overlay_histogram.*\n")
cat("  2. Specificity_classification.*\n")
cat("  3. Tissue_vs_subgroup_scatter.*\n")
cat("  4. Pattern_comparison.*\n")
cat("  5. Global_correlation_distribution.*\n")
if(nrow(dual_specific) >= 10) {
  cat("  6. Dual_specific_heatmap.*\n")
}

cat("\n- Data files (Corrected):\n")
cat("  1. Tissue_variability_analysis.txt\n")
cat("  2. Subgroup_variability_analysis.txt\n")
cat("  3. Specificity_comparison_results_corrected.txt\n")
cat("  4. Specificity_classification_summary_corrected.txt\n")
cat("  5. Pattern_classification_results.txt\n")
if(nrow(dual_specific) > 0) {
  cat("  6. Dual_specific_transcripts_corrected.txt\n")
  cat("  7. Top_dual_specific_details.txt\n")
}
cat("  8. Comprehensive_analysis_summary_corrected.txt\n")

cat("\nBiological Interpretation:\n")
cat("1. Tissue-specific transcripts reflect tissue-specific APA regulation\n")
cat("2. Subgroup-specific transcripts reflect developmental/physiological state differences\n")
cat("3. Dual-specific transcripts represent genes under complex, multi-level regulation\n")
cat("4. Negative correlations dominate, suggesting distal APA usage often reduces expression\n")

cat("\nCorrection Note: Previous report of 0 dual-specific transcripts was incorrect.\n")
cat("Correct analysis shows", sum(combined_variability$specificity_type == "Dual-specific"),
    "dual-specific transcripts.\n")

cat("\nDone!\n")
