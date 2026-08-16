library(tidyverse)

# 读取数据
global_corr <- read.table("PDUI_PAC_global_correlation.txt", header = TRUE, sep = "\t")
group_corr_tissue <- read.table("PDUI_PAC_tissue_correlation.txt", header = TRUE, sep = "\t")
group_corr_subgroup <- read.table("PDUI_PAC_subgroup_correlation.txt", header = TRUE, sep = "\t")

cat("=== Comprehensive Global Correlation Analysis ===\n\n")

# 1. 全局相关性统计
cat("1. Global Correlation Statistics:\n")
cat("Total transcripts:", nrow(global_corr), "\n")
cat("Significant transcripts (FDR<0.05):", sum(global_corr$FDR < 0.05, na.rm = TRUE),
    sprintf("(%.1f%%)", mean(global_corr$FDR < 0.05, na.rm = TRUE)*100), "\n\n")

# 分类全局相关性
global_corr_classified <- global_corr %>%
  mutate(
    correlation_type = case_when(
      Correlation > 0.5 ~ "Strong Positive",
      Correlation > 0.3 ~ "Moderate Positive",
      Correlation > 0.1 ~ "Weak Positive",
      Correlation > -0.1 ~ "Neutral",
      Correlation > -0.3 ~ "Weak Negative",
      Correlation > -0.5 ~ "Moderate Negative",
      TRUE ~ "Strong Negative"
    ),
    correlation_type = factor(correlation_type, levels = c(
      "Strong Positive", "Moderate Positive", "Weak Positive",
      "Neutral", "Weak Negative", "Moderate Negative", "Strong Negative"
    )),
    sign_type = ifelse(Correlation > 0, "Positive", "Negative"),
    strength_type = case_when(
      abs(Correlation) > 0.5 ~ "Strong",
      abs(Correlation) > 0.3 ~ "Moderate",
      abs(Correlation) > 0.1 ~ "Weak",
      TRUE ~ "Neutral"
    )
  )

# 统计
cat("Correlation Type Distribution:\n")
type_counts <- global_corr_classified %>%
  group_by(correlation_type) %>%
  summarise(
    count = n(),
    percentage = n()/nrow(global_corr)*100,
    mean_cor = mean(Correlation, na.rm = TRUE),
    .groups = "drop"
  )

print(type_counts)
cat("\n")

cat("Summary by Sign and Strength:\n")
sign_summary <- global_corr_classified %>%
  group_by(sign_type, strength_type) %>%
  summarise(
    count = n(),
    percentage = n()/nrow(global_corr)*100,
    .groups = "drop"
  ) %>%
  arrange(desc(sign_type), desc(strength_type))

print(sign_summary)
cat("\n")

# 2. 提取强全局正相关和负相关转录本
cat("2. Extracting Strong Global Correlations:\n")

# 强正相关
strong_positive <- global_corr_classified %>%
  filter(correlation_type == "Strong Positive") %>%
  arrange(desc(Correlation))

cat("Strong Positive (r > 0.5):", nrow(strong_positive),
    sprintf("(%.1f%%)", nrow(strong_positive)/nrow(global_corr)*100), "\n")

# 强负相关
strong_negative <- global_corr_classified %>%
  filter(correlation_type == "Strong Negative") %>%
  arrange(Correlation)  # 从最负到最不負

cat("Strong Negative (r < -0.5):", nrow(strong_negative),
    sprintf("(%.1f%%)", nrow(strong_negative)/nrow(global_corr)*100), "\n")

# 3. 计算变异性
cat("\n3. Calculating Variability:\n")

tissue_range <- group_corr_tissue %>%
  group_by(Transcript_ID) %>%
  summarise(
    tissue_range = max(Group_Cor, na.rm = TRUE) - min(Group_Cor, na.rm = TRUE),
    tissue_mean = mean(Group_Cor, na.rm = TRUE),
    tissue_median = median(Group_Cor, na.rm = TRUE),
    tissue_max = max(Group_Cor, na.rm = TRUE),
    tissue_min = min(Group_Cor, na.rm = TRUE),
    .groups = "drop"
  )

subgroup_range <- group_corr_subgroup %>%
  group_by(Transcript_ID) %>%
  summarise(
    subgroup_range = max(Group_Cor, na.rm = TRUE) - min(Group_Cor, na.rm = TRUE),
    subgroup_mean = mean(Group_Cor, na.rm = TRUE),
    subgroup_median = median(Group_Cor, na.rm = TRUE),
    subgroup_max = max(Group_Cor, na.rm = TRUE),
    subgroup_min = min(Group_Cor, na.rm = TRUE),
    .groups = "drop"
  )

# 4. 识别双重特异性（包括正负相关考虑）
cat("\n4. Identifying Dual-Specific Transcripts:\n")

dual_specific <- tissue_range %>%
  inner_join(subgroup_range, by = "Transcript_ID") %>%
  filter(tissue_range > 0.5 & subgroup_range > 0.5) %>%
  mutate(
    combined_range = tissue_range + subgroup_range,
    max_range = pmax(tissue_range, subgroup_range),
    # 根据平均相关性分类
    tissue_direction = ifelse(tissue_mean > 0, "Positive", "Negative"),
    subgroup_direction = ifelse(subgroup_mean > 0, "Positive", "Negative"),
    direction_consistency = ifelse(tissue_direction == subgroup_direction, "Consistent", "Inconsistent"),
    overall_direction = case_when(
      tissue_mean > 0 & subgroup_mean > 0 ~ "Overall Positive",
      tissue_mean < 0 & subgroup_mean < 0 ~ "Overall Negative",
      TRUE ~ "Mixed Direction"
    )
  ) %>%
  arrange(desc(combined_range))

cat("Dual-specific transcripts (range > 0.5 in both):", nrow(dual_specific), "\n")
cat("Direction consistency:", table(dual_specific$direction_consistency), "\n")
cat("Overall direction:", table(dual_specific$overall_direction), "\n\n")

# 5. 交叉分析：全局相关性与特异性的关系
cat("5. Cross Analysis: Global Correlation vs Specificity\n")

# 合并所有数据
combined_analysis <- global_corr_classified %>%
  select(Transcript_ID, Global_Correlation = Correlation, Global_FDR = FDR, 
         global_type = correlation_type, global_sign = sign_type) %>%
  left_join(
    tissue_range %>% select(Transcript_ID, tissue_range, tissue_mean),
    by = "Transcript_ID"
  ) %>%
  left_join(
    subgroup_range %>% select(Transcript_ID, subgroup_range, subgroup_mean),
    by = "Transcript_ID"
  ) %>%
  mutate(
    is_tissue_specific = tissue_range > 0.5,
    is_subgroup_specific = subgroup_range > 0.5,
    is_dual_specific = is_tissue_specific & is_subgroup_specific,
    specificity_type = case_when(
      is_dual_specific ~ "Dual-specific",
      is_tissue_specific & !is_subgroup_specific ~ "Tissue-specific",
      !is_tissue_specific & is_subgroup_specific ~ "Subgroup-specific",
      TRUE ~ "Non-specific"
    )
  )

# 分析全局相关性与特异性的关系
cat("Relationship between global correlation and specificity:\n")

correlation_specificity <- combined_analysis %>%
  group_by(global_type, specificity_type) %>%
  summarise(
    count = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = specificity_type,
    values_from = count,
    values_fill = 0
  )

print(correlation_specificity)
cat("\n")

# 6. 重点分析：强全局负相关的双重特异性转录本
cat("6. Special Analysis: Strong Negative Global & Dual-Specific\n")

negative_dual <- combined_analysis %>%
  filter(global_type == "Strong Negative" & specificity_type == "Dual-specific")

cat("Strong negative global AND dual-specific transcripts:", nrow(negative_dual), "\n")

if(nrow(negative_dual) > 0) {
  cat("\nThese transcripts (global r < -0.5 and high variability in both tissue and subgroup):\n")
  print(negative_dual %>% 
          select(Transcript_ID, Global_Correlation, tissue_range, subgroup_range) %>%
          arrange(Global_Correlation))
  
  # 保存结果
  write.table(
    negative_dual,
    file = "Strong_negative_dual_specific_transcripts.txt",
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}

# 7. 比较强正相关和强负相关的双重特异性
cat("\n7. Comparing Strong Positive vs Strong Negative Dual-Specific Transcripts\n")

positive_dual <- combined_analysis %>%
  filter(global_type == "Strong Positive" & specificity_type == "Dual-specific")

cat("Strong positive global & dual-specific:", nrow(positive_dual), "\n")
cat("Strong negative global & dual-specific:", nrow(negative_dual), "\n\n")

if(nrow(positive_dual) > 0 && nrow(negative_dual) > 0) {
  comparison <- bind_rows(
    positive_dual %>% mutate(group = "Positive"),
    negative_dual %>% mutate(group = "Negative")
  ) %>%
    group_by(group) %>%
    summarise(
      count = n(),
      mean_global_cor = mean(Global_Correlation, na.rm = TRUE),
      mean_tissue_range = mean(tissue_range, na.rm = TRUE),
      mean_subgroup_range = mean(subgroup_range, na.rm = TRUE),
      mean_tissue_mean = mean(tissue_mean, na.rm = TRUE),
      mean_subgroup_mean = mean(subgroup_mean, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("Comparison between positive and negative dual-specific transcripts:\n")
  print(comparison)
}

# 8. 保存所有结果
cat("\n8. Saving Results...\n")

# 8.1 全局相关性分类
write.table(
  global_corr_classified,
  file = "Global_correlation_classified.txt",
  sep = "\t", row.names = FALSE, quote = FALSE
)

# 8.2 强正相关转录本
write.table(
  strong_positive,
  file = "Strong_positive_global_transcripts_detailed.txt",
  sep = "\t", row.names = FALSE, quote = FALSE
)

# 8.3 强负相关转录本
write.table(
  strong_negative,
  file = "Strong_negative_global_transcripts_detailed.txt",
  sep = "\t", row.names = FALSE, quote = FALSE
)

# 8.4 双重特异性转录本
write.table(
  dual_specific,
  file = "Dual_specific_transcripts_comprehensive.txt",
  sep = "\t", row.names = FALSE, quote = FALSE
)

# 8.5 综合分析结果
write.table(
  combined_analysis,
  file = "Combined_analysis_results.txt",
  sep = "\t", row.names = FALSE, quote = FALSE
)

# 8.6 总结报告
summary_report <- data.frame(
  Category = c(
    "Total transcripts",
    "Significant global (FDR<0.05)",
    "Strong positive global (r > 0.5)",
    "Strong negative global (r < -0.5)",
    "Dual-specific (range>0.5 in both)",
    "Strong positive & dual-specific",
    "Strong negative & dual-specific",
    "Tissue-specific only",
    "Subgroup-specific only",
    "Non-specific"
  ),
  Count = c(
    nrow(global_corr),
    sum(global_corr$FDR < 0.05, na.rm = TRUE),
    nrow(strong_positive),
    nrow(strong_negative),
    nrow(dual_specific),
    nrow(positive_dual),
    nrow(negative_dual),
    sum(combined_analysis$specificity_type == "Tissue-specific", na.rm = TRUE),
    sum(combined_analysis$specificity_type == "Subgroup-specific", na.rm = TRUE),
    sum(combined_analysis$specificity_type == "Non-specific", na.rm = TRUE)
  ),
  Percentage = c(
    100,
    sum(global_corr$FDR < 0.05, na.rm = TRUE)/nrow(global_corr)*100,
    nrow(strong_positive)/nrow(global_corr)*100,
    nrow(strong_negative)/nrow(global_corr)*100,
    nrow(dual_specific)/nrow(global_corr)*100,
    nrow(positive_dual)/nrow(global_corr)*100,
    nrow(negative_dual)/nrow(global_corr)*100,
    sum(combined_analysis$specificity_type == "Tissue-specific", na.rm = TRUE)/nrow(global_corr)*100,
    sum(combined_analysis$specificity_type == "Subgroup-specific", na.rm = TRUE)/nrow(global_corr)*100,
    sum(combined_analysis$specificity_type == "Non-specific", na.rm = TRUE)/nrow(global_corr)*100
  )
)

write.table(
  summary_report,
  file = "Comprehensive_summary_report.txt",
  sep = "\t", row.names = FALSE, quote = FALSE
)

cat("\nSummary Report:\n")
print(summary_report)

# 9. 可视化：全局相关性分布
cat("\n9. Creating Visualizations...\n")

library(ggplot2)

# 全局相关性分布图（包括正负）
p1 <- ggplot(global_corr_classified, aes(x = Correlation, fill = sign_type)) +
  geom_histogram(bins = 50, alpha = 0.7, color = "white") +
  geom_vline(xintercept = c(-0.5, 0, 0.5), linetype = "dashed", 
             color = c("red", "black", "red"), linewidth = 0.8) +
  scale_fill_manual(values = c("Positive" = "#377EB8", "Negative" = "#E41A1C")) +
  labs(
    title = "Distribution of Global PDUI-PAC Correlations",
    subtitle = sprintf("Positive: %d (%.1f%%) | Negative: %d (%.1f%%)",
                      sum(global_corr$Correlation > 0, na.rm = TRUE),
                      mean(global_corr$Correlation > 0, na.rm = TRUE)*100,
                      sum(global_corr$Correlation < 0, na.rm = TRUE),
                      mean(global_corr$Correlation < 0, na.rm = TRUE)*100),
    x = "Spearman Correlation Coefficient",
    y = "Number of Transcripts",
    fill = "Correlation Sign"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

ggsave("Global_correlation_distribution_with_sign.png", plot = p1, width = 10, height = 6, dpi = 300)
ggsave("Global_correlation_distribution_with_sign.pdf", plot = p1, width = 10, height = 6)

# 全局相关性与变异性的关系图
p2 <- ggplot(combined_analysis, aes(x = Global_Correlation, y = tissue_range + subgroup_range)) +
  geom_point(aes(color = specificity_type), alpha = 0.6, size = 2) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "red", alpha = 0.5) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "blue", alpha = 0.5) +
  scale_color_manual(values = c(
    "Dual-specific" = "#E41A1C",
    "Tissue-specific" = "#377EB8",
    "Subgroup-specific" = "#FF7F00",
    "Non-specific" = "#999999"
  )) +
  labs(
    title = "Global Correlation vs Total Variability",
    subtitle = "X-axis: Global correlation | Y-axis: Tissue range + Subgroup range",
    x = "Global Correlation Coefficient",
    y = "Total Variability (Tissue + Subgroup Range)",
    color = "Specificity Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave("Global_vs_variability_scatter.png", plot = p2, width = 12, height = 8, dpi = 300)
ggsave("Global_vs_variability_scatter.pdf", plot = p2, width = 12, height = 8)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("\nKey Findings:\n")
cat("1. Strong positive global correlations:", nrow(strong_positive), "\n")
cat("2. Strong negative global correlations:", nrow(strong_negative), "\n")
cat("3. Dual-specific transcripts:", nrow(dual_specific), "\n")
cat("4. Strong positive & dual-specific:", nrow(positive_dual), "\n")
cat("5. Strong negative & dual-specific:", nrow(negative_dual), "\n")
cat("\nNote: Negative correlations indicate that higher distal APA usage is associated with LOWER expression.\n")
