#!/usr/bin/env Rscript

# 用法：Rscript summary_susie_final_v2.R <input_file> [output_file]
# 示例：Rscript summary_susie_final_v2.R susieR_res.all_genes.txt susie_stats_final.txt

library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript summary_susie.R <input_file> [output_file]")
}
input_file <- args[1]
output_file <- ifelse(length(args) >= 2, args[2], "susie_summary_final.txt")

cat("Reading", input_file, "...\n")
data <- read.table(input_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
cat("Total rows (associations):", nrow(data), "\n")

# ---------- 1. 基础统计 ----------
total_assoc <- nrow(data)
unique_pheno <- length(unique(data$locus_id))

# ---------- 2. Credible set 级别的统计 ----------
# 创建唯一的 cs_id = locus_id + cs（排除无效cs）
data <- data %>%
  mutate(cs_id = ifelse(!is.na(cs) & cs != "" & cs != "NA",
                        paste(locus_id, cs, sep = ":"),
                        NA_character_))

# 只保留属于某个 CS 的行
cs_data <- data %>% filter(!is.na(cs_id))

# 2.1 unique credible sets 数量
unique_cs <- n_distinct(cs_data$cs_id)

# 2.2 每个 CS 的 size（cs_size）去重（因为同一 CS 所有行 cs_size 相同）
cs_summary <- cs_data %>%
  distinct(cs_id, locus_id, cs, cs_size)   # 每个 CS 一行

# 2.3 CS size 的统计：mean, median, IQR (Q3-Q1)
cs_size_mean <- mean(cs_summary$cs_size, na.rm = TRUE)
cs_size_median <- median(cs_summary$cs_size, na.rm = TRUE)
cs_size_q1 <- quantile(cs_summary$cs_size, 0.25, na.rm = TRUE)
cs_size_q3 <- quantile(cs_summary$cs_size, 0.75, na.rm = TRUE)
cs_size_iqr <- cs_size_q3 - cs_size_q1

# 2.4 每个 phenotype 的 CS 数量
cs_per_pheno <- cs_summary %>% count(locus_id)
avg_cs_per_pheno <- mean(cs_per_pheno$n, na.rm = TRUE)

# ---------- 新增：有CS和无CS的phenotype统计 ----------
phenotype_with_cs <- n_distinct(cs_data$locus_id)
phenotype_without_cs <- unique_pheno - phenotype_with_cs
percent_with_cs <- phenotype_with_cs / unique_pheno * 100   # 百分比

# ---------- 3. Lead SNP 统计 ----------
# 每个 CS 中取 PIP 最大的 SNP（若并列，取第一个）
lead_snps <- cs_data %>%
  group_by(cs_id) %>%
  slice_max(pip, n = 1, with_ties = FALSE) %>%
  ungroup()

num_lead_snps <- nrow(lead_snps)   # 应等于 unique_cs

# lead SNP 的 PIP 分布
lead_pip_mean <- mean(lead_snps$pip, na.rm = TRUE)
lead_pip_median <- median(lead_snps$pip, na.rm = TRUE)

# ---------- 4. Candidate SNP (PIP > 0.1, 去重) ----------
candidate_snps <- cs_data %>%
  filter(pip > 0.1) %>%
  distinct(variant_id)   # 去重

num_candidate_snps <- nrow(candidate_snps)

# ---------- 5. 结果汇总 ----------
stats <- data.frame(
  Metric = c(
    "Total associations (rows)",
    "Unique phenotypes (locus_id)",
    "Phenotypes with at least one CS",
    "Phenotypes without any CS",
    "Percent phenotypes with CS (%)",
    "Unique credible sets (cs_id: locus_id+cs)",
    "Candidate SNPs (PIP>0.1, distinct variant_id)",
    "Lead SNPs (one per CS)",
    "CS size - Mean",
    "CS size - Median",
    "CS size - IQR (Q3-Q1)",
    "Average credible sets per resolved phenotype",
    "Lead SNP PIP - Mean",
    "Lead SNP PIP - Median"
  ),
  Value = c(
    total_assoc,
    unique_pheno,
    phenotype_with_cs,
    phenotype_without_cs,
    round(percent_with_cs, 2),
    unique_cs,
    num_candidate_snps,
    num_lead_snps,
    round(cs_size_mean, 2),
    round(cs_size_median, 2),
    round(cs_size_iqr, 2),
    round(avg_cs_per_pheno, 2),
    round(lead_pip_mean, 4),
    round(lead_pip_median, 4)
  )
)

# 打印到屏幕
print(stats, row.names = FALSE)

# 保存到文件
write.table(stats, file = output_file, row.names = FALSE, sep = "\t", quote = FALSE)

cat("\nSummary saved to", output_file, "\n")

# 可选：输出 lead SNP 列表和 candidate SNP 列表（供后续分析）
write.table(lead_snps, file = "lead_snps_list.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(candidate_snps, file = "candidate_snps_list.txt", sep = "\t", row.names = FALSE, quote = FALSE)
cat("Lead SNP list saved to lead_snps_list.txt\n")
cat("Candidate SNP list saved to candidate_snps_list.txt\n")

phenotypes_with_cs <- cs_data %>% distinct(locus_id)
write.table(phenotypes_with_cs, file = "phenotypes_with_cs_list.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
cat("Phenotypes with at least one CS saved to phenotypes_with_cs_list.txt\n")
