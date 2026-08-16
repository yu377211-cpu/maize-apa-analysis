# 加载必要的库
library(coloc)
library(dplyr)

# --------------------------
# 1. 读取 GWAS 和 3aQTL 数据
# --------------------------
# 假设输入文件是 TSV 格式（如果不是，请调整 read.table 参数）
gwas <- read.table("/share/pub/xingsl/shilai/project/APA/3aQTL_link/coloc/ZEAMAP_GWAS_signals_addmaf.input.modify.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)  # 列：chr, variantID, pos, ref, alt, pvalue, maf
qtl <- read.table("/share/pub/xingsl/shilai/project/APA/3aQTL_link/coloc/Seedling.cis_3aQTL.susieR.addpm.input.modify.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)  # 列：chr, variantID, pos, ref, alt, pvalue, maf

# --------------------------
# 2. 设置样本量
# --------------------------
N_gwas <- 744   # GWAS 样本量
N_3aqtl <- 941     # 3aQTL 样本量

# --------------------------
# 3. 数据清洗和验证
# --------------------------
# 确保MAF是数值型且在有效范围内
clean_data <- function(df) {
  df <- df[as.numeric(df$maf) > 0 & as.numeric(df$maf) < 1, ]
  df <- df[!is.na(df$pvalue) & df$pvalue > 0 & df$pvalue <= 1, ]
  return(df)
}

gwas <- clean_data(gwas)
qtl <- clean_data(qtl)

# --------------------------
# 3. 计算 beta 和 se
# --------------------------
# GWAS 计算
gwas <- gwas %>%
  mutate(
    Z = qnorm(pvalue / 2, lower.tail = FALSE),  # Z-score
    se = 1 / sqrt(2 * maf * (1 - maf) * N_gwas),  # 标准误
    beta = Z * se  # 效应量
  )

# 3aQTL 计算
qtl <- qtl %>%
  mutate(
    Z = qnorm(pvalue / 2, lower.tail = FALSE),
    se = 1 / sqrt(2 * maf * (1 - maf) * N_3aqtl),
    beta = Z * se
  )

# 检查清洗后的数据
print("GWAS dimensions after cleaning:")
print(dim(gwas))
print("QTL dimensions after cleaning:")
print(dim(qtl))

# 检查是否有共同SNP
common_snps <- intersect(gwas$variantID, qtl$variantID)
message("Number of common SNPs: ", length(common_snps))

# --------------------------
# 4. 合并数据（确保 SNP 匹配）
# --------------------------
# 按 variantID 或 chr:pos 合并
merged_data <- inner_join(
  gwas %>% select(variantID, chr, pos, beta_gwas = beta, se_gwas = se, pvalue_gwas = pvalue, maf),
  qtl %>% select(variantID, chr, pos, beta_qtl = beta, se_qtl = se, pvalue_qtl = pvalue, maf),
  by = c("variantID", "chr", "pos")
) %>%
  #filter(maf > 0 & maf < 1)

# --------------------------
# 5. 转换为 COLOC 输入格式
# --------------------------
dataset_gwas <- list(
  beta = merged_data$beta_gwas,
  varbeta = merged_data$se_gwas^2,
  pvalues = merged_data$pvalue_gwas,
  type = "quant",  # GWAS 通常是定量性状
  snp = as.character(merged_data$variantID),
  MAF = merged_data$maf,
  N = N_gwas
)

dataset_qtl <- list(
  beta = merged_data$beta_qtl,
  varbeta = merged_data$se_qtl^2,
  pvalues = merged_data$pvalue_qtl,
  type = "quant",  # 3aQTL 也是定量性状
  snp = as.character(merged_data$variantID),
  MAF = merged_data$maf,
  N = N_3aqtl
)

# --------------------------
# 6. 运行 COLOC 分析
# --------------------------
coloc_result <- coloc.abf(dataset_gwas, dataset_qtl)

# 查看结果
head(coloc_result$summary)

# --------------------------
# 7. 输出结果到文件
# --------------------------
# 创建输出目录（如果不存在）
output_dir <- "/share/pub/xingsl/shilai/project/APA/3aQTL_link/coloc/results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 保存完整结果
saveRDS(coloc_result, file = file.path(output_dir, "coloc_result.rds"))

# 保存摘要结果为文本文件
sink(file = file.path(output_dir, "coloc_summary.txt"))
print(coloc_result$summary)
sink()

# 保存详细结果为CSV文件
write.csv(coloc_result$results, file = file.path(output_dir, "coloc_detailed_results.csv"), row.names = FALSE)

# 保存过滤后的输入数据
write.csv(merged_data, file.path(output_dir, "filtered_input_data.csv"), 
          row.names = FALSE)

message("COLOC analysis completed successfully!")
message("Results saved to: ", output_dir)
message("Number of SNPs analyzed: ", nrow(merged_data))
