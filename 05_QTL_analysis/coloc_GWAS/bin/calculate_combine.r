library(dplyr)
library(purrr)

# --------------------------
# 1. 读取 GWAS 和 3aQTL 数据
# --------------------------
gwas <- read.table("/share/pub/xingsl/shilai/project/APA/3aQTL_link/coloc/ZEAMAP_GWAS_signals_addmaf.input.modify.tsv", 
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)
qtl <- read.table("/share/pub/xingsl/shilai/project/APA/3aQTL_link/coloc/Seedling.cis_3aQTL.susieR.addpm.input.modify.tsv", 
                  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# --------------------------
# 2. 设置样本量
# --------------------------
N_gwas <- 744
N_3aqtl <- 941

# --------------------------
# 3. 数据清洗和验证
# --------------------------
clean_data <- function(df) {
  df %>%
    mutate(maf = as.numeric(maf)) %>%
    filter(maf > 0 & maf < 1, !is.na(pvalue), pvalue > 0, pvalue <= 1)
}

gwas <- clean_data(gwas)
qtl <- clean_data(qtl)

# --------------------------
# 4. 计算 beta 和 se 并追加到原文件
# --------------------------
# GWAS计算
gwas <- gwas %>%
  mutate(
    Z = qnorm(pvalue / 2, lower.tail = FALSE),
    se = 1 / sqrt(2 * maf * (1 - maf) * N_gwas),
    beta = Z * se
  ) %>%
  select(-Z)

# QTL计算
qtl <- qtl %>%
  mutate(
    Z = qnorm(pvalue / 2, lower.tail = FALSE),
    se = 1 / sqrt(2 * maf * (1 - maf) * N_3aqtl),
    beta = Z * se
  ) %>%
  select(-Z)

# 输出目录
output_dir <- "/share/pub/xingsl/shilai/project/APA/3aQTL_link/coloc/results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 输出带beta/se的完整文件
write.table(gwas, file.path(output_dir, "GWAS_with_beta_se.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(qtl, file.path(output_dir, "QTL_with_beta_se.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# --------------------------
# 5. 改进的最近位点匹配合并
# --------------------------
find_nearest_match <- function(query, target) {
  # 按染色体分组处理
  chromosomes <- unique(query$chr)
  
  map_dfr(chromosomes, function(chr) {
    q_sub <- query %>% filter(chr == !!chr) %>% arrange(pos)
    t_sub <- target %>% filter(chr == !!chr) %>% arrange(pos)
    
    if (nrow(t_sub) == 0) return(tibble())
    
    # 使用findInterval快速定位最近位点
    idx <- findInterval(q_sub$pos, t_sub$pos)
    idx <- pmin(pmax(idx, 1), nrow(t_sub))  # 确保索引在范围内
    
    # 比较前后两个候选位点
    dist_prev <- abs(q_sub$pos - t_sub$pos[idx])
    dist_next <- abs(q_sub$pos - t_sub$pos[pmin(idx + 1, nrow(t_sub))])
    
    # 选择更近的位点
    nearest_idx <- ifelse(dist_prev <= dist_next, idx, pmin(idx + 1, nrow(t_sub)))
    
    # 构造结果
    q_sub %>%
      mutate(
        nearest_pos = t_sub$pos[nearest_idx],
        distance = abs(pos - nearest_pos),
        beta_gwas = t_sub$beta[nearest_idx],
        se_gwas = t_sub$se[nearest_idx],
        pvalue_gwas = t_sub$pvalue[nearest_idx],
        maf_gwas = t_sub$maf[nearest_idx],
        nearest_gwas_id = t_sub$variantID[nearest_idx]
      )
  })
}

# 执行匹配（QTL为主，找最近的GWAS位点）
matched_data <- find_nearest_match(qtl, gwas)

# 输出匹配结果
write.table(matched_data, file.path(output_dir, "nearest_match_results.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# 计算并输出匹配统计信息
match_stats <- matched_data %>%
  group_by(chr) %>%
  summarise(
    n_matches = n(),
    avg_distance = mean(distance),
    min_distance = min(distance),
    max_distance = max(distance)
  ) %>%
  ungroup()

write.table(match_stats, file.path(output_dir, "match_statistics.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

message("Analysis completed successfully!")
message(paste0("Matched ", nrow(matched_data), " GWAS-QTL pairs"))
message(paste0("Average distance: ", round(mean(matched_data$distance), 2), " bp"))
message(paste0("Results saved to: ", output_dir))
