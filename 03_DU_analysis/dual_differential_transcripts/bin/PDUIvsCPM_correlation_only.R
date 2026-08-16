#!/usr/bin/env Rscript

# ================================================================
# 计算 PDUI 与 TMM‑CPM 总表达（求和）的 Spearman 相关性
# 过滤条件严格沿用 PDUIvsCPM_bothTissueSubpop.R：
#   - TPM 中 0 转为 NA
#   - 按转录本聚合时对 PAC 的 TPM 求和（na.rm=TRUE）
#   - 仅保留共同转录本和共同样本
#   - 长格式整合后删除 Tissue/Subgroup/PDUI/TPM 为 NA 的行
#   - 每个转录本（或转录本×分组）有效观测数 >= 3 才计算相关性
#   - 最终输出仅包含成功计算的记录（非 NA）
# 输出：
#   - PDUI_PAC_global_correlation.txt
#   - PDUI_PAC_tissue_correlation.txt
#   - PDUI_PAC_subgroup_correlation.txt
# ================================================================

library(data.table)

# ------------------------- 1. 读取数据 -------------------------
cat("读取 PDUI 矩阵...\n")
pdui_mat <- fread("/share/pub/xingsl/shilai/project/APA/AFreview/PAC_3seq_diff_usage/sta_plot/dual_differential_transcripts/subpop_inter_tissue.sig.transcript.PDUI", header = TRUE, sep = "\t")
transcript_ids <- pdui_mat[[1]]
pdui_mat <- as.matrix(pdui_mat[, -1])
rownames(pdui_mat) <- transcript_ids

cat("读取 CPM 矩阵...\n")
cpm_mat <- fread("/share/pub/xingsl/shilai/project/APA/AFreview/PAC_3seq_diff_usage/sta_plot/dual_differential_transcripts/subpop_inter_tissue.sig.transcript.TMM_CPM", header = TRUE, sep = "\t")
pac_ids <- cpm_mat[[1]]
cpm_mat <- as.matrix(cpm_mat[, -1])
rownames(cpm_mat) <- pac_ids

# 将 0 转为 NA（与原脚本一致）
#cpm_mat[cpm_mat == 0] <- NA

cat("读取 PAC-Transcript 映射...\n")
map_df <- fread("/share/pub/xingsl/shilai/project/APA/AFreview/PAC_3seq_diff_usage/sta_plot/dual_differential_transcripts/subpop_inter_tissue.sig.transcript_pac", header = TRUE, sep = "\t")
setnames(map_df, c("PAC_ID", "Transcript_ID"))

cat("读取分组信息...\n")
group_df <- fread("/share/pub/xingsl/shilai/project/APA/PAC_3seq_diff_usage/sta_plot/specific_transcript/S1960_without_L3Mid.group.txt", header = TRUE, sep = "\t")
setnames(group_df, c("Sample_ID", "Tissue", "Subgroup"))

# 获取所有样本（以 PDUI 矩阵的列名为准）
all_samples <- colnames(pdui_mat)
group_df <- group_df[Sample_ID %in% all_samples, ]
group_df <- group_df[match(all_samples, group_df$Sample_ID), ]
tissue_vec <- group_df$Tissue
subgroup_vec <- group_df$Subgroup

tissues <- unique(tissue_vec)
subgroups <- unique(subgroup_vec)

# ------------------------- 2. 按转录本聚合 TPM（求和） -------------------------
cat("按转录本聚合 TPM（求和）...\n")
# 将 cpm_mat 转为 data.table 方便操作
cpm_dt <- as.data.table(cpm_mat, keep.rownames = "PAC_ID")
# 合并映射
cpm_dt <- merge(map_df, cpm_dt, by = "PAC_ID", all.y = FALSE)  # 只保留有映射的 PAC
# 按转录本分组，对每个样本列求和（na.rm=TRUE）
cpm_agg <- cpm_dt[, lapply(.SD, sum, na.rm = TRUE), 
                  by = Transcript_ID, .SDcols = all_samples]
# 将求和结果转为矩阵，行名为转录本
cpm_agg_mat <- as.matrix(cpm_agg[, ..all_samples])
rownames(cpm_agg_mat) <- cpm_agg$Transcript_ID

cat("精细聚合（区分真0和全缺失）...\n")
all_transcripts <- unique(map_df$Transcript_ID)
# 初始化聚合矩阵
cpm_agg2 <- matrix(NA, nrow = length(all_transcripts), ncol = length(all_samples),
                   dimnames = list(all_transcripts, all_samples))

# 提前按转录本分组 PAC 列表
pac_list <- split(map_df$PAC_ID, map_df$Transcript_ID)

for (tr in all_transcripts) {
  pacs <- pac_list[[tr]]
  if (is.null(pacs) || length(pacs) == 0) next
  # 提取这些 PAC 的 TPM 矩阵
  sub_cpm <- cpm_mat[pacs, , drop = FALSE]
  # 对每个样本，检查是否全部为 NA
  for (samp in all_samples) {
    vals <- sub_cpm[, samp]
    if (all(is.na(vals))) {
      cpm_agg2[tr, samp] <- NA
    } else {
      cpm_agg2[tr, samp] <- sum(vals, na.rm = TRUE)
    }
  }
}

# 只保留与 PDUI 矩阵共有的转录本
common_transcripts <- intersect(rownames(pdui_mat), rownames(cpm_agg2))
cpm_agg_final <- cpm_agg2[common_transcripts, , drop = FALSE]
pdui_sub <- pdui_mat[common_transcripts, , drop = FALSE]

# 共同样本（已确保一致）
common_samples <- intersect(colnames(pdui_sub), colnames(cpm_agg_final))
pdui_sub <- pdui_sub[, common_samples, drop = FALSE]
cpm_agg_final <- cpm_agg_final[, common_samples, drop = FALSE]

# ------------------------- 3. 整合成长格式并过滤 -------------------------
cat("整合长格式并过滤...\n")
# 转成长格式
pdui_long <- as.data.table(pdui_sub, keep.rownames = "Transcript_ID")
pdui_long <- melt(pdui_long, id.vars = "Transcript_ID", variable.name = "Sample_ID", value.name = "PDUI")

cpm_long <- as.data.table(cpm_agg_final, keep.rownames = "Transcript_ID")
cpm_long <- melt(cpm_long, id.vars = "Transcript_ID", variable.name = "Sample_ID", value.name = "TPM")

# 合并
combined <- merge(pdui_long, cpm_long, by = c("Transcript_ID", "Sample_ID"))
# 添加分组信息
combined <- merge(combined, group_df, by = "Sample_ID", all.x = FALSE)  # 仅保留有分组信息的样本

# 过滤（与原脚本一致）
combined <- combined[!is.na(Tissue) & !is.na(Subgroup), ]
combined <- combined[!is.na(PDUI), ]
combined <- combined[!is.na(TPM), ]

cat("有效观测对数: ", nrow(combined), "\n")

# ------------------------- 4. 计算全局相关性 -------------------------
cat("计算全局相关性...\n")
global_res <- combined[, .(
  obs_count = .N,
  Correlation = ifelse(.N >= 3, 
                       suppressWarnings(cor(PDUI, TPM, method = "spearman", use = "complete.obs")),
                       NA_real_),
  P_value = ifelse(.N >= 3,
                   suppressWarnings(cor.test(PDUI, TPM, method = "spearman", exact = FALSE)$p.value),
                   NA_real_)
), by = Transcript_ID]

# 过滤掉相关性计算失败的行（与原脚本一致）
global_res <- global_res[!is.na(Correlation) & !is.na(P_value), ]
global_res[, FDR := p.adjust(P_value, method = "fdr")]
setorder(global_res, P_value)

write.table(global_res, file = "PDUI_PAC_global_correlation.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ------------------------- 5. 计算组织相关性 -------------------------
cat("计算组织相关性...\n")
tissue_res_list <- list()
for (tiss in tissues) {
  sub_data <- combined[Tissue == tiss, ]
  if (nrow(sub_data) == 0) next
  # 计算每个转录本在该组织内的相关性
  res <- sub_data[, .(
    obs_count = .N,
    Group_Cor = ifelse(.N >= 3,
                       suppressWarnings(cor(PDUI, TPM, method = "spearman", use = "complete.obs")),
                       NA_real_),
    Group_P = ifelse(.N >= 3,
                     suppressWarnings(cor.test(PDUI, TPM, method = "spearman", exact = FALSE)$p.value),
                     NA_real_)
  ), by = Transcript_ID]
  res <- res[!is.na(Group_Cor) & !is.na(Group_P), ]
  if (nrow(res) > 0) {
    res[, Tissue := tiss]
    tissue_res_list[[tiss]] <- res
  }
}
tissue_res <- rbindlist(tissue_res_list)
# 按组织分组 FDR
tissue_res[, Group_FDR := p.adjust(Group_P, method = "fdr"), by = Tissue]
setorder(tissue_res, Tissue, Group_P)

write.table(tissue_res, file = "PDUI_PAC_tissue_correlation.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ------------------------- 6. 计算亚群相关性 -------------------------
cat("计算亚群相关性...\n")
subgroup_res_list <- list()
for (subg in subgroups) {
  sub_data <- combined[Subgroup == subg, ]
  if (nrow(sub_data) == 0) next
  res <- sub_data[, .(
    obs_count = .N,
    Group_Cor = ifelse(.N >= 3,
                       suppressWarnings(cor(PDUI, TPM, method = "spearman", use = "complete.obs")),
                       NA_real_),
    Group_P = ifelse(.N >= 3,
                     suppressWarnings(cor.test(PDUI, TPM, method = "spearman", exact = FALSE)$p.value),
                     NA_real_)
  ), by = Transcript_ID]
  res <- res[!is.na(Group_Cor) & !is.na(Group_P), ]
  if (nrow(res) > 0) {
    res[, Subgroup := subg]
    subgroup_res_list[[subg]] <- res
  }
}
subgroup_res <- rbindlist(subgroup_res_list)
subgroup_res[, Group_FDR := p.adjust(Group_P, method = "fdr"), by = Subgroup]
setorder(subgroup_res, Subgroup, Group_P)

write.table(subgroup_res, file = "PDUI_PAC_subgroup_correlation.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n分析完成！\n")
cat("输出文件:\n")
cat("  - PDUI_PAC_global_correlation.txt\n")
cat("  - PDUI_PAC_tissue_correlation.txt\n")
cat("  - PDUI_PAC_subgroup_correlation.txt\n")
