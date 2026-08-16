library(tidyverse)
library(readr)

# 定义文件列表
file_list <- c(
  "GRoot.updown_longshort.list",
  "GShoot.updown_longshort.list",
  "Kern.updown_longshort.list",
  "L3Base.updown_longshort.list",
  "L3Tip.updown_longshort.list",
  "LMAD.updown_longshort.list",
  "LMAN.updown_longshort.list"
)

# 读取并合并所有文件
all_data <- data.frame()

for (file in file_list) {
  if (file.exists(file)) {
    cat("正在读取文件:", file, "\n")
    
    # 先读取原始文件内容
    lines <- readLines(file)
    
    # 跳过表头行和包含"---"的分隔线
    data_lines <- lines[!grepl("^--|^---|^----", lines) & 
                        !grepl("^组.*比较组", lines) &
                        lines != ""]
    
    if (length(data_lines) > 0) {
      # 将数据行转换为数据框
      # 使用str_split_fixed按空格分割，最多分成5列
      split_data <- stringr::str_split_fixed(data_lines, "\\s+", 5)
      
      # 转换为数据框
      file_data <- as.data.frame(split_data, stringsAsFactors = FALSE)
      colnames(file_data) <- c("group", "comp", "feature", "type", "mark")
      
      # 转换mark列为数值
      file_data$mark <- as.numeric(file_data$mark)
      
      # 添加到总数据
      all_data <- bind_rows(all_data, file_data)
      
      cat("  读取了", nrow(file_data), "行数据\n")
    } else {
      cat("  文件没有有效数据\n")
    }
  } else {
    warning(paste("文件不存在:", file))
  }
}

# 查看数据结构
cat("\n总数据行数:", nrow(all_data), "\n")
cat("列名:", colnames(all_data), "\n")
cat("\n前几行数据:\n")
print(head(all_data))

# 清理数据
all_data_clean <- all_data %>%
  # 移除空特征或无效特征
  filter(!is.na(feature), 
         feature != "", 
         feature != "特征ID",
         !grepl("^-----", feature)) %>%
  # 移除无效的比较组
  filter(!is.na(comp), 
         comp != "", 
         comp != "比较组",
         !grepl("^-----", comp)) %>%
  # 移除无效的组
  filter(!is.na(group), 
         group != "", 
         group != "组",
         !grepl("^-----", group)) %>%
  # 创建唯一的样本标识符
  mutate(
    sample_id = paste(group, comp, sep = "_"),
    # 确保mark是数值
    mark = as.numeric(mark)
  )

cat("\n清理后数据行数:", nrow(all_data_clean), "\n")

# 获取所有唯一的特征ID
all_features <- unique(all_data_clean$feature[all_data_clean$feature != "" & 
                                              !is.na(all_data_clean$feature)])

# 获取所有唯一的样本ID
all_samples <- unique(all_data_clean$sample_id)

cat("\n唯一特征数:", length(all_features), "\n")
cat("唯一样本数:", length(all_samples), "\n")

# 创建矩阵
result_matrix <- matrix(0, 
                       nrow = length(all_samples),
                       ncol = length(all_features),
                       dimnames = list(all_samples, all_features))

# 填充矩阵
for (i in 1:nrow(all_data_clean)) {
  sample_id <- all_data_clean$sample_id[i]
  feature <- all_data_clean$feature[i]
  mark_val <- all_data_clean$mark[i]
  
  if (!is.na(sample_id) && !is.na(feature) && !is.na(mark_val)) {
    if (sample_id %in% rownames(result_matrix) && feature %in% colnames(result_matrix)) {
      result_matrix[sample_id, feature] <- mark_val
    }
  }
}

# 查看结果
cat("\n矩阵维度:", dim(result_matrix), "\n")
cat("\n前5行前5列:\n")
print(result_matrix[1:min(5, nrow(result_matrix)), 1:min(5, ncol(result_matrix))])

# 统计信息
cat("\n统计信息:\n")
cat("总样本数:", nrow(result_matrix), "\n")
cat("总特征数:", ncol(result_matrix), "\n")
cat("标记为1的数量:", sum(result_matrix == 1, na.rm = TRUE), "\n")
cat("标记为-1的数量:", sum(result_matrix == -1, na.rm = TRUE), "\n")
cat("标记为0的数量:", sum(result_matrix == 0, na.rm = TRUE), "\n")

# 保存结果
write.csv(result_matrix, "feature_matrix.csv", row.names = TRUE)

# 创建包含分组信息的数据框
result_df <- as.data.frame(result_matrix) %>%
  rownames_to_column("sample_id") %>%
  separate(sample_id, into = c("group", "comp"), sep = "_", remove = FALSE, fill = "right")

write.csv(result_df, "feature_matrix_with_groups.csv", row.names = FALSE)

# 转置版本
transposed_matrix <- t(result_matrix)
write.csv(transposed_matrix, "feature_matrix_transposed.csv", row.names = TRUE)

cat("\n文件已保存:\n")
cat("  feature_matrix.csv - 原始矩阵文件\n")
cat("  feature_matrix_with_groups.csv - 包含分组信息\n")
cat("  feature_matrix_transposed.csv - 转置矩阵\n")

library(dplyr)
library(tidyr)
library(tibble)

# =========================
# 1. 亚群内二次聚合（关键修复）
# =========================

main_by_subpop <- all_data_clean %>%
  filter(mark != 0) %>%                # 只保留有耦合的记录
  group_by(group, feature) %>%
  summarise(
    n_comp = n_distinct(comp),         # 该亚群中出现于多少比较组合
    mean_mark = mean(mark),             # 方向一致性
    sum_mark  = sum(mark),              # 强度
    .groups = "drop"
  ) %>%
  mutate(
    consensus = case_when(
      mean_mark >= 0.6  ~  1,
      mean_mark <= -0.6 ~ -1,
      TRUE              ~  0
    )
  )

# =========================
# 2. 转为矩阵（此时一对 group-feature 只有一行）
# =========================

main_matrix_by_subpop <- main_by_subpop %>%
  filter(consensus != 0) %>%            # 可选：主图只画有一致方向的
  select(group, feature, consensus) %>%
  distinct() %>%                        # ★保险：彻底去重
  pivot_wider(
    names_from  = feature,
    values_from = consensus,
    values_fill = 0
  ) %>%
  column_to_rownames("group")

write.csv(main_matrix_by_subpop,
          "main_matrix_by_tissue.csv")

cat("✔ 成功生成 main_matrix_by_tissue.csv - 聚合矩阵文件\n")

# =========================
# 3. 补图：跨亚群共享转录本
# =========================

# 先统计每个转录本出现于多少亚群
feature_freq <- main_by_subpop %>%
  filter(consensus != 0) %>%
  distinct(group, feature) %>%
  count(feature, name = "n_subpop")

# 设定补图阈值（≥2 个亚群）
selected_features <- feature_freq %>%
  filter(n_subpop >= 2) %>%
  pull(feature)

# ⚠️ 关键：再次按 (group, feature) 聚合
supp_by_subpop <- main_by_subpop %>%
  filter(feature %in% selected_features) %>%
  group_by(group, feature) %>%
  summarise(
    consensus = mean(consensus),
    .groups = "drop"
  ) %>%
  mutate(
    consensus = ifelse(consensus > 0, 1, -1)
  )

# 转成补图矩阵
supp_matrix_by_subpop <- supp_by_subpop %>%
  select(group, feature, consensus) %>%
  distinct() %>%                      # 双保险
  pivot_wider(
    names_from  = feature,
    values_from = consensus,
    values_fill = 0
  ) %>%
  column_to_rownames("group")

write.csv(supp_matrix_by_subpop,
          "supp_matrix_by_tissue_shared.csv")

cat("✔ 成功生成 supp_matrix_by_tissue_shared.csv - 比较组的矩阵文件\n")

