#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript GO_plot_and_filter_final.R input_file output_prefix")
}

df <- read.table(args[1], sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", fill = TRUE)

# 1. 检查并提取所需列（按列名）
required_cols <- c("Term", "term_type", "queryitem", "bgitem", "pvalue", "FDR")
missing <- setdiff(required_cols, colnames(df))
if (length(missing) > 0) {
  stop(paste("输入文件缺少以下列:", paste(missing, collapse = ", "),
             "\n实际列名:", paste(colnames(df), collapse = ", ")))
}

df_sel <- df[, required_cols, drop = FALSE]

# 2. 强制转换为数值（处理因子或字符）
df_sel$queryitem <- as.numeric(as.character(df_sel$queryitem))
df_sel$bgitem <- as.numeric(as.character(df_sel$bgitem))
df_sel$pvalue <- as.numeric(as.character(df_sel$pvalue))
df_sel$FDR <- as.numeric(as.character(df_sel$FDR))

# 检查转换结果
if (any(is.na(df_sel$queryitem)) || any(is.na(df_sel$bgitem))) {
  warning("部分 queryitem 或 bgitem 转换为 NA，请检查原始数据是否有非数字字符")
}

# 计算 Rich Factor
df_sel$bi <- df_sel$queryitem / df_sel$bgitem

# 3. 筛选 FDR < 0.05
a <- df_sel[df_sel$FDR < 0.05 & !is.na(df_sel$FDR), ]
if (nrow(a) == 0) {
  stop("没有 FDR < 0.05 的条目，请检查数据")
}

# 重命名 term_type（若为 P/C/F 则改为 BP/CC/MF）
a$term_type[a$term_type == "P"] <- "BP"
a$term_type[a$term_type == "C"] <- "CC"
a$term_type[a$term_type == "F"] <- "MF"

cat("过滤后各类型数量:\n")
print(table(a$term_type))

# 4. 提取每种类型的前10个（按 queryitem 降序）
get_top_10 <- function(data, type) {
  sub <- data[data$term_type == type, ]
  if (nrow(sub) == 0) return(NULL)
  sub <- sub[order(sub$queryitem, decreasing = TRUE), ]
  head(sub, 10)
}

BP_top <- get_top_10(a, "BP")
CC_top <- get_top_10(a, "CC")
MF_top <- get_top_10(a, "MF")

GO_data <- rbind(
  if (!is.null(BP_top)) BP_top,
  if (!is.null(CC_top)) CC_top,
  if (!is.null(MF_top)) MF_top
)

if (nrow(GO_data) == 0) {
  stop("没有数据用于绘图，请检查是否有至少一种类型有显著条目")
}

GO_data$term_type <- factor(GO_data$term_type, levels = c("BP", "CC", "MF"))
cat("\n最终绘图数据:\n")
print(table(GO_data$term_type))

# 5. 绘图函数
make_GO_bubble <- function(go_data) {
  ggplot(go_data, aes(x = bi, y = reorder(Term, queryitem))) +
    geom_point(aes(size = queryitem, color = -log10(pvalue))) +
    scale_color_gradient(low = "#377eb8", high = "#e41a1c") +
    scale_size_continuous(range = c(3, 8)) +
    facet_grid(term_type ~ ., scales = "free_y", space = "free_y") +
    labs(x = "Rich Factor", y = "", size = "Gene Count", color = "-log10(pvalue)") +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold", size = 12),
      axis.text.y = element_text(size = 10, color = "black"),
      panel.grid.major.y = element_line(linetype = "dotted"),
      legend.position = "right",
      panel.spacing = unit(0.6, "lines"),
      panel.background = element_rect(fill = NA, color = "gray30", size = 0.6),
      strip.background = element_rect(fill = "gray90", color = "gray30")
    )
}

# 6. 输出
output_prefix <- args[2]

pdf(paste0(output_prefix, "_bubble.pdf"), width = 9, height = 10)
print(make_GO_bubble(GO_data))
dev.off()

png(paste0(output_prefix, "_bubble.png"), width = 9, height = 10, units = "in", res = 300)
print(make_GO_bubble(GO_data))
dev.off()

# 7. 输出所有显著条目（排序后）
a$term_type <- factor(a$term_type, levels = c("BP", "CC", "MF"))
a <- a[order(a$term_type, -a$queryitem), ]
outfile <- paste0(output_prefix, "_all_significant.txt")
write.table(a, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n显著富集GO条目已保存至: ", outfile, "\n")
