#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(data.table)
})

# 参数解析
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1) {
  stop("Usage: Rscript draw_plots.R <limma_results.txt> [pvalue_cutoff] [logFC_cutoff]")
}

# 文件名处理
input_file <- args[1]
file_parts <- unlist(strsplit(basename(input_file), "\\."))
prefix <- paste(file_parts[1:(length(file_parts)-2)], collapse = ".")
#expr_file <- args[4]
output_prefix <- file.path(dirname(input_file), prefix)

# 参数设置
p_cutoff <- ifelse(length(args) > 1, as.numeric(args[2]), 0.05)
fc_cutoff <- ifelse(length(args) > 2, as.numeric(args[3]), 0.2)

# 1. 火山图函数
volcano_plot <- function() {
  res <- read.delim(input_file, stringsAsFactors = FALSE)

  # 创建差异类型列（保持原有逻辑）
  res$Type <- with(res, ifelse(
    adj.P.Val < p_cutoff & logFC > fc_cutoff, "Lengthening",
    ifelse(
      adj.P.Val < p_cutoff & logFC < -fc_cutoff, "Shortening",
      "Non-significant"
    )
  ))

  # 绘图（关键修改：颜色命名与Type值完全匹配）
  p <- ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), color = Type)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(
      values = c(
        "Lengthening"   = "#E64B35",   # 红色
        "Shortening"    = "#3182bd",   # 蓝色
        "Non-significant" = "gray60"
      ),
      # 可选：调整图例顺序和标题
      breaks = c("Lengthening", "Shortening", "Non-significant"),
      name = "APA Direction"
    ) +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
    labs(
      title = "Volcano Plot of Differential Poly(A) Site Usage",
      x = "log2 Fold Change (Usage)",
      y = "-log10(Adjusted P-value)"
    ) +
    theme_minimal()

  ggsave(paste0(output_prefix, "_volcano.pdf"), plot = p, width = 8, height = 6)
}

# 主程序
message("\n=== 开始绘图 ===")
tryCatch({ volcano_plot() }, error = function(e) message("火山图错误: ", e$message))
message("=== 绘图完成 ===\n")
