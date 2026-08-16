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
#expr_file <- file.path(dirname(input_file), paste0(prefix, ".txt"))
expr_file <- args[4]
output_prefix <- file.path(dirname(input_file), prefix)

# 参数设置
p_cutoff <- ifelse(length(args) > 1, as.numeric(args[2]), 0.05)
fc_cutoff <- ifelse(length(args) > 2, as.numeric(args[3]), 0.2)

# 1. 火山图函数
volcano_plot <- function() {
  res <- read.delim(input_file, stringsAsFactors = FALSE)
  
  # 创建差异类型列
  res$Type <- with(res, ifelse(
    adj.P.Val < p_cutoff & logFC > fc_cutoff, "Upregulated",
    ifelse(
      adj.P.Val < p_cutoff & logFC < -fc_cutoff, "Downregulated",
      "Non-significant"
    )
  ))
  
  # 绘图
  p <- ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), color = Type)) + 
    geom_point(aes(color = Type), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Upregulated" = "#E64B35", 
                                "Downregulated" = "#3182bd",
                                "Non-significant" = "gray60")) +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
    labs(title = "Volcano Plot", 
         x = "log2 Fold Change", 
         y = "-log10(Adjusted P-value)") +
    theme_minimal()
  
  # 保存（确保这一行是完整的）
  ggsave(paste0(output_prefix, "_volcano.pdf"), plot = p, width = 8, height = 6)
}

# 2. 热图函数
heatmap_plot <- function() {
  res <- read.delim(input_file, stringsAsFactors = FALSE)
  sig <- res[res$adj.P.Val < p_cutoff & abs(res$logFC) > fc_cutoff, ]
  
  if(nrow(sig) == 0) {
    message("没有显著差异转录本可绘制热图")
    return()
  }
  
  if(!file.exists(expr_file)) {
    stop(paste("表达矩阵文件不存在:", expr_file))
  }
  
  expr <- fread(expr_file)
  expr_mat <- as.matrix(expr[, -1])
  rownames(expr_mat) <- expr[[1]]
  
  heatdata <- expr_mat[sig$Transcript, , drop = FALSE]
  heatdata_scaled <- t(scale(t(heatdata)))
  
  # 热图参数
  show_ids <- ifelse(nrow(heatdata) <= 50, TRUE, FALSE)
  font_size <- ifelse(nrow(heatdata) > 100, 6, 10)
  
  pheatmap(
    heatdata_scaled,
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    breaks = seq(-3, 3, length.out = 101),
    show_rownames = show_ids,
    fontsize_row = font_size,
    filename = paste0(output_prefix, "_heatmap.pdf"),
    width = 10,
    height = max(6, nrow(heatdata)*0.3)
  )
}

# 主程序
message("\n=== 开始绘图 ===")
tryCatch({ volcano_plot() }, error = function(e) message("火山图错误: ", e$message))
tryCatch({ heatmap_plot() }, error = function(e) message("热图错误: ", e$message))
message("=== 绘图完成 ===\n")
