#!/usr/bin/env Rscript
library(data.table)
library(gplots)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 3) {
  stop("Usage: Rscript heatmap.R expression_matrix.txt sig_genes.txt group_info.txt [output_prefix]")
}

# 1. 读取数据
read_data <- function(expr_file, sig_file, group_file) {
  # 读取表达矩阵
  expr <- fread(expr_file)
  expr_mat <- as.matrix(expr[, -1])
  rownames(expr_mat) <- expr[[1]]
  
  # 读取显著基因列表
  sig_genes <- fread(sig_file)
  sig_ids <- sig_genes[[1]]  # 第一列为基因ID
  
  # 读取分组信息
  groups <- read.delim(group_file, header=FALSE, row.names=1)
  colnames(groups) <- "group"
  
  # 确保样本一致
  common_samples <- intersect(colnames(expr_mat), rownames(groups))
  list(
    expr = expr_mat[, common_samples],
    sig_ids = sig_ids,
    groups = groups[common_samples, , drop=FALSE]
  )
}

# 2. 绘制热图
draw_heatmap <- function(data, outfile) {
  # 提取显著基因表达数据
  sig_expr <- data$expr[rownames(data$expr) %in% data$sig_ids, ]
  
  if(nrow(sig_expr) == 0) {
    stop("没有找到匹配的显著差异基因，请检查基因ID是否一致")
  }
  
  # 数据标准化（行标准化）
  sig_expr_scaled <- t(scale(t(sig_expr)))
  
  # 设置分组颜色
  group_levels <- unique(data$groups$group)
  n_groups <- length(group_levels)
  
  # 处理分组颜色不足的情况
  if(n_groups < 3) {
    group_colors <- c("#E41A1C", "#377EB8")[1:n_groups]
  } else {
    group_colors <- brewer.pal(n_groups, "Set1")
  }
  names(group_colors) <- group_levels
  sample_colors <- group_colors[data$groups$group]
  
  # 热图颜色
  heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)
  
  # 是否显示ID
  show_colnames <- ifelse(ncol(sig_expr) <= 50, TRUE, FALSE)
  show_rownames <- ifelse(nrow(sig_expr) <= 50, TRUE, FALSE)
  
  # 动态调整图形大小
  pdf_width <- 10
  pdf_height <- max(6, min(20, nrow(sig_expr)/3))  # 限制最大高度为20英寸
  
  # 调整边距
  row_margin <- ifelse(show_rownames, 10, 2)
  col_margin <- ifelse(show_colnames, 8, 2)
  
  # 输出PDF
  pdf(paste0(outfile, ".heatmap.pdf"), width=pdf_width, height=pdf_height)
  
  # 设置图形参数
  par(mar=c(1,1,1,1))
  
  # 绘制热图
  heatmap.2(sig_expr_scaled,
            col = heatmap_colors,
            scale = "row",
            trace = "none",
            margins = c(col_margin, row_margin),
            cexRow = ifelse(show_rownames, 0.6 + 1/log10(nrow(sig_expr)), 0.1),
            cexCol = ifelse(show_colnames, 0.5 + 1/log10(ncol(sig_expr)), 0.1),
            Colv = TRUE,
            Rowv = TRUE,
            dendrogram = "both",
            key = TRUE,
            keysize = 1,
            density.info = "none",
            labCol = if(show_colnames) colnames(sig_expr) else NA,
            labRow = if(show_rownames) rownames(sig_expr) else NA,
            main = paste("Significant Transcripts (n=", nrow(sig_expr), ")", sep=""),
            ColSideColors = sample_colors)
  
  # 添加图例
  legend("topright", 
         legend = names(group_colors),
         fill = group_colors,
         border = group_colors,
         bty = "n",
         x.intersp = 0.5,
         y.intersp = 0.8,
         cex = 0.8)
  
  dev.off()
}

# 3. 主程序
tryCatch({
  # 参数处理
  expr_file <- args[1]
  sig_file <- args[2]
  group_file <- args[3]
  outfile <- ifelse(length(args) > 3, args[4], "heatmap")
  
  # 读取数据
  data <- read_data(expr_file, sig_file, group_file)
  
  # 绘制热图
  draw_heatmap(data, outfile)
  
  cat("热图绘制完成！\n")
  cat("显著基因数量:", length(data$sig_ids), "\n")
  cat("匹配到的基因数量:", sum(data$sig_ids %in% rownames(data$expr)), "\n")
  cat("样本数量:", ncol(data$expr), "\n")
  cat("分组数量:", length(unique(data$groups$group)), "\n")
  
}, error=function(e) {
  stop(paste("错误:", e$message))
})
