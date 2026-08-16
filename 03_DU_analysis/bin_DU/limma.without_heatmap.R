#!/usr/bin/env Rscript
library(data.table)
library(limma)

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 2) {
  stop("Usage: Rscript script.R PDUI_matrix.txt condition.txt [output_prefix]\n",
       "  condition.txt format (with header):\n",
       "    sample_id\\tgroup (e.g., Tropical, Temperate, NSS, Mixed)\n",
       "  The script will automatically infer comparison direction from the condition file name.\n",
       "  Example: NSSvsMixed.condition -> logFC = Mixed - NSS\n")
}

# ===============================
# 0. 从文件名解析比较方向（与block脚本完全一致）
# ===============================
parse_comparison <- function(cond_file) {
  base_name <- basename(cond_file)
  base_name <- sub("\\.(condition|txt|csv|group)$", "", base_name)
  
  # 使用正则 [Vv][Ss] 匹配 "vs"（不区分大小写），保留分割后的大小写
  if(grepl("[Vv][Ss]", base_name)) {
    parts <- strsplit(base_name, "[Vv][Ss]")[[1]]
    if(length(parts) == 2) {
      group1 <- trimws(parts[1])
      group2 <- trimws(parts[2])
      cat("Auto-detected comparison:", group1, "vs", group2, "\n")
      cat("  logFC will be calculated as:", group2, "-", group1, "\n")
      return(list(group1 = group1, group2 = group2))
    }
  }
  
  # 如果包含 _vs_ 也尝试
  if(grepl("_[Vv][Ss]_", base_name)) {
    parts <- strsplit(base_name, "_[Vv][Ss]_")[[1]]
    if(length(parts) == 2) {
      group1 <- trimws(parts[1])
      group2 <- trimws(parts[2])
      cat("Auto-detected comparison (underscore):", group1, "vs", group2, "\n")
      cat("  logFC will be calculated as:", group2, "-", group1, "\n")
      return(list(group1 = group1, group2 = group2))
    }
  }
  
  warning("Cannot parse comparison direction from file name: ", cond_file, 
          ". Using alphabetical order.")
  return(NULL)
}

# ===============================
# 1. 读取数据（统一使用带header的格式）
# ===============================
read_data <- function(expr_file, cond_file) {
  # 读取 PDUI 矩阵
  expr <- fread(expr_file)
  mat <- as.matrix(expr[, -1])
  rownames(mat) <- expr[[1]]
  
  # 读取 condition 文件（带header：sample_id  group）
  cond <- read.delim(cond_file, header=TRUE, row.names=1)
  
  # 确保样本一致（按表达矩阵列顺序保留，保持顺序稳定）
  common <- colnames(mat)[colnames(mat) %in% rownames(cond)]
  if(length(common) == 0) {
    stop("No common samples between expression matrix and condition file!")
  }
  
  # 提取分组信息
  groups <- data.frame(group = cond[common, 1])
  rownames(groups) <- common
  colnames(groups) <- "group"
  
  return(list(
    expr = mat[, common, drop=FALSE],
    groups = groups
  ))
}

# ===============================
# 2. 构建设计矩阵
# ===============================
build_design <- function(data, group_levels_ordered = NULL) {
  group_factor <- factor(data$groups$group)
  
  if(!is.null(group_levels_ordered)) {
    if(all(group_levels_ordered %in% levels(group_factor))) {
      group_factor <- factor(data$groups$group, levels = group_levels_ordered)
    } else {
      warning("Provided levels do not match data. Using alphabetical order.")
    }
  }
  
  if(length(levels(group_factor)) != 2) {
    stop("This script only supports two-group comparisons. Found: ",
         paste(levels(group_factor), collapse=", "))
  }
  
  design <- model.matrix(~0 + group_factor)
  colnames(design) <- levels(group_factor)
  
  return(list(design = design, group_levels = levels(group_factor)))
}

# ===============================
# 3. 执行差异分析（标准limma，不加block）
# ===============================
run_analysis <- function(data, group_levels_ordered = NULL) {
  # ---- 3.1 处理NA值（与block脚本一致） ----
  expr_clean <- t(apply(data$expr, 1, function(x) {
    if(all(is.na(x))) {
      return(rep(NA, length(x)))
    }
    x[is.na(x)] <- median(x, na.rm=TRUE)
    x
  }))
  
  # ---- 3.2 构建设计矩阵 ----
  design_info <- build_design(data, group_levels_ordered)
  design <- design_info$design
  group_levels <- design_info$group_levels
  
  # ---- 3.3 标准limma拟合 ----
  fit <- lmFit(expr_clean, design)
  
  # ---- 3.4 构建对比矩阵（levels[2] - levels[1]） ----
  cont.matrix <- makeContrasts(
    contrasts = paste0(group_levels[2], " - ", group_levels[1]),
    levels = design
  )
  
  cat("\nContrast:", group_levels[2], "-", group_levels[1], "\n")
  
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  
  # ---- 3.5 提取结果 ----
  result <- topTable(fit2, number=Inf)
  final <- data.frame(
    Transcript = rownames(result),
    result[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")],
    row.names = NULL
  )
  final <- final[complete.cases(final), ]
  
  return(final)
}

# ===============================
# 4. 主程序
# ===============================
tryCatch({
  expr_file <- args[1]
  cond_file <- args[2]
  out_prefix <- if(length(args) >= 3) args[3] else "output"
  
  cat("Reading data...\n")
  data <- read_data(expr_file, cond_file)
  cat("Samples:", ncol(data$expr), "\n")
  cat("PACs/Transcripts:", nrow(data$expr), "\n")
  cat("Groups:", paste(unique(data$groups$group), collapse=" vs "), "\n")
  
  # ---- 解析比较方向 ----
  group_order <- parse_comparison(cond_file)
  if(!is.null(group_order)) {
    group_levels_ordered <- c(group_order$group1, group_order$group2)
    if(all(group_levels_ordered %in% unique(data$groups$group))) {
      cat("Using auto-detected order:", group_levels_ordered[1], "vs", group_levels_ordered[2], "\n")
    } else {
      warning("Auto-detected groups not found in data. Using alphabetical order.")
      group_levels_ordered <- NULL
    }
  } else {
    group_levels_ordered <- NULL
  }
  
  cat("\nRunning limma analysis (standard, no blocking)...\n")
  results <- run_analysis(data, group_levels_ordered)
  
  # ---- 筛选显著结果 ----
  sig <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 0.2, ]
  
  # ---- 输出文件 ----
  out_cols <- c("Transcript", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
  write.table(results[, out_cols],
              paste0(out_prefix, ".full.txt"),
              sep="\t", row.names=FALSE, quote=FALSE, na="")
  write.table(sig[, out_cols],
              paste0(out_prefix, ".sig.txt"),
              sep="\t", row.names=FALSE, quote=FALSE, na="")
  
  cat("\n=== 结果验证 ===\n")
  cat("输入PAC数:", nrow(data$expr), "\n")
  cat("有效结果数:", nrow(results), "\n")
  cat("显著差异数 (adj.P.Val<0.05 & |logFC|>0.2):", nrow(sig), "\n")
  if(nrow(sig) > 0) {
    cat("示例显著结果:\n")
    print(head(sig, 3))
  }
  cat("===============\n")
  
}, error=function(e) stop(paste("错误:", e$message)))
