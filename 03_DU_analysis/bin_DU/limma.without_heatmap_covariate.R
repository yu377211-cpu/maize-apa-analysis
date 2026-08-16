#!/usr/bin/env Rscript
library(data.table)
library(limma)

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 2) {
  stop("Usage: Rscript script.R PDUI_matrix.txt condition.txt [output_prefix]\n",
       "  condition.txt format:\n",
       "    2 columns (no covariate): sample_id\\tgroup (e.g., Tropical, Temperate)\n",
       "    3 columns (with covariate): sample_id\\ttissue\\tgenotype\n",
       "    - column 2 (tissue) will be treated as the grouping variable\n",
       "    - column 3 (genotype) will be treated as a covariate\n")
}

# ===============================
# 1. 读取数据
# ===============================
read_data <- function(expr_file, cond_file) {
  # 读取 PDUI 矩阵
  expr <- fread(expr_file)
  mat <- as.matrix(expr[, -1])
  rownames(mat) <- expr[[1]]
  
  # 读取 condition 文件
  cond <- read.delim(cond_file, header=TRUE, row.names=1)
  # cond 列名: 第1列是 sample (已作为行名), 第2列是 group, 第3列是可选的协变量
  
  # 自动检测是否有协变量
  has_cov <- ncol(cond) >= 2  # 如果至少有2列非行名列，说明有协变量
  
  # 确保样本一致
  common <- intersect(colnames(mat), rownames(cond))
  if(length(common) == 0) {
    stop("No common samples between expression matrix and condition file!")
  }
  
  # 提取分组信息（第1列）
  groups <- data.frame(group = cond[common, 1])
  rownames(groups) <- common
  colnames(groups) <- "group"
  
  # 提取协变量（如果存在）
  cov <- NULL
  if(has_cov) {
    cov <- data.frame(covariate = cond[common, 2])
    rownames(cov) <- common
    colnames(cov) <- "covariate"
  }
  
  return(list(
    expr = mat[, common, drop=FALSE],
    groups = groups,
    cov = cov,
    has_cov = has_cov
  ))
}

# ===============================
# 2. 构建模型矩阵
# ===============================
build_design <- function(data) {
  group_factor <- factor(data$groups$group)
  
  if(data$has_cov) {
    # 有协变量：设计公式 ~ covariate + group
    # 这样会先扣除基因型效应，再检验组织差异
    cov_factor <- factor(data$cov$covariate)
    # 确保分组和协变量不混淆
    design <- model.matrix(~ cov_factor + group_factor)
    # 此时列名例如: (Intercept), cov_factorB73, group_factorRoot
    # group_factor 的最后一个水平会作为差异检验的基准
    return(list(design=design, group_levels=levels(group_factor), has_cov=TRUE))
  } else {
    # 无协变量：标准双分组比较 ~0 + group
    design <- model.matrix(~0 + group_factor)
    colnames(design) <- levels(group_factor)
    return(list(design=design, group_levels=levels(group_factor), has_cov=FALSE))
  }
}

# ===============================
# 3. 执行差异分析
# ===============================
run_analysis <- function(data) {
  # 处理NA值（行中位数填充）
  expr_clean <- t(apply(data$expr, 1, function(x) {
    x[is.na(x)] <- median(x, na.rm=TRUE)
    x
  }))
  
  # 构建模型
  design_info <- build_design(data)
  design <- design_info$design
  has_cov <- design_info$has_cov
  group_levels <- design_info$group_levels
  
  if(length(group_levels) != 2) {
    stop("This script only supports two-group comparisons. Found: ", 
         paste(group_levels, collapse=", "))
  }
  
  # 线性拟合
  fit <- lmFit(expr_clean, design)
  
  # 构建对比矩阵
  if(has_cov) {
    # 有协变量时，分组相关的列名为 "group_factorLevel2"
    # 构建对比：group_levels[2] vs group_levels[1]
    group_cols <- grep("^group_factor", colnames(design))
    if(length(group_cols) < 2) {
      # 如果只有一个分组列，说明编码方式为 treatment contrast
      # 此时 intercept 对应第一个水平，group_factor[2] 对应第二个水平 vs 第一个
      if("(Intercept)" %in% colnames(design)) {
        # 标准 treatment contrast: 对比系数为第2个分组列 - 0
        group_col <- grep("^group_factor", colnames(design))
        if(length(group_col) == 1) {
          cont_vector <- rep(0, ncol(design))
          cont_vector[group_col] <- 1
          cont.matrix <- matrix(cont_vector, ncol=1)
          colnames(cont.matrix) <- paste0(group_levels[2], "_vs_", group_levels[1])
        } else {
          stop("Cannot determine contrast structure with covariates.")
        }
      } else {
        stop("Cannot find group columns in design matrix.")
      }
    } else {
      # 有两个分组列（使用 ~0 + cov + group 编码时可能出现）
      cont_vector <- rep(0, ncol(design))
      cont_vector[group_cols[2]] <- 1
      cont_vector[group_cols[1]] <- -1
      cont.matrix <- matrix(cont_vector, ncol=1)
      colnames(cont.matrix) <- paste0(group_levels[2], "_vs_", group_levels[1])
    }
  } else {
    # 无协变量：标准对比
    cont.matrix <- makeContrasts(
      contrasts = paste0(group_levels[2], " - ", group_levels[1]),
      levels = design
    )
  }
  
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  
  # 提取结果
  result <- topTable(fit2, number=Inf)
  final <- data.frame(
    Transcript = rownames(result),
    result[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")],
    row.names = NULL
  )
  final[complete.cases(final), ]
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
  if(data$has_cov) {
    cat("Covariate (genotype) included\n")
    cat("Genotypes:", paste(unique(data$cov$covariate), collapse=", "), "\n")
  } else {
    cat("No covariate included\n")
  }
  
  cat("Running limma analysis...\n")
  results <- run_analysis(data)
  
  # 筛选显著结果
  sig <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 0.2, ]
  
  # 输出文件
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
