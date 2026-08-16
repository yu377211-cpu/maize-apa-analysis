#!/usr/bin/env Rscript
# limma 差异分析（design 固定 + contrast 方向控制）
#
# 核心改动：design 矩阵永远用 alphabetical levels，不依赖文件名；
#           对比方向（logFC 正负）只通过 makeContrasts 决定。
# 效果：P 值 100% 不变（只换 contrast 方向），结果完全可复现。
#
# 用法:
#   Rscript limma_de_block_v2.R <expr_matrix> <condition_file> [output_prefix]
#
# condition.txt 格式:
#   2 列 (无 block):    sample_id<TAB>group
#   3 列 (有 block):    sample_id<TAB>tissue<TAB>genotype
#   contrast 方向从文件名解析: KernvsLMAN -> logFC = LMAN - Kern

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

# ===============================
# 0. 命令行参数
# ===============================
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2) {
  stop("Usage: Rscript limma_de_block_v2.R PDUI_matrix.txt condition.txt [output_prefix]\n",
       "  condition.txt format:\n",
       "    2 columns (no block): sample_id\\tgroup\n",
       "    3 columns (with block): sample_id\\ttissue\\tgenotype\n",
       "  contrast direction is inferred from file name: KernvsLMAN -> LMAN - Kern\n")
}

expr_file <- args[1]
cond_file <- args[2]
out_prefix <- if(length(args) >= 3) args[3] else "output"

# ===============================
# 1. 解析文件名中的对比组（不影响 design，只用于决定 contrast 方向）
# ===============================
parse_groups <- function(cond_file) {
  base_name <- basename(cond_file)
  base_name <- sub("\\.(condition|txt|csv)$", "", base_name)

  # 支持两种格式: "KernvsLMAN" 和 "Kern_vs_LMAN"
  for(pat in c("[Vv][Ss]", "_[Vv][Ss]_")) {
    if(grepl(pat, base_name)) {
      parts <- strsplit(base_name, pat)[[1]]
      if(length(parts) == 2) {
        return(list(group1 = trimws(parts[1]), group2 = trimws(parts[2])))
      }
    }
  }
  return(NULL)
}

# ===============================
# 2. 读取数据
# ===============================
read_data <- function(expr_file, cond_file) {
  expr <- fread(expr_file)
  mat <- as.matrix(expr[, -1, with=FALSE])
  rownames(mat) <- expr[[1]]

  cond <- read.delim(cond_file, header=TRUE, row.names=1, check.names=FALSE)

  # 2 列 (sample + group) -> ncol=1, no block
  # 3 列 (sample + group + covariate) -> ncol=2, has block
  has_cov <- ncol(cond) >= 2

  common <- intersect(colnames(mat), rownames(cond))
  if(length(common) == 0) {
    stop("No common samples between expression matrix and condition file!")
  }

  groups <- data.frame(group = cond[common, 1], stringsAsFactors = FALSE)
  rownames(groups) <- common

  cov <- NULL
  if(has_cov) {
    cov <- data.frame(covariate = cond[common, 2], stringsAsFactors = FALSE)
    rownames(cov) <- common
  }

  return(list(
    expr = mat[, common, drop=FALSE],
    groups = groups,
    cov = cov,
    has_cov = has_cov
  ))
}

# ===============================
# 3. 构建设计矩阵（永远用 alphabetical levels，固定不变）
# ===============================
build_fixed_design <- function(data) {
  # 关键：无论 data$groups$group 是 character 还是 factor，
  # 都强制先 as.character()，再用 sort(unique(...)) 显式设 levels
  # 这样 design 永远按字母顺序、完全固定
  # （如果直接 factor(factor_obj)，R 会继承原 levels，不重新排序）
  group_char <- as.character(data$groups$group)
  group_levels <- sort(unique(group_char))

  if(length(group_levels) != 2) {
    stop("Only 2-group comparisons supported. Found: ",
         paste(group_levels, collapse=", "))
  }

  group_factor <- factor(group_char, levels = group_levels)

  design <- model.matrix(~0 + group_factor)
  colnames(design) <- group_levels  # 用显式的 group_levels，不用 levels(group_factor)

  return(list(
    design = design,
    group_levels = group_levels
  ))
}

# ===============================
# 4. 根据文件名决定 contrast 方向（只决定公式，不动 design）
# ===============================
build_contrast <- function(design_info, cond_file) {
  g_levels <- design_info$group_levels  # alphabetical, 不变

  parsed <- parse_groups(cond_file)
  if(is.null(parsed)) {
    warning("Cannot parse direction from file name, using alphabetical order")
    formula <- paste0(g_levels[2], " - ", g_levels[1])
    cat("  (fallback) logFC =", g_levels[2], "-", g_levels[1], "\n")
  } else {
    g1 <- parsed$group1
    g2 <- parsed$group2

    if(!(g1 %in% g_levels) || !(g2 %in% g_levels)) {
      warning(sprintf("Groups '%s' or '%s' not in design levels (%s). Using alphabetical.",
                      g1, g2, paste(g_levels, collapse=", ")))
      formula <- paste0(g_levels[2], " - ", g_levels[1])
    } else if(g1 == g_levels[1]) {
      # 文件名 KernvsLMAN, alphabetical = c("Kern", "LMAN")
      # group1 在 alphabetical 中排第一 -> 方向 +1
      formula <- paste0(g_levels[2], " - ", g_levels[1])
      cat("  File says:", g1, "vs", g2, "-> logFC =", g2, "-", g1, "\n")
    } else {
      # 文件名 LMANvsKern, alphabetical = c("Kern", "LMAN")
      # group1 在 alphabetical 中排第二 -> 方向 -1
      formula <- paste0(g_levels[1], " - ", g_levels[2])
      cat("  File says:", g1, "vs", g2, "-> logFC =", g2, "-", g1, "\n")
    }
  }

  return(makeContrasts(contrasts = formula, levels = design_info$design))
}

# ===============================
# 5. 主分析
# ===============================
tryCatch({
  cat("========================================\n")
  cat("Step 1: Reading data\n")
  cat("========================================\n")
  data <- read_data(expr_file, cond_file)
  cat("  Samples:", ncol(data$expr), "\n")
  cat("  PACs/Transcripts:", nrow(data$expr), "\n")
  cat("  Groups:", paste(unique(data$groups$group), collapse=" vs "), "\n")
  if(data$has_cov) {
    cat("  Block: genotype information detected\n")
  } else {
    cat("  No block information (standard limma)\n")
  }

  # NA 填充
  cat("\nStep 2: Filling NA with per-row median\n")
  na_before <- sum(is.na(data$expr))
  expr_clean <- t(apply(data$expr, 1, function(x) {
    if(all(is.na(x))) {
      return(rep(NA, length(x)))
    }
    x[is.na(x)] <- median(x, na.rm=TRUE)
    x
  }))
  cat("  NA filled:", na_before - sum(is.na(expr_clean)), "/", na_before, "\n")

  # 固定 design
  cat("\nStep 3: Building design (FIXED, alphabetical)\n")
  design_info <- build_fixed_design(data)
  design <- design_info$design
  cat("  Design columns (will NOT change between runs):",
      paste(design_info$group_levels, collapse=", "), "\n")

  # contrast 方向
  cat("\nStep 4: Determining contrast direction from file name\n")
  cont.matrix <- build_contrast(design_info, cond_file)

  # lmFit + (可选) duplicateCorrelation（只跑一次，design 固定）
  cat("\nStep 5: Fitting model (lmFit runs ONCE)\n")
  if(data$has_cov) {
    genotype_factor <- factor(data$cov$covariate)

    cat("  Block summary:\n")
    print(table(genotype_factor))
    cat("  Total genotypes:", length(unique(genotype_factor)), "\n")
    cat("  Total samples:", length(genotype_factor), "\n")

    single_blocks <- names(which(table(genotype_factor) < 2))
    if(length(single_blocks) > 0) {
      cat("  WARNING: Single-sample blocks (will not contribute to correlation):\n")
      cat("    ", paste(single_blocks, collapse=", "), "\n")
    }

    cat("  Estimating intra-block correlation...\n")
    corfit <- duplicateCorrelation(expr_clean, design, block = genotype_factor)
    cat("  consensus correlation:", corfit$consensus, "\n")

    write.table(data.frame(Parameter = "consensus_correlation",
                           Value = corfit$consensus),
                file = paste0(out_prefix, ".correlation.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)

    fit <- lmFit(expr_clean, design,
                 block = genotype_factor,
                 correlation = corfit$consensus)
  } else {
    fit <- lmFit(expr_clean, design)
  }

  # contrasts.fit + eBayes（只根据 contrast 公式换方向）
  cat("\nStep 6: Computing contrast and eBayes\n")
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)

  # 提取结果
  result <- topTable(fit2, number = Inf)
  final <- data.frame(
    Transcript = rownames(result),
    result[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")],
    row.names = NULL
  )
  final <- final[complete.cases(final), ]

  # 显著筛选
  sig <- final[final$adj.P.Val < 0.05 & abs(final$logFC) > 0.2, ]

  # 输出
  out_cols <- c("Transcript", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
  write.table(final[, out_cols], paste0(out_prefix, ".full.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  write.table(sig[, out_cols], paste0(out_prefix, ".sig.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE, na = "")

  cat("\n========================================\n")
  cat("Summary\n")
  cat("========================================\n")
  cat("  Input PACs:", nrow(data$expr), "\n")
  cat("  Valid results:", nrow(final), "\n")
  cat("  Significant (adj.P<0.05 & |logFC|>0.2):", nrow(sig), "\n")
  if(nrow(sig) > 0) {
    cat("  Top 3:\n")
    print(head(sig[, out_cols], 3), row.names = FALSE)
  }
  cat("========================================\n")

}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  quit(status = 1)
})
