#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(pheatmap)
})

## -----------------------------
## 参数
## -----------------------------
option_list <- list(
  make_option(c("-i", "--input"), type="character", help="final_score matrix (rows=subpop, cols=transcripts)"),
  make_option(c("-o", "--output"), type="character", default="APA_coupling"),
  make_option(c("-t", "--top"), type="integer", default=30, help="Top N transcripts for main heatmap"),
  make_option(c("-m", "--min_abs"), type="double", default=1, help="Min abs(final_score) to consider")
)

opt <- parse_args(OptionParser(option_list=option_list))

## -----------------------------
## 读取数据
## -----------------------------
mat <- read.table(opt$input, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
mat <- as.matrix(mat)

# NA → 0（表示无耦合或未出现）
mat[is.na(mat)] <- 0

## -----------------------------
## 计算 strength & contrast
## -----------------------------
strength <- apply(mat, 2, function(x) sum(abs(x)))
contrast <- apply(mat, 2, function(x) {
  max(abs(x)) - median(abs(x))
})

score_df <- data.frame(
  transcript = colnames(mat),
  strength = strength,
  contrast = contrast,
  final_rank_score = strength * contrast
)

# 基本过滤
score_df <- score_df %>%
  filter(strength >= opt$min_abs, contrast > 0)

if (nrow(score_df) == 0) {
  stop("No transcript passed min_abs / contrast filter")
}

## -----------------------------
## 主图：Top transcripts
## -----------------------------
top_n <- min(opt$top, nrow(score_df))

top_transcripts <- score_df %>%
  arrange(desc(final_rank_score)) %>%
  slice(1:top_n) %>%
  pull(transcript)

mat_main <- mat[, top_transcripts, drop=FALSE]

## 主图绘制
pdf(paste0(opt$output, ".main_heatmap.pdf"), width=10, height=4)
pheatmap(
  mat_main,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 12,
  fontsize_col = 6,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  main = "Subpopulation-specific APA–expression coupling (Top transcripts)"
)
dev.off()

png(paste0(opt$output, ".main_heatmap.png"), width=1200, height=400)
pheatmap(
  mat_main,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 12,
  fontsize_col = 6,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  main = "Subpopulation-specific APA–expression coupling (Top transcripts)"
)
dev.off()

## -----------------------------
## 补图：全量（不画列名）
## -----------------------------
pdf(paste0(opt$output, ".supp_heatmap.pdf"), width=10, height=4)
pheatmap(
  mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  main = "Global APA–expression coupling landscape"
)
dev.off()

## -----------------------------
## 输出表格
## -----------------------------
write.csv(score_df,
          paste0(opt$output, ".transcript_strength_contrast.csv"),
          row.names = FALSE)

write.csv(
  score_df %>% arrange(desc(final_rank_score)) %>% slice(1:opt$top),
  paste0(opt$output, ".top_transcripts.csv"),
  row.names = FALSE
)

cat("Done.\n")

