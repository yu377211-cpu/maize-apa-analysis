#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(optparse)
})

## =========================
## 参数
## =========================
option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "Input coupling table (space separated)"),
  make_option(c("-o", "--out_prefix"), type = "character",
              default = "Coupling_diverging_bar"),
  make_option(c("--width"), type = "numeric", default = 6),
  make_option(c("--height"), type = "numeric", default = 4)
)

opt <- parse_args(OptionParser(option_list = option_list))

## =========================
## 读入数据（空格分隔 & 列不齐）
## =========================
df <- read.table(
  opt$input,
  header = TRUE,
  sep = "",
  fill = TRUE,
  stringsAsFactors = FALSE
)

## 必要列检查
required_cols <- c("group1", "flag")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

## =========================
## 清洗 + 定义 coupling
## =========================
df2 <- df %>%
  filter(!is.na(flag), flag != 0) %>%
  mutate(
    coupling = ifelse(flag == 1, "negative", "positive")
  )

## =========================
## 按 group1 汇总比例
## =========================
sum_df <- df2 %>%
  count(group1, coupling) %>%
  group_by(group1) %>%
  mutate(
    total = sum(n),
    proportion = n / total
  ) %>%
  ungroup() %>%
  mutate(
    y = ifelse(coupling == "positive",
               -proportion,
               proportion)
  )

## =========================
## 作图
## =========================
p <- ggplot(sum_df,
            aes(x = group1, y = y, fill = coupling)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      negative = "#7AAAE0",
      positive = "#E07B7B"
    ),
    labels = c(
      negative = "Elongation-associated upregulation",
      positive = "Elongation-associated downregulation"
    ),
    name = NULL  # 移除图例标题
  ) +
  scale_y_continuous(
    labels = abs,
    limits = c(-1, 1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(
    x = NULL,
    y = "Proportion of transcripts",
    title = "Directional bias of APA–expression coupling"
  )

## =========================
## 输出
## =========================
ggsave(
  paste0(opt$out_prefix, ".pdf"),
  p,
  width = opt$width,
  height = opt$height
)

ggsave(
  paste0(opt$out_prefix, ".png"),
  p,
  width = opt$width,
  height = opt$height,
  dpi = 300
)
