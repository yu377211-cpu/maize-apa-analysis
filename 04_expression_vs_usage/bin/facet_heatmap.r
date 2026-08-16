#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character"),
  make_option(c("-o", "--out_prefix"), type="character",
              default="APA_coupling_heatmap"),
  make_option(c("--width"), type="numeric", default=14),
  make_option(c("--height"), type="numeric", default=4)
)

opt <- parse_args(OptionParser(option_list = option_list))

## 读入数据（自动识别空格）
df <- read.table(opt$input,
                 header = TRUE,
                 sep = "\t",
                 stringsAsFactors = FALSE,
                 check.names = FALSE)

## 转换百分号
cols <- c("long∩up","long∩down","short∩up","short∩down")

df[cols] <- lapply(df[cols], function(x){
  as.numeric(gsub("%","",x))
})

## 转为长格式
df_long <- df %>%
  pivot_longer(cols = all_of(cols),
               names_to = "type",
               values_to = "ratio")

## 画图
p <- ggplot(df_long,
            aes(x = type,
                y = compare,
                fill = ratio)) +
  geom_tile(color="white") +
  facet_wrap(~group1, nrow = 1, scales="free_y") +
  scale_fill_gradient(
    low = "white",
    high = "#c53030", #红色
    #high = "#3070B0", #蓝色
    limits = c(0, 60), oob = scales::squish, #仅仅针对蓝色的分组织情况，有极端值拉伸色阶，这里限制一下色阶上限
    name = "Ratio (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1),
    strip.text = element_text(face="bold")
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Overlap between differential APA usage and expression"
  )

ggsave(paste0(opt$out_prefix,".pdf"),
       p, width=opt$width, height=opt$height)

ggsave(paste0(opt$out_prefix,".png"),
       p, width=opt$width, height=opt$height, dpi=300)

