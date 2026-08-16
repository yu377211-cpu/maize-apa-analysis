# 加载必要的包
library(dplyr)
library(readr)

# 定义主函数
process_files <- function(input_file1, input_file2, output_file, max_distance) {
  # 读取文件
  file1 <- read_tsv(input_file1)
  file2 <- read_tsv(input_file2)

# 为文件1的列名添加后缀 .x
colnames(file1) <- paste0(colnames(file1), ".x")
# 为文件2的列名添加后缀 .y
colnames(file2) <- paste0(colnames(file2), ".y")

  # 按染色体进行内连接，并计算距离，过滤距离小于等于 max_distance 的记录
  results <- file1 %>%
    inner_join(file2, by = c("chr.x" = "chr.y"), relationship = "many-to-many") %>%
    mutate(distance = abs(pos.x - pos.y)) %>%
    filter(distance <= max_distance) %>%
    arrange(pos.x, distance)

  # 写入输出文件
  write_tsv(results, output_file)
}

# 从命令行获取参数
args <- commandArgs(trailingOnly = TRUE)

# 检查参数数量是否正确
if (length(args) != 4) {
  stop("Usage: Rscript script_name.R <input_file1> <input_file2> <output_file> <max_distance>")
}

# 解析参数
input_file1 <- args[1]
input_file2 <- args[2]
output_file <- args[3]
max_distance <- as.numeric(args[4])

# 检查 max_distance 是否为正数
if (is.na(max_distance) || max_distance <= 0) {
  stop("max_distance must be a positive number")
}

# 调用主函数
process_files(input_file1, input_file2, output_file, max_distance)
