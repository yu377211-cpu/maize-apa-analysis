library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
ratio_file <- args[1]
out_file <- args[2]

# 每个数字占一行
numbers <- scan(ratio_file, what = numeric(), sep = "\n")

# 检查是否读取成功
if(length(numbers) == 0) {
  stop("没有读取到数字，请检查文件路径和格式。")
}
# 将大于2000的数值替换为2000
numbers[numbers > 100] <- 100
numbers[numbers < 0] <- 0

# 绘制密度分布图
p <- ggplot(data.frame(x = numbers), aes(x = x)) +
  geom_density(fill = 'skyblue', alpha = 0.7) +
  labs(x = "Ratio of the distance from PA to transcriptional terminal to 3UTR length (%)", y = "Density") +
  theme_classic()

# 显示图形
print(p)
p <- p + xlim(0, 100) # 设置x轴的范围，排除大于2000的值
ggsave(filename = paste0(out_file, ".png"), plot=p, width = 8, height = 6, dpi = 300)
ggsave(filename = paste0(out_file, ".pdf"), plot=p, width=8, height=6)
