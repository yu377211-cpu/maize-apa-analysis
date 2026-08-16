library(ggplot2)
library(readr)
library(dplyr)

# 读取数据，假设文件中只有一列名为"length"的数值
lengths1 <- read_delim("ERP011069_CS.3UTR.length",delim="\t") %>% select(length)
lengths2 <- read_delim("ERP011069_DP2.3UTR.length",delim="\t") %>% select(length)

# 合并数据集，添加一个新列来区分两个数据集
data_combined <- bind_rows(
  tibble(length = lengths1$length, source = rep("CS", nrow(lengths1))),
  tibble(length = lengths2$length, source = rep("DP2", nrow(lengths2)))
)

# 绘制分布密度图
density_plot <- ggplot(data_combined, aes(x = log2(length), fill = source)) +
  geom_density(alpha = 0.5) +  # 添加半透明效果以便重叠部分可见
  labs(title = "Length Distribution Density Across Two Files", x = "PAC Length(log2)", y = "Density(%)")

# 保存为PNG文件
ggsave("length_distribution_density.png", plot = density_plot, width = 8, height = 6, dpi = 300)
