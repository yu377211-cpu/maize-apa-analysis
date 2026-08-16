# 加载必要的包
library(ggplot2)
library(dplyr)

# 读取数据（假设文件是制表符分隔的）
data <- read.table("Traits.stat", header = FALSE, sep = "\t", col.names = c("Count", "Name"), colClasses = c("numeric", "character"))

# 合并计数 <20 的条目为 "others"
data_processed <- data %>%
  mutate(Name = ifelse(Count < 30, "Others", as.character(Name))) %>%
  group_by(Name) %>%
  summarise(Count = sum(Count)) %>%
  arrange(desc(Count))  # 按计数降序排列

# 计算比例
data_processed$Percentage <- round(data_processed$Count / sum(data_processed$Count) * 100, 1)

# 绘制饼图
pie_chart <- ggplot(data_processed, aes(x = "", y = Count, fill = factor(Name))) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(Percentage, "%")), 
            position = position_stack(vjust = 0.5), 
            size = 4) +
  labs(fill = "Traits") +  # 图例标题
  theme_void() +  # 移除背景、坐标轴等
  theme(legend.position = "right")  # 图例放在右侧

# 显示图形
print(pie_chart)

# 保存图形（可选）
ggsave("traits.pie_chart.png", pie_chart, width = 6, height = 6, dpi = 300)
