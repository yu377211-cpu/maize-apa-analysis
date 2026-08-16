ATGC_1_plot <- function(
  x_fu_file,
  x_zheng_file,
  zhongjian_file,
  output_png,
  output_pdf
){
  library(ggplot2)

  width = 6
  height = 4
  line_size=1.5
  legend_size=1.2

  print("读取不变文件")
  # 读取文件反向数据并进行处理
  data1 <- read.table(x_zheng_file, header = TRUE)
  data1_cols=ncol(data1)
  data1 <- data1[, -ncol(data1)]
  col_names1 <- names(data1)

  print("读取x坐标反向文件")
  # 读取文件正向数据并进行处理
  data2 <- read.table(x_fu_file, header = TRUE)
  data2_cols=ncol(data2)
  data2 <- data2[, -ncol(data2)]
  #print(data2)
  print("合并数据")
  # 合并两个数据集
  combined_data <- merge(data2, data1, by = col_names1[1])
  combined_data[, -1] <- t(t(combined_data[, -1]) / colSums(combined_data[, -1]))
  #print(colSums(combined_data[,-1]))
  colnames(combined_data) <- c("motif",-(data2_cols-3):(data1_cols-2))
  # 将数据表导出为txt文件
  write.table(combined_data, file = zhongjian_file, sep = "\t", row.names = FALSE)

  #-----------------------------------------------------------------------------------------------------------------------------------------
  print("开始画图")
  # 创建空的数据框来存储重塑后的数据
  data_melted <- data.frame()

  # 遍历每一列进行重塑操作
  for (i in 2:length(combined_data)) {
    # 创建一个临时数据框用于存储当前列的重塑结果
    temp <- data.frame(motif = combined_data$motif, x =names(combined_data)[i] , y = combined_data[,i])
    # 将当前列的重塑结果添加到数据框中
    data_melted <- rbind(data_melted, temp)
  }

  # 定义颜色向量
  colors <- c(A = "blue", U = "yellow", C = "red", G = "green")

  # 定义标签的顺序和名称
  labels <- c("A", "U", "C", "G")

  #data_melted$x <- as.factor(data_melted$x)

  #nums<-c(-(data1_cols-2),0,(data2_cols-2))
  # 绘制折线图

  p <- ggplot(data_melted, aes(x = as.numeric(x), y = y, color = motif, group = motif)) +
    geom_line(linewidth = line_size) +
    scale_color_manual(values = colors, breaks = labels, labels = labels) +
    labs(x = "Position", y = "Proportion", color = "motif") +scale_x_continuous(limits = c(-(data2_cols-2),(data1_cols-2)),breaks=seq((floor((data2_cols-2)/10) * (-10)),(ceiling((data1_cols-2)/10) * 10),50)) +theme_bw() + 
    theme(legend.text = element_text(face = "bold"),
    legend.title = element_blank(),
    legend.key.size = unit(legend_size, "lines"),
    panel.grid = element_line(color = "gray"))

  # 保存图像为_png格式
  ggsave(output_png, plot = p, width = width, height = height, dpi = 300)

  # 保存图像为_pdf格式
  ggsave(output_pdf, plot = p, width = width, height = height, device = "pdf")
  print("完成画图")
}
