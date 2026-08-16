# ============================
# 1. 加载包
# ============================
library(ggplot2)
library(dplyr)

# ============================
# 2. 读取输入文件（修正文件名）
# ============================
tau_data <- read.table("tissue_specific_Tau_results.txt", 
                       header = TRUE, row.names = 1, stringsAsFactors = FALSE)
tau_col <- grep("tau|Tau", colnames(tau_data), value = TRUE, ignore.case = TRUE)[1]
if (is.na(tau_col)) stop("No Tau column found")
tau_data <- tau_data[, tau_col, drop = FALSE]
colnames(tau_data) <- "Tau"

# 高特异性列表（全部组织）
high_list <- read.table("tissue_specific_allTissue_specific_genes.txt", 
                        header = FALSE, stringsAsFactors = FALSE)[, 1]

# 保守列表
conserved_list <- read.table("tissue_specific_conserved_across_subgroups.txt", 
                             header = FALSE, stringsAsFactors = FALSE)[, 1]

# ============================
# 3. 构建数据框
# ============================
df <- data.frame(
  Transcript = rownames(tau_data),
  Tau = tau_data$Tau,
  stringsAsFactors = FALSE
)
df$is_conserved <- df$Transcript %in% conserved_list

n_total <- nrow(df)
n_high <- sum(df$Tau >= 0.5)
n_conserved <- sum(df$is_conserved)

cat("Total:", n_total, "\nHigh-specificity (Tau>=0.5):", n_high, 
    sprintf("(%.1f%%)", n_high/n_total*100), "\nConserved:", n_conserved, "\n")

# ============================
# 4. 绘图（修正版）
# ============================
# 定义配色（可在此处修改为您的其他图配色）
# 直方图填充色、密度线色、阈值线色、地毯线色
col_hist <- "gray70"       # 灰色直方图
col_density <- "black"     # 黑色密度曲线
col_threshold <- "#D55E00" # 橙色阈值（也可改为"darkred"）
col_rug <- "#0072B2"       # 蓝色地毯线（也可改为"#56B4E9"等）

p <- ggplot(df, aes(x = Tau)) +
  
  # 灰色直方图（设置boundary=0确保第一个bin从0开始，避免溢出左边界）
  geom_histogram(aes(y = after_stat(density)), 
                 binwidth = 0.02, 
                 boundary = 0,          # 关键：对齐到0
                 fill = col_hist, 
                 color = "white", 
                 size = 0.2,
                 alpha = 0.9) +
  
  # 黑色密度曲线
  geom_density(color = col_density, size = 0.9, bw = 0.02) +
  
  # 阈值线
  geom_vline(xintercept = 0.5, linetype = "dashed", color = col_threshold, size = 0.9) +
  
  # 蓝色地毯线（保守转录本）
  geom_rug(data = subset(df, is_conserved),
           aes(x = Tau), 
           sides = "b", 
           color = col_rug, 
           size = 1.5, 
           length = unit(0.2, "cm")) +
  
  # 左上角注释（纯英文，无特殊字符）
  annotate("text", x = 0.02, y = Inf, 
           label = paste0("High-specificity: ", n_high, " (", 
                          round(n_high/n_total*100, 1), "%)\nConserved: ", n_conserved), 
           hjust = 0, vjust = 1, size = 4.5, color = "black", fontface = "bold",
           family = "sans") +
  
  labs(x = "Tau specificity score", y = "Density") +
  
  # 主题：纯白背景，无网格，保留坐标轴框线
  theme_classic(base_size = 14) +   # 经典主题自带白色背景和框线
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.8),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    legend.position = "none"
  ) +
  
  # x轴严格从0开始，防止左侧溢出
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  # y轴从0开始，不扩展
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

# ============================
# 5. 保存 PNG + PDF
# ============================
ggsave("Tau_distribution_final.png", 
       plot = p, width = 6, height = 5, dpi = 600, units = "in")

ggsave("Tau_distribution_final.pdf", 
       plot = p, width = 6, height = 5, units = "in")

print(p)
