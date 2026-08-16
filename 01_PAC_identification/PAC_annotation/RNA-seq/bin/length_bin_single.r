library(data.table)
library(ggplot2)
library(ggpubr)
library(readr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
cs_file <- args[1]
output_file <- args[2]


# 读取数据，假设文件中只有一列名为"length"的数值
lengths1 <- read_delim(cs_file,delim="\t") %>% select(length)

# 合并数据集，添加一个新列来区分两个数据集
data_combined <- bind_rows(
  tibble(length = lengths1$length, source = rep("CS", nrow(lengths1))),
)


data_combined$length_bin <- cut(data_combined$length,c(0,1,10,50,100,1000))


cs <- data.frame(table(data_combined[data_combined$source=="CS","length_bin"]))
cs$source <- "CS" 
cs$ratio <- cs$Freq/sum(cs$Freq) * 100

p <- ggplot(cs,aes(x=length_bin,y=ratio,fill= length_bin)) +  # fill=source 的话就是单色
  geom_col() +
  scale_x_discrete(labels=c("1","2-10","11-50","51-100",">100"))+
  labs(x = "PAC Length", y = "Percentage(%)")+
  theme_classic()+
  theme(legend.position = "none") # 添加这一行来去掉图例

ggsave(filename = paste0(output_file, ".png"), plot=p, width=8, height=6,dpi=300)
ggsave(filename = paste0(output_file, ".pdf"), plot=p, width=8, height=6)
