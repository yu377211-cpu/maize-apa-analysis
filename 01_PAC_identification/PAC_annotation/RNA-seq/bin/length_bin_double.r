library(data.table)
library(ggplot2)
library(ggpubr)
library(readr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
cs_file <- args[1]
dp2_file <- args[2]
output_file <- args[3]


# 读取数据，假设文件中只有一列名为"length"的数值
lengths1 <- read_delim(cs_file,delim="\t") %>% select(length)
lengths2 <- read_delim(dp2_file,delim="\t") %>% select(length)

# 合并数据集，添加一个新列来区分两个数据集
data_combined <- bind_rows(
  tibble(length = lengths1$length, source = rep("CS", nrow(lengths1))),
  tibble(length = lengths2$length, source = rep("DP2", nrow(lengths2)))
)


data_combined$length_bin <- cut(data_combined$length,c(0,1,2,10,100,300))


cs <- data.frame(table(data_combined[data_combined$source=="CS","length_bin"]))
cs$source <- "CS" 
cs$ratio <- cs$Freq/sum(cs$Freq)
dp2 <- data.frame(table(data_combined[data_combined$source=="DP2","length_bin"]))
dp2$source <- "DP2"
dp2$ratio <- dp2$Freq/sum(dp2$Freq)
sp <- rbind(cs,dp2)

p <- ggplot(sp,aes(x=length_bin,y=ratio,fill=source)) + 
  geom_col(position = position_dodge(width = 0.9,preserve = "single")) +
  scale_x_discrete(labels=c("1","2","3-10","11-100",">100"))+
  labs(x = "PAC Length", y = "Ratio")
ggsave(output_file,p)

