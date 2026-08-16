library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)

#Rscript *.r num.txt num.png

args=commandArgs(trailingOnly=T)

pac <- fread(args[1],header=T,sep='\t')
pac[, PAC_num := factor(PAC_num, labels = c("1", "2", "3", "4", "5", "6", rep(">=7",nrow(pac)-6)))]

pac <- pac %>%
  mutate(Total_CS = sum(CS),
         Total_DP2 = sum(DP2), # 计算所有数值列的总和
         CS_percentage = CS / Total_CS * 100, # 计算Group1的百分比
         DP2_percentage = DP2 / Total_DP2 * 100) # 计算Group2的百分比
pac_long <- pac %>%
  select(PAC_num, CS_percentage, DP2_percentage) %>% 
  gather(key = "group", value = "percentage", -PAC_num) %>%
  mutate(group = gsub("CS_percentage", "CS", group),
         group = gsub("DP2_percentage", "DP2", group))

bar_plot <- ggplot(pac_long, aes(x = PAC_num, y = percentage, fill = group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  labs(x = "Number of PACs", y = "Percentage of Genes (%)") +
  theme_classic()

ggsave(args[2], plot = bar_plot, width = 7, height = 6, dpi = 300)
