library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
dis_file <- args[1]
out_file <- args[2]

rt=read.table(dis_file,sep="\t",header=F,check.names=F)
rt <- rt %>%
  mutate(V2 = recode(V2,
    "<0" = "<0",
    "0"  = "0",
    "1"  = "1",
    ._default = as.character(V2),  # 保留原始值
    ">100" = ifelse(V2 > 100, ">100", as.character(V2))
  ))

#x_axis_order <- c("<0", "0", "1", as.character(2:100), ">100")
x_axis_order <- c("<=0", "1", as.character(2:99), ">=100")
rt$V2 <- factor(rt$V2, levels = x_axis_order)

g <- ggplot(rt, aes(x=factor(V2),y = V1))+
  geom_col(fill='grey') + labs(x="Ratio of the distance from PA to transcriptional terminal to 3UTR length (%)",y = "Abundance")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=0.5))
ggsave(g,filename = out_file,width = 15, height = 10, dpi = 300)
