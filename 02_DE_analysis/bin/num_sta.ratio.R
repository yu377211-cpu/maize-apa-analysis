library(data.table)
library(ggplot2)

args=commandArgs(trailingOnly=T)
mydata <- read.table(args[1],header=F,sep='\t')
mydata$V3 = mydata$V2/sum(mydata$V2)

p2 = ggplot(data = mydata,aes(x=V1,y=V3))+
  geom_point()+
  geom_line()+
  xlab("sample_mun")+#横坐标名称
  ylab("PAC_count%")+#纵坐标名称
  theme_bw()+#去掉背景灰色
  labs(title = args[3])+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(0,300),breaks = seq(0,300,10))+
  scale_y_continuous(limits = c(0,0.15),breaks = seq(0,0.15,0.05))
  
ggsave(args[2],p2)
