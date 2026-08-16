library(data.table)
library(ggplot2)
library(ggpubr)

args <- commandArgs(trailingOnly = TRUE)
gene_file <- args[1]
pac_file <- args[2]
pos_file <- args[3]
output_file <- args[4]

gene <- fread(gene_file)
colnames(gene) <- c("gene","group")
pac <- fread(pac_file)

pos <- fread(pos_file)
pos  <- pos[,c("V5","V7")]
colnames(pos) <- c("PACid","gene")

sp <- merge(gene,pos)
sp <- merge(sp,pac,by="PACid",all.x=T)
sp_melt <- melt(sp[,-c(1,2)],id.vars = "group")
sp_melt[sp_melt$group>=7,"group"] <- 7
sp_melt <- sp_melt[value>0,]



p <- ggboxplot(sp_melt,x="group",y="value",color="group",outlier.shape = 1, bxp.errorbar = T, legend="none")+
   ylim(0,0.4) +  scale_x_discrete(labels=c("1","2","3","4","5","6",">=7"))+
  xlab("Number of PACs in gene") + ylab("Raletive abundance")
   
ggsave(output_file,p)
