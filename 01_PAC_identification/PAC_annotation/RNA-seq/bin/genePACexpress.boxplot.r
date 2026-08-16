library(data.table)
library(ggplot2)
library(ggpubr)

args <- commandArgs(trailingOnly = TRUE)
gene_file <- args[1] #ERP011069.PAC.CS.3UTR.geneSta.condition
pac_file <- args[2] #ERP011069.3UTR.PAC.hq.TPM
pos_file <- args[3] #ERP011069.PAC.CS.3UTR.pos_strand.bed
output_file <- args[4] #ERP011069.PAC.CS.3UTR.TPM.box.png

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

 
p <- ggboxplot(sp_melt,x="group",y="value",color="group",legend="none", bxp.errorbar = T)+
   ylim(0,500) +  scale_x_discrete(labels=c("1","2","3","4","5","6",">=7"))+
  xlab("Number of PACs in gene") + ylab("PAC TPM")
   
ggsave(output_file, p)
