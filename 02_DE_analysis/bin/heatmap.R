library(pheatmap)
library(data.table)
args=commandArgs(trailingOnly=T)

A <- as.matrix(fread(args[1],header=T),rownames=1)
groupcol = read.table(args[2],sep="\t",header=T,row.names=1,check.names=F)

pheatmap(A,cluster_cols = T, #根据聚类需求，可修改
         show_colnames = T,border_color = NA,scale = "row",show_rownames =F,
         annotation_col = groupcol,filename = args[3])
