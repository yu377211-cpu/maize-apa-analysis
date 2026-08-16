library(argparser)
library(reshape2)
library(ggplot2)
library(gplots)  # 用于heatmap.2函数
library(RColorBrewer)  # 用于颜色调色板

argv <- arg_parser('')
argv <- add_argument(argv,"--fpkm", help="the normalized readcount or fpkm file")
argv <- add_argument(argv,"--condition", help="the condition file")
argv <- add_argument(argv,"--outdir", help="the directory of output file")
argv <- parse_args(argv)

fpkm <- argv$fpkm
condition <- argv$condition
outdir <- argv$outdir

setwd(outdir)
rawdata <- as.matrix(read.delim(fpkm, header=T, row.names=1,check.names=F,stringsAsFactors=FALSE))
class(rawdata)
groupdata <- read.delim(condition, header=T,check.names=F)

# 处理样本名称匹配
condition_samples <- as.character(groupdata[,'sample'])
fpkm_samples <- colnames(rawdata)

# 检查fpkm样本名是否包含.bed后缀
if(all(grepl("\\.bed$", fpkm_samples))) {
  # 创建匹配的样本名映射
  clean_samples <- sub("\\.bed$", "", fpkm_samples)
  matched_indices <- which(clean_samples %in% condition_samples)
  
  if(length(matched_indices) == 0) {
    stop("Error: No matching samples found between condition file and fpkm file")
  }
  data <- rawdata[, matched_indices, drop=FALSE]  # 只保留condition文件中指定的样本
  colnames(data) <- clean_samples[matched_indices]
}else{
  # 如果没有.bed后缀，直接匹配
  matched_indices <- which(fpkm_samples %in% condition_samples)
  
  if(length(matched_indices) == 0) {
    stop("Error: No matching samples found between condition file and fpkm file")
  }
  
  data <- rawdata[, matched_indices, drop=FALSE]
}

data <- log2(data+1)
correlation <- cor(data, method='pearson', use='pairwise.complete.obs')
rownames(correlation) <- colnames(correlation) <- colnames(data)  # 确保行列名一致且不带.bed
correlation <- round(correlation,3)
correlations <- melt(correlation)

if(ncol(data)<10){
        size=4
}else if(ncol(data)<15){
        size=3
}else if(ncol(data)<20){
        size=2
}else if(ncol(data)<30){
        size=1.5
}else{
        size=1
}

p <- ggplot(correlations,aes(Var1,Var2,label=value)) +
        geom_tile(aes(fill=value),colour="white") +
        scale_fill_gradient(name=expression(R^2),low="white",high="#4876FF") +
        theme(panel.background = element_rect(fill='white', colour='white')) +
        labs(x="",y="", title="Pearson correlation coefficient") +
        theme(legend.position="right",axis.text.x=element_text(angle=45,vjust=1,hjust=1),plot.title = element_text(hjust = 0.5)) +
        coord_fixed() +
        geom_text(size=size)

ggsave(filename='correlation.pdf',plot=p)
ggsave(filename='correlation.png',plot=p,type="cairo-png")
write.table(correlation,file='correlation.xls',sep='\t',quote=F)

# 创建带聚类图的热图
pdf("correlation_with_clustering.pdf", width=8, height=8)
heatmap.2(correlation,
          trace="none",
          col=brewer.pal(9, "Blues"),
          main="Pearson Correlation with Clustering",
          margins=c(10,10),
          cexRow=0.8,
          cexCol=0.8,
          key=TRUE,
          keysize=1.5,
          density.info="none")
dev.off()

png("correlation_with_clustering.png", width=800, height=800, res=150)
heatmap.2(correlation,
          trace="none",
          col=brewer.pal(9, "Blues"),
          main="Pearson Correlation with Clustering",
          margins=c(10,10),
          cexRow=0.8,
          cexCol=0.8,
          key=TRUE,
          keysize=1.5,
          density.info="none")
dev.off()
