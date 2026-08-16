library(argparser)
library(reshape2)
library(gplots)
library(RColorBrewer)

argv <- arg_parser('')
argv <- add_argument(argv,"--fpkm", help="the normalized readcount or fpkm file")
argv <- add_argument(argv,"--condition", help="the condition file with sample and groups columns")
argv <- add_argument(argv,"--outdir", help="the directory of output file")
argv <- parse_args(argv)

fpkm <- argv$fpkm
condition <- argv$condition
outdir <- argv$outdir

setwd(outdir)

# 读取数据
rawdata <- as.matrix(read.delim(fpkm, header=T, row.names=1, check.names=F, stringsAsFactors=FALSE))
groupdata <- read.delim(condition, header=T, check.names=F)

# 处理样本名称匹配
condition_samples <- as.character(groupdata[,'sample'])
fpkm_samples <- colnames(rawdata)

# 检查fpkm样本名是否包含.bed后缀
if(all(grepl("\\.bed$", fpkm_samples))) {
  clean_samples <- sub("\\.bed$", "", fpkm_samples)
  matched_indices <- which(clean_samples %in% condition_samples)
  # 筛选数据并处理样本名
  data <- rawdata[, matched_indices, drop=FALSE]
  colnames(data) <- sub("\\.bed$", "", colnames(data))
} else {
  matched_indices <- which(fpkm_samples %in% condition_samples)
  data <- rawdata[, matched_indices, drop=FALSE]
}

if(length(matched_indices) == 0) {
  stop("Error: No matching samples found between condition file and fpkm file")
}

# 对数转换
data <- log2(data + 1)

# 计算相关系数
correlation <- cor(data, method='pearson', use='pairwise.complete.obs')
rownames(correlation) <- colnames(correlation) <- colnames(data)
correlation <- round(correlation, 3)

# 保存相关系数矩阵
write.table(correlation, file='correlation.xls', sep='\t', quote=F)

# 准备分组信息
group_info <- groupdata[match(colnames(data), groupdata$sample), "groups", drop=FALSE]
group_levels <- unique(group_info$groups)
group_colors <- brewer.pal(length(group_levels), "Set1")[1:length(group_levels)]
names(group_colors) <- group_levels

# 创建带分组信息的综合热图
pdf("combined_correlation_heatmap.pdf", width=12, height=10)
heatmap.2(correlation,
          trace="none",
          col=brewer.pal(9, "Blues"),
          main="Pearson Correlation with Clustering and Group Info",
          margins=c(12,12),
          cexRow=1,
          cexCol=1,
          key=TRUE,
          keysize=1.2,
          density.info="none",
          cellnote=correlation,  # 显示相关系数值
          notecol=ifelse(correlation > 0.7, "white", "black"),      # 数值颜色
          notecex=0.7,          # 数值大小
          RowSideColors=group_colors[as.character(group_info$groups)],  # 行分组颜色
          ColSideColors=group_colors[as.character(group_info$groups)],  # 列分组颜色
          hclustfun=function(x) hclust(x, method="complete"),  # 聚类方法
          lwid=c(1.5,4),
          lhei=c(1.5,4))

# 添加分组图例
legend("topright", 
       legend=names(group_colors),
       col=group_colors,
       pch=15,
       pt.cex=1.5,
       cex=0.8,
       title="Sample Groups")
dev.off()

png("combined_correlation_heatmap.png", width=1200, height=1000, res=150)
heatmap.2(correlation,
          trace="none",
          col=brewer.pal(9, "Blues"),
          main="Pearson Correlation with Clustering and Group Info",
          margins=c(12,12),
          cexRow=0.8,
          cexCol=0.8,
          key=TRUE,
          keysize=1.2,
          density.info="none",
          cellnote=correlation,
          notecol=ifelse(correlation > 0.7, "white", "black"),
          notecex=0.7,
          RowSideColors=group_colors[as.character(group_info$groups)],
          ColSideColors=group_colors[as.character(group_info$groups)],
          hclustfun=function(x) hclust(x, method="complete"),
          lwid=c(1.5,4),
          lhei=c(1.5,4))

legend("topright", 
       legend=names(group_colors),
       col=group_colors,
       pch=15,
       pt.cex=1.5,
       cex=0.8,
       title="Sample Groups")
dev.off()
