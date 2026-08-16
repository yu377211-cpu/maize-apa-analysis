library(argparser)
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

# 检查列名
if(!all(c("sample", "groups") %in% colnames(groupdata))) {
  stop("Condition file must contain 'sample' and 'groups' columns")
}

# 样本匹配处理
condition_samples <- as.character(groupdata$sample)
fpkm_samples <- colnames(rawdata)

if(all(grepl("\\.bed$", fpkm_samples))) {
  clean_samples <- sub("\\.bed$", "", fpkm_samples)
  matched_indices <- which(clean_samples %in% condition_samples)
  data <- rawdata[, matched_indices, drop=FALSE]
  colnames(data) <- sub("\\.bed$", "", colnames(data))
} else {
  matched_indices <- which(fpkm_samples %in% condition_samples)
  data <- rawdata[, matched_indices, drop=FALSE]
}

if(length(matched_indices) == 0) {
  stop("No matching samples found between condition file and fpkm file")
}

# 数据转换
data <- log2(data + 1)
correlation <- cor(data, method='pearson')
rownames(correlation) <- colnames(correlation) <- colnames(data)
correlation <- round(correlation, 3)
write.table(correlation, file='correlation.xls', sep='\t', quote=F)

# 动态参数计算函数
calc_params <- function(n_samples, n_groups) {
  params <- list()
  
  # 基础尺寸计算
  base_size <- max(8, min(12, 12 - 0.1*(n_samples-40)))
  
  # 热图参数
  params$margins <- c(
    max(8, min(15, 15 - 0.2*(n_samples-40))),
    max(8, min(15, 15 - 0.2*(n_samples-40)))
  )
  params$cex <- max(0.4, min(1, 1 - 0.015*(n_samples-40)))
  params$notecex <- max(0.3, params$cex*0.7)
  
  # 图例参数
  params$legend_cex <- max(0.5, min(0.8, 0.8 - 0.01*(n_groups-20)))
  params$legend_cols <- ifelse(n_groups > 30, 2, 1)
  
  # 输出尺寸
  params$pdf_width <- max(10, min(20, 10 + 0.2*(n_samples-40)))
  params$pdf_height <- params$pdf_width
  params$png_size <- max(1500, min(3000, 1500 + 20*(n_samples-40)))
  params$png_res <- ifelse(n_samples > 100, 200, 150)
  
  return(params)
}

# 获取参数
n_samples <- ncol(data)
n_groups <- length(unique(groupdata$groups))
params <- calc_params(n_samples, n_groups)

# 优化分组颜色（支持更多分组）
if(n_groups <= 12) {
  group_colors <- brewer.pal(max(3, n_groups), "Set3")
} else {
  group_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_groups)
}
names(group_colors) <- unique(groupdata$groups)

# 生成热图（PDF）
pdf("optimized_heatmap.pdf", width=params$pdf_width, height=params$pdf_height)
cor_range <- range(correlation, na.rm = TRUE)
# 智能颜色设置
if(diff(cor_range) < 0.3) {
  # 窄范围使用较少颜色
  n_colors <- max(3, min(7, ceiling(diff(cor_range)/0.05)))
  color_palette <- colorRampPalette(brewer.pal(n_colors, "Blues"))(n_colors)
  color_breaks <- seq(cor_range[1], cor_range[2], length.out = n_colors + 1)
} else {
  # 宽范围使用较多颜色但不超过9种
  n_colors <- min(9, ceiling(diff(cor_range)/0.1))
  color_palette <- colorRampPalette(brewer.pal(n_colors, "Blues"))(n_colors)
  color_breaks <- seq(cor_range[1], cor_range[2], length.out = n_colors + 1)
}
heatmap.2(correlation,
          trace="none",
          col=brewer.pal(9, "Blues"),
          main=paste("Pearson Correlation (", n_samples, "samples)"),
          margins=params$margins,
          cexRow=params$cex,
          cexCol=params$cex,
          key=TRUE,
          keysize=1,
          density.info="none",
          cellnote=if(n_samples <= 50) correlation else NA,
          notecol=ifelse(correlation > 0.7, "white", "black"),
          notecex=params$notecex,
          RowSideColors=group_colors[groupdata$groups[match(colnames(data), groupdata$sample)]],
          hclustfun=function(x) hclust(x, method="average"),
          lwid=c(1,4),
          lhei=c(1,4),
          labRow=NA,   #去除右侧的样本id
          labCol=if(n_samples > 50) NA else colnames(data),
          Rowv=TRUE,
          Colv=TRUE,
          dendrogram="row",
          symkey=FALSE,  # Disable symmetric key scaling
          key.xlab="Correlation",
          key.ylab="Density",
          key.title=NA,
)

# 智能图例布局
legend_pos <- "right"
legend(legend_pos,
       legend=names(group_colors),
       col=group_colors,
       pch=15,
       pt.cex=1.5,
       cex=params$legend_cex,
       ncol=params$legend_cols,
       title="Sample Groups")

dev.off()

# 生成热图（PNG）
png("optimized_heatmap.png", width=params$png_size, height=params$png_size, res=params$png_res)
heatmap.2(correlation,
          trace="none",
          col=brewer.pal(9, "Blues"),
          main=paste("Pearson Correlation (", n_samples, "samples)"),
          margins=params$margins,
          cexRow=params$cex,
          cexCol=params$cex,
          key=TRUE,
          keysize=1,
          density.info="none",
          cellnote=if(n_samples <= 50) correlation else NA,
          notecol=ifelse(correlation > 0.7, "white", "black"),
          notecex=params$notecex,
          RowSideColors=group_colors[groupdata$groups[match(colnames(data), groupdata$sample)]],
          hclustfun=function(x) hclust(x, method="average"),
          lwid=c(1,4),
          lhei=c(1,4),
          labRow=NA,
          labCol=if(n_samples > 50) NA else colnames(data),
          Rowv=TRUE,
          Colv=TRUE,
          dendrogram="row",
          symkey=FALSE,  # Disable symmetric key scaling
          key.xlab="Correlation",
          key.ylab="Density",
          key.title=NA,
)
legend_pos <- "right"

legend(legend_pos,
       legend=names(group_colors),
       col=group_colors,
       pch=15,
       pt.cex=1.5,
       cex=params$legend_cex,
       ncol=params$legend_cols,
       title="Sample Groups")
dev.off()
