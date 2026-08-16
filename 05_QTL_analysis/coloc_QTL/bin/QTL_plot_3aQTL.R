library(optparse)

option_list <- list(
    make_option(c("-s","--snp"), type="character", help="SNP ID"),
    make_option(c("-g","--transcript"), type="character", help="Transcript ID"),
    make_option(c("-G","--genotype"), type="character",
                default="./Matrix_eQTL/Genotype_matrix.txt"),
    make_option(c("-P","--phenotype"), type="character",
                default="./Matrix_eQTL/Phenotype_matrix.txt"),
    make_option(c("-o","--out"), type="character",
                default="3aQTL")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Load matrices
gt <- read.table(opt$genotype, header=TRUE, sep="\t", check.names=FALSE)
pt <- read.table(opt$phenotype, header=TRUE, sep="\t", check.names=FALSE)

rownames(gt) <- gt[,1]
rownames(pt) <- pt[,1]
gt <- gt[,-1]
pt <- pt[,-1]

snp  <- opt$snp
gene <- opt$transcript

e1 <- as.numeric(pt[gene, ])
s1 <- as.numeric(gt[snp, ])

pdf(paste(opt$out, snp, gene, "pdf", sep="."), width=5, height=5)

boxplot(
    e1 ~ s1,
    lwd = 2,
    xaxt = "n",
    xlab = "Genotype",
    ylab = "Normalized PDUI",
    #ylim = c(0,1),
    main = paste(snp, gene, sep="  |  "),
    outline = FALSE
)

axis(1, at = 1:3, labels = c("REF", "HET", "ALT"))

stripchart(
    e1 ~ s1,
    vertical = TRUE,
    method = "jitter",
    add = TRUE,
    pch = 16,
    cex = 0.7,
    #col = rgb(0,0,0,0.5)
    col = c(rgb(102,194,165,max=255),rgb(252,141,98,max=255),rgb(141,160,203,max=255))
)

dev.off()

