args <- commandArgs(trailingOnly = TRUE)

fu_file <- args[1]
zheng_file <-args[2]
output_txt <- args[3]
output_png <- args[4]
output_pdf <- args[5]

source("/share2/pub/xingsl/xingsl/data_download/MZAPA/ZY/post_analysis/PAC_CS_RNAseq_feature/bin/ATGC_4.2.2.r")
ATGC_1_plot(x_fu_file = fu_file, x_zheng_file = zheng_file,zhongjian_file = output_txt,output_png=output_png,output_pdf=output_pdf)

