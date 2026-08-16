suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--ae", type="character", help="3'aeQTL susieR result (PAC-level)"),
  make_option("--a",  type="character", help="3'aQTL susieR result (transcript-level)"),
  make_option("--map",type="character", help="PAC-transcript mapping"),
  make_option("--out",type="character", help="output prefix")
)

opt <- parse_args(OptionParser(option_list=option_list))

ae  <- fread(opt$ae)
a   <- fread(opt$a)
map <- fread(opt$map, col.names=c("PAC","transcript"))

# CS 过滤
filter_cs <- function(df,
                      purity_cut=0.8, max_cs_size=10000) {
  df[!is.na(cs) &
     cs_purity >= purity_cut &
     cs_size <= max_cs_size]
}

ae_f <- filter_cs(ae)
a_f  <- filter_cs(a)

# 按 locus 和 SNP 汇总 cs
collapse_cs <- function(df) {
  df[, .(variants = list(unique(variant_id))),
     by = .(locus_id, cs)]
}

ae_cs <- collapse_cs(ae_f)   # locus_id = PAC
a_cs  <- collapse_cs(a_f)    # locus_id = transcript

# 把 transcript cs 映射到PAC
a_cs_map <- merge(
  a_cs,
  map,
  by.x = "locus_id",
  by.y = "transcript",
  allow.cartesian = TRUE
)

setnames(a_cs_map, "locus_id", "transcript")

# 共定位判定函数
has_overlap <- function(v1, v2) {
  length(intersect(v1, v2)) > 0
}

# PAC 共定位扫描
hits <- list()

for (pac in intersect(unique(ae_cs$locus_id),
                       unique(a_cs_map$PAC))) {

  ae_sub <- ae_cs[locus_id == pac]
  a_sub  <- a_cs_map[PAC == pac]

  for (i in 1:nrow(ae_sub)) {
    for (j in 1:nrow(a_sub)) {

      if (has_overlap(ae_sub$variants[[i]],
                      a_sub$variants[[j]])) {

        hits[[length(hits)+1]] <- data.table(
          PAC = pac,
          ae_cs = ae_sub$cs[i],
          transcript = a_sub$transcript[j],
          a_cs = a_sub$cs[j],
          shared_variants = paste(
            intersect(ae_sub$variants[[i]],
                      a_sub$variants[[j]]),
            collapse = ","
          )
        )
      }
    }
  }
}

# 输出结果
detail <- if (length(hits) > 0) rbindlist(hits) else data.table()
fwrite(detail, paste0(opt$out, ".detail.txt"), sep="\t")

summary <- unique(detail[, .(PAC)])
summary[, colocalized := TRUE]

fwrite(summary, paste0(opt$out, ".summary.txt"), sep="\t")

#final <- merge(
#  pac_coloc,
#  diff_APA_PAC,
#  by="PAC",
#  all.x=TRUE
#)

