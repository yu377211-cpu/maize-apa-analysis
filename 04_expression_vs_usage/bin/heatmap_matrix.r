library(ComplexHeatmap)
library(circlize)

mat_main <- as.matrix(read.csv(
  "main_matrix_by_subpop.csv",
  row.names = 1,
  check.names = FALSE
))

# 颜色：蓝 = short-up / long-down，红 = long-up / short-down
col_fun <- colorRamp2(
  c(-1, 0, 1),
  c("#4575b4", "white", "#d73027")
)

ht_main <- Heatmap(
  mat_main,
  name = "APA–Expr\nCoupling",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_side = "left",
  column_title = "APA usage–expression coupling across subpopulations",
  heatmap_legend_param = list(
    at = c(-1, 0, 1),
    labels = c("Short–up / Long–down", "No coupling", "Long–up / Short–down")
  )
)

pdf("subpop_main_APA_expr_heatmap.pdf", width = 10, height = 4)
draw(ht_main)
dev.off()

mat_shared <- as.matrix(read.csv(
  "supp_matrix_by_subpop_shared.csv",
  row.names = 1,
  check.names = FALSE
))

ht_shared <- Heatmap(
  mat_shared,
  name = "Shared\nCoupling",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_title = "Shared APA usage–expression coupling events"
)

pdf("subpop_shared_APA_expr_heatmap.pdf", width = 10, height = 4)
draw(ht_shared)
dev.off()

