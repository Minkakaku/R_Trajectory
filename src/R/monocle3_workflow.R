# 公共入口可被项目 workflow 调用，也可直接作为 CLI 运行。
main <- function(cli_args = commandArgs(trailingOnly = TRUE)) {
# Monocle3 独立轨迹入口：从主 Seurat RDS 学习 principal graph 并计算 pseudotime。
suppressPackageStartupMessages({
  library(monocle3)
  library(Seurat)
  library(SingleCellExperiment)
  library(ggplot2)
})

# 功能：读取 --name value 形式的命令行参数。
arg_value <- function(name, default = NULL) {
  args <- cli_args
  index <- tail(which(args == name), 1L)
  if (!length(index) || index == length(args)) default else args[[index + 1L]]
}

# 功能：解析逗号分隔的 metadata 列或 root 标签。
split_values <- function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1]])

input_rds <- arg_value("--input")
output_dir <- arg_value("--output")
if (is.null(input_rds) || is.null(output_dir)) {
  stop("用法：Rscript monocle3_workflow.R --input annotated.rds --output results/monocle3 --root-col cell_type --root-groups RootCellType")
}
object_dir <- arg_value("--object-dir", output_dir)
assay <- arg_value("--assay", "RNA")
reduction <- arg_value("--reduction", "umap")
cluster_col <- arg_value("--cluster-col", "cell_type")
group_cols <- split_values(arg_value("--group-cols", "cell_type,group,sample"))
root_col <- arg_value("--root-col", "cell_type")
root_groups <- split_values(arg_value("--root-groups", "RootCellType"))
use_partition <- as.logical(arg_value("--use-partition", "FALSE"))
close_loop <- as.logical(arg_value("--close-loop", "FALSE"))
workers <- as.integer(arg_value("--workers", "4"))
seed <- as.integer(arg_value("--seed", "123"))

# 功能：绘制并整理 Monocle3 principal graph 分析结果（函数：save_plot）。
save_plot <- function(plot, stem, width = 8, height = 7) {
  plot <- plot + ggplot2::labs(subtitle = NULL) + ggplot2::theme(plot.subtitle = ggplot2::element_blank())
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height)
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 600)
}

# 功能：根据指定起始细胞在 Monocle3 主图上选择 root 节点（函数：get_root_principal_node）。
get_root_principal_node <- function(cds, column, groups) {
  if (!(column %in% colnames(colData(cds)))) stop("root_col not found: ", column, call. = FALSE)
  root_cells <- colnames(cds)[colData(cds)[[column]] %in% groups]
  if (length(root_cells) == 0) stop("No cells matched root_groups in ", column, call. = FALSE)

  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), , drop = FALSE])
  root_vertex <- names(which.max(table(closest_vertex[root_cells, 1])))
  vertex_names <- igraph::V(principal_graph(cds)[["UMAP"]])$name
  if (root_vertex %in% vertex_names) return(root_vertex)
  vertex_names[as.integer(root_vertex)]
}

set.seed(seed)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

seu <- readRDS(input_rds)
seu <- SeuratObject::JoinLayers(seu, assay = assay)
DefaultAssay(seu) <- assay
counts <- SeuratObject::LayerData(seu, assay = assay, layer = "counts")
cell_metadata <- seu@meta.data[colnames(counts), , drop = FALSE]
gene_metadata <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))

cds <- new_cell_data_set(
  counts,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)
cds <- estimate_size_factors(cds)

if (reduction %in% Reductions(seu)) {
  embedding <- Embeddings(seu, reduction = reduction)[colnames(cds), 1:2, drop = FALSE]
  colnames(embedding) <- c("UMAP_1", "UMAP_2")
  reducedDims(cds)[["UMAP"]] <- embedding
  cds <- cluster_cells(cds, reduction_method = "UMAP")
} else {
  cds <- preprocess_cds(cds, num_dim = 50)
  cds <- reduce_dimension(cds, reduction_method = "UMAP")
  cds <- cluster_cells(cds, reduction_method = "UMAP")
}

cds <- learn_graph(
  cds,
  use_partition = use_partition,
  close_loop = close_loop
)
root_pr_node <- get_root_principal_node(cds, root_col, root_groups)
cds <- order_cells(cds, root_pr_nodes = root_pr_node)

plot_columns <- unique(c(cluster_col, group_cols))
plot_columns <- plot_columns[plot_columns %in% colnames(colData(cds))]
for (column in plot_columns) {
  p <- plot_cells(
    cds,
    color_cells_by = column,
    label_cell_groups = FALSE,
    label_leaves = TRUE,
    label_branch_points = TRUE
  ) + coord_fixed()
  save_plot(p, file.path(output_dir, paste0("trajectory_", column)))
}

p <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE
) + coord_fixed()
save_plot(p, file.path(output_dir, "trajectory_pseudotime"))

pseudotime_table <- data.frame(
  cell = colnames(cds),
  pseudotime = pseudotime(cds),
  as.data.frame(colData(cds))[colnames(cds), plot_columns, drop = FALSE],
  check.names = FALSE
)
write.csv(pseudotime_table, file.path(output_dir, "pseudotime_metadata.csv"), row.names = FALSE)

graph_test_result <- graph_test(cds, neighbor_graph = "principal_graph", cores = workers)
graph_test_result$gene <- rownames(graph_test_result)
graph_test_result <- graph_test_result[order(graph_test_result$q_value), , drop = FALSE]
write.csv(graph_test_result, file.path(output_dir, "trajectory_genes_graph_test.csv"), row.names = FALSE)

object_path <- file.path(object_dir, "monocle3_cds.rds")
saveRDS(cds, object_path)
writeLines(c(
  "# Monocle3 轨迹分析结果",
  "",
  paste0("- 输入对象：", input_rds),
  paste0("- 聚类列：", cluster_col),
  paste0("- root 条件：", root_col, " = ", paste(root_groups, collapse = ", ")),
  paste0("- root principal node：", root_pr_node),
  paste0("- use_partition：", use_partition, "；close_loop：", close_loop),
  paste0("- 结果目录：", output_dir),
  paste0("- Monocle3 CDS：", object_path),
  "- root 必须由已知生物学起点指定；pseudotime 表示图上的相对顺序，不等于真实时间。"
), file.path(output_dir, "README_results_CN.md"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

}
if (sys.nframe() == 0L) main()
