# 公共入口可被项目 workflow 调用，也可直接作为 CLI 运行。
main <- function(cli_args = commandArgs(trailingOnly = TRUE)) {
# Slingshot 独立轨迹入口：使用主 Seurat RDS 的低维坐标和生物学注释拟合 lineage。
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
})

# 功能：读取 --name value 形式的命令行参数。
arg_value <- function(name, default = NULL) {
  args <- cli_args
  index <- tail(which(args == name), 1L)
  if (!length(index) || index == length(args)) default else args[[index + 1L]]
}

input_rds <- arg_value("--input")
output_dir <- arg_value("--output")
start_cluster <- arg_value("--start-cluster")
if (is.null(input_rds) || is.null(output_dir) || is.null(start_cluster)) {
  stop("用法：Rscript slingshot_workflow.R --input annotated.rds --output results/slingshot --start-cluster RootCellType")
}
object_dir <- arg_value("--object-dir", output_dir)
assay <- arg_value("--assay", "RNA")
reduction <- arg_value("--reduction", "umap")
cluster_col <- arg_value("--cluster-col", "cell_type")
seed <- as.integer(arg_value("--seed", "123"))

set.seed(seed)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

seu <- readRDS(input_rds)
DefaultAssay(seu) <- assay
counts <- GetAssayData(seu, assay = assay, layer = "counts")
metadata <- seu@meta.data[colnames(counts), , drop = FALSE]
embedding <- Embeddings(seu, reduction = reduction)[colnames(counts), 1:2, drop = FALSE]
metadata[[cluster_col]] <- factor(metadata[[cluster_col]])
if (!(start_cluster %in% levels(metadata[[cluster_col]]))) {
  stop("start_cluster not found in ", cluster_col, ": ", start_cluster, call. = FALSE)
}

sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = metadata
)
reducedDims(sce)[[toupper(reduction)]] <- embedding
sce <- slingshot(
  sce,
  clusterLabels = cluster_col,
  reducedDim = toupper(reduction),
  start.clus = start_cluster
)

pseudotime <- slingPseudotime(sce, na = FALSE)
weights <- slingCurveWeights(sce)
write.csv(
  data.frame(cell = colnames(sce), pseudotime, check.names = FALSE),
  file.path(output_dir, "slingshot_pseudotime.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(cell = colnames(sce), weights, check.names = FALSE),
  file.path(output_dir, "slingshot_lineage_weights.csv"),
  row.names = FALSE
)

cluster_factor <- colData(sce)[[cluster_col]]
NPG_COLORS <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"
)
n_levels <- nlevels(cluster_factor)
cluster_colors <- if (n_levels <= length(NPG_COLORS)) {
  NPG_COLORS[seq_len(n_levels)]
} else {
  extra_n <- n_levels - length(NPG_COLORS)
  c(NPG_COLORS, colorRampPalette(NPG_COLORS)(2L * extra_n + 1L)[seq.int(2L, 2L * extra_n, by = 2L)])
}
cluster_colors <- setNames(cluster_colors, levels(cluster_factor))
# 功能：绘制 Slingshot 降维坐标、细胞标签和 lineage 曲线（函数：plot_trajectory）。
plot_trajectory <- function() {
  plot(
    embedding,
    col = cluster_colors[cluster_factor],
    pch = 16,
    cex = 0.5,
    xlab = colnames(embedding)[1],
    ylab = colnames(embedding)[2]
  )
  lines(SlingshotDataSet(sce), lwd = 2, col = "black")
  legend("topright", legend = levels(cluster_factor), col = cluster_colors, pch = 16, cex = 0.7, bty = "n")
}
pdf(file.path(output_dir, "slingshot_trajectory.pdf"), width = 8, height = 7)
plot_trajectory()
dev.off()
png(file.path(output_dir, "slingshot_trajectory.png"), width = 8, height = 7, units = "in", res = 600)
plot_trajectory()
dev.off()

object_path <- file.path(object_dir, "slingshot_sce.rds")
saveRDS(sce, object_path)
writeLines(c(
  "# Slingshot 轨迹分析结果",
  "",
  paste0("- 输入对象：", input_rds),
  paste0("- 聚类列：", cluster_col),
  paste0("- 起始 cluster：", start_cluster),
  paste0("- lineage 数量：", ncol(pseudotime)),
  paste0("- Slingshot SCE：", object_path),
  "- 起始 cluster 必须由已知生物学起点指定；lineage 和 pseudotime 是低维结构的计算推断。"
), file.path(output_dir, "README_results_CN.md"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

}
if (sys.nframe() == 0L) main()
