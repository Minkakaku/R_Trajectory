# 公共入口可被项目 workflow 调用，也可直接作为 CLI 运行。
main <- function(cli_args = commandArgs(trailingOnly = TRUE)) {
#!/usr/bin/env Rscript

# tradeSeq 独立轨迹差异分析入口
# 可直接接收主 Seurat RDS 并在内部运行 Slingshot，也可接收已生成的 Slingshot SCE。

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
  library(tradeSeq)
  library(ggplot2)
})

# 功能：读取 --name value 形式的命令行参数。
arg_value <- function(name, default = NULL) {
  args <- cli_args
  index <- tail(which(args == name), 1L)
  if (!length(index) || index == length(args)) default else args[[index + 1L]]
}

# 功能：把 ggplot 同时保存为 PDF 和 600 dpi PNG。
save_plot <- function(plot, stem, width = 8, height = 6) {
  plot <- plot + ggplot2::labs(subtitle = NULL) + ggplot2::theme(plot.subtitle = ggplot2::element_blank())
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height)
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 600)
}

# 功能：从主 Seurat RDS 构建 Slingshot SingleCellExperiment。
build_slingshot_sce <- function(seu, assay, reduction, cluster_col, start_cluster) {
  seu <- SeuratObject::JoinLayers(seu, assay = assay)
  counts <- SeuratObject::LayerData(seu, assay = assay, layer = "counts")
  metadata <- seu[[]][colnames(counts), , drop = FALSE]
  embedding <- Seurat::Embeddings(seu, reduction = reduction)[colnames(counts), 1:2, drop = FALSE]
  metadata[[cluster_col]] <- factor(metadata[[cluster_col]])
  if (!start_cluster %in% levels(metadata[[cluster_col]])) {
    stop("start_cluster 不在 ", cluster_col, " 中：", start_cluster)
  }
  sce <- SingleCellExperiment(
    assays = list(counts = counts),
    colData = metadata
  )
  reducedDims(sce)[[toupper(reduction)]] <- embedding
  slingshot(
    sce,
    clusterLabels = cluster_col,
    reducedDim = toupper(reduction),
    start.clus = start_cluster
  )
}

# 功能：为 tradeSeq 四类检验添加 BH qvalue 并稳定排序。
add_qvalue <- function(table) {
  table$gene <- rownames(table)
  table$qvalue <- p.adjust(table$pvalue, method = "BH")
  table[order(table$qvalue, table$pvalue), , drop = FALSE]
}

# 功能：拟合 NB-GAM，并输出沿拟时序变化、起终点和 lineage 差异检验。
run_tradeseq <- function(
  input_rds,
  output_dir,
  input_type = "auto",
  assay = "RNA",
  reduction = "umap",
  cluster_col = "cell_type",
  start_cluster = NULL,
  nknots = 6L,
  workers = 4L,
  top_genes = 6L
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  input <- readRDS(input_rds)
  if (input_type == "auto") {
    input_type <- if (inherits(input, "Seurat")) "seurat" else "slingshot"
  }
  if (input_type == "seurat") {
    if (is.null(start_cluster)) stop("Seurat 输入必须通过 --start-cluster 指定生物学起点。")
    sce <- build_slingshot_sce(input, assay, reduction, cluster_col, start_cluster)
  } else {
    sce <- input
  }

  counts <- assay(sce, "counts")
  pseudotime <- slingPseudotime(sce, na = FALSE)
  cell_weights <- slingCurveWeights(sce)
  bpparam <- if (.Platform$OS.type == "windows") {
    BiocParallel::SnowParam(workers = workers)
  } else {
    BiocParallel::MulticoreParam(workers = workers)
  }

  set.seed(123)
  gam <- fitGAM(
    counts = counts,
    pseudotime = pseudotime,
    cellWeights = cell_weights,
    nknots = nknots,
    parallel = workers > 1,
    BPPARAM = bpparam,
    verbose = TRUE
  )

  association <- add_qvalue(associationTest(gam))
  start_end <- add_qvalue(startVsEndTest(gam))
  write.csv(association, file.path(output_dir, "association_test.csv"), row.names = FALSE)
  write.csv(start_end, file.path(output_dir, "start_vs_end_test.csv"), row.names = FALSE)

  if (ncol(cell_weights) > 1) {
    write.csv(add_qvalue(patternTest(gam)), file.path(output_dir, "pattern_test.csv"), row.names = FALSE)
    write.csv(add_qvalue(diffEndTest(gam)), file.path(output_dir, "differential_end_test.csv"), row.names = FALSE)
  }

  # 对 associationTest 最显著基因绘制平滑表达曲线，展示而非替代统计检验。
  selected <- head(association$gene[!is.na(association$qvalue)], top_genes)
  for (gene in selected) {
    plot <- plotSmoothers(gam, counts, gene = gene) + ggtitle(gene) + theme_classic()
    save_plot(plot, file.path(output_dir, paste0("smoother_", make.names(gene))))
  }

  saveRDS(sce, file.path(object_dir, "slingshot_sce.rds"))
  saveRDS(gam, file.path(object_dir, "tradeseq_gam.rds"))
  writeLines(c(
    "# tradeSeq 轨迹差异表达结果",
    "",
    paste0("- 输入对象：", normalizePath(input_rds, winslash = "/", mustWork = FALSE)),
    paste0("- 输入类型：", input_type, "；lineage 数量：", ncol(cell_weights), "。"),
    paste0("- nknots：", nknots, "。"),
    "- associationTest：任一 lineage 上随拟时序变化的基因。",
    "- startVsEndTest：各 lineage 起点与终点的表达差异。",
    "- patternTest/diffEndTest：多 lineage 时比较表达模式或终点。",
    "- 各检验分别使用 BH qvalue；pseudotime 是相对顺序，不等于真实时间。"
  ), file.path(output_dir, "README_results_CN.md"))
  capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
  invisible(list(sce = sce, gam = gam))
}

input_rds <- arg_value("--input")
output_dir <- arg_value("--output")
object_dir <- arg_value("--object-dir", file.path("data/processed", "tradeseq"))
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
if (is.null(input_rds) || is.null(output_dir)) {
  stop("用法：Rscript tradeseq_workflow.R --input annotated.rds --output results/tradeseq --input-type seurat --start-cluster RootCellType")
}
run_tradeseq(
  input_rds = input_rds,
  output_dir = output_dir,
  input_type = arg_value("--input-type", "auto"),
  assay = arg_value("--assay", "RNA"),
  reduction = arg_value("--reduction", "umap"),
  cluster_col = arg_value("--cluster-col", "cell_type"),
  start_cluster = arg_value("--start-cluster", NULL),
  nknots = as.integer(arg_value("--nknots", "6")),
  workers = as.integer(arg_value("--workers", "4")),
  top_genes = as.integer(arg_value("--top-genes", "6"))
)

}
if (sys.nframe() == 0L) main()
