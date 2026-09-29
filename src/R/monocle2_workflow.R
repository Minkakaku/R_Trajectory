# 公共入口可被项目 workflow 调用，也可直接作为 CLI 运行。
main <- function(cli_args = commandArgs(trailingOnly = TRUE)) {
# rm(list = ls())
suppressPackageStartupMessages({
  library(monocle)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyverse)
  library(ggsci)
  library(patchwork)
  library(stringr)
  library(getopt)
})

# dplyr compatibility shims for packages that still call legacy underscore verbs
# Implemented with modern dplyr/rlang APIs so no environment downgrade is required.
# 功能：为 Monocle2 需要的旧版 dplyr 接口提供兼容实现（函数：patch_dplyr_legacy）。
patch_dplyr_legacy <- function() {
  ns <- asNamespace('dplyr')
  # 功能：向 dplyr 命名空间写入一个兼容函数（函数：patch_one）。
  patch_one <- function(fname, fdef) {
    if (exists(fname, envir = ns, inherits = FALSE) && bindingIsLocked(fname, ns)) {
      unlockBinding(fname, ns)
    }
    assign(fname, fdef, envir = ns)
    lockBinding(fname, ns)
  }

  patch_one('group_by_', function(.data, ..., .dots = list(), add = FALSE) {
    quos_from_dots <- rlang::enquos(...)
    quos_from_chr <- lapply(.dots, rlang::sym)
    dplyr::group_by(.data, !!!quos_from_dots, !!!quos_from_chr, .add = add)
  })

  patch_one('select_', function(.data, ..., .dots = list()) {
    quos_from_dots <- rlang::enquos(...)
    quos_from_chr <- lapply(.dots, rlang::sym)
    dplyr::select(.data, !!!quos_from_dots, !!!quos_from_chr)
  })

  patch_one('rename_', function(.data, ..., .dots = list()) {
    if (length(.dots) > 0) {
      rename_exprs <- lapply(names(.dots), function(nm) rlang::expr(!!rlang::sym(nm) := !!rlang::sym(.dots[[nm]])))
      dplyr::rename(.data, !!!rename_exprs)
    } else {
      dplyr::rename(.data, ...)
    }
  })

  patch_one('summarise_', function(.data, ..., .dots = list()) {
    dplyr::summarise(.data, !!!.dots, ...)
  })

  patch_one('mutate_', function(.data, ..., .dots = list()) {
    dplyr::mutate(.data, !!!.dots, ...)
  })

  patch_one('filter_', function(.data, ..., .dots = list()) {
    dplyr::filter(.data, !!!.dots, ...)
  })

  patch_one('arrange_', function(.data, ..., .dots = list()) {
    dplyr::arrange(.data, !!!lapply(.dots, rlang::sym), ...)
  })

  invisible(TRUE)
}

patch_dplyr_legacy()

# 功能：检查 Monocle2 输入是否为非负、近似整数的原始 counts（函数：validate_raw_counts）。
validate_raw_counts <- function(mat) {
  if (!inherits(mat, c('matrix', 'dgCMatrix', 'dgTMatrix', 'dgRMatrix'))) {
    warning('counts matrix class is unusual: ', paste(class(mat), collapse = ', '))
  }
  values <- if (inherits(mat, 'matrix')) as.vector(mat) else mat@x
  rng <- range(values)
  if (isTRUE(rng[1] < 0)) {
    stop(paste0(
      'RNA@counts contains negative values (min = ', signif(rng[1], 6), '). ',
      'Monocle2 with negbinomial.size() requires raw non-negative count matrix. ',
      'Please rebuild the Seurat/RDS object from raw counts or h5ad raw/layers["counts"].'
    ))
  }
  if (any(abs(values - round(values)) > 1e-8, na.rm = TRUE)) {
    warning('RNA@counts contains non-integer values; Monocle2 expects raw counts. Please verify the input object.')
  }
}

command = matrix(c(
  'rds','r',1,'character',
  'type','t',2,'character',
  'downcluster','D',0,'logical',
  'maxcells','m',2,'numeric',
  'cluster','c',2,'character',
  'ptsize','p',2,'numeric',
  'BEAM','B',0,'logical',
  'output','o',1,'character',
  'root_col','C',2,'character',
  'root_group','G',2,'character',
  'object_dir','O',2,'character'
), byrow=TRUE, ncol=4)
args = getopt(command, opt = cli_args)

# 功能：按 cluster 分层下采样，避免大群体完全主导轨迹（函数：sample_cells_per_cluster）。
sample_cells_per_cluster <- function(seurat_object, colnames='leiden', num_cells=1000) {
  Idents(seurat_object) = colnames
  clusters <- Idents(seurat_object)
  cluster_ids <- levels(clusters)
  sampled_cell_ids <- c()
  for (cluster_id in cluster_ids) {
    cell_ids_in_cluster <- names(clusters[clusters == cluster_id])
    num_cells_to_sample <- min(num_cells, length(cell_ids_in_cluster))
    sampled_ids <- sample(cell_ids_in_cluster, size = num_cells_to_sample, replace = FALSE)
    sampled_cell_ids <- c(sampled_cell_ids, sampled_ids)
  }
  subset(seurat_object, cells = sampled_cell_ids)
}

# 功能：保存 Monocle2 轨迹图及按标签分面的轨迹图（函数：monocleplot）。
monocleplot <- function(cds, colnames, ptsize=1) {
  p = plot_cell_trajectory(cds, color_by=colnames, size=ptsize, show_backbone=TRUE)
  ggsave(paste0('monocle.trajectory.', colnames, '.png'), p, width=7, height=7.5, dpi=600)
  ggsave(paste0('monocle.trajectory.', colnames, '.pdf'), p, width=7, height=7.5)

  ncols = ceiling(sqrt(length(unique(cds@phenoData@data[, colnames]))))
  nrows = ceiling(length(unique(cds@phenoData@data[, colnames])) / ncols)
  p = plot_cell_trajectory(cds, color_by=colnames, size=ptsize, show_backbone=TRUE) +
    facet_wrap(vars(!!sym(colnames)), nrow=nrows)
  ggsave(paste0('monocle.trajectory.', colnames, '.faceted.png'), p, width=4*ncols, height=4.2*nrows, dpi=600)
  ggsave(paste0('monocle.trajectory.', colnames, '.faceted.pdf'), p, width=4*ncols, height=4.2*nrows)
}

rdsfile = args$rds
if (!is.null(args$cluster)) { clustercolname = args$cluster } else { clustercolname = 'leiden' }
if (is.null(args$downcluster)) { args$downcluster = FALSE }
if (!is.null(args$downcluster)) { downcluster = args$downcluster } else { downcluster = TRUE }
if (is.null(args$BEAM)) { args$BEAM = FALSE }
if (!is.null(args$BEAM)) { BEAM_flag = args$BEAM } else { BEAM_flag = TRUE }
if (!is.null(args$ptsize)) { ptsize = args$ptsize } else { ptsize = 1 }
root_col = if (!is.null(args$root_col)) args$root_col else NA
root_group = if (!is.null(args$root_group)) args$root_group else NA

pbmc <- readRDS(rdsfile)
if (!is.null(args$type)) typeuse = str_split(args$type, pattern=',')[[1]]

project_dir <- getwd()
outpath <- args$output
object_dir <- if (!is.null(args$object_dir)) args$object_dir else file.path(project_dir, 'data', 'processed')
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
outpath <- normalizePath(outpath, winslash = '/', mustWork = FALSE)
object_dir <- normalizePath(object_dir, winslash = '/', mustWork = FALSE)
previous_dir <- setwd(outpath)
on.exit(setwd(previous_dir), add = TRUE)

if (downcluster) {
  if (!is.null(args$maxcells)) { maxcell = args$maxcells } else { maxcell = 1000 }
  set.seed(123)
  pbmc = sample_cells_per_cluster(pbmc, num_cells=maxcell, colnames=clustercolname)
  print('Cell number of each cluster:')
  print(table(pbmc@meta.data[, clustercolname]))
}

pbmc <- SeuratObject::JoinLayers(pbmc, assay = 'RNA')
counts_mat <- SeuratObject::LayerData(pbmc, assay = 'RNA', layer = 'counts')
p_data <- pbmc@meta.data[colnames(counts_mat), , drop = FALSE]
print('Meta data columns:')
print(colnames(p_data))

f_data <- data.frame(gene_short_name=rownames(counts_mat), row.names=rownames(counts_mat))
pd <- new('AnnotatedDataFrame', data = p_data)
fd <- new('AnnotatedDataFrame', data = f_data)
validate_raw_counts(counts_mat)
cds <- newCellDataSet(counts_mat,
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 0.5,
  expressionFamily = negbinomial.size()
)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr=0.1)
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed > 10))

print('Analysis dpFeature:')
diff <- differentialGeneTest(cds[expressed_genes, ], fullModelFormulaStr = paste0('~', clustercolname))
deg <- subset(diff, qval < 0.01)
deg <- deg %>% arrange(qval, desc(num_cells_expressed))
if (nrow(deg) == 0) stop('No ordering genes passed qval < 0.01.')
deg_output <- data.frame(Genes = rownames(deg), deg[, c('pval','qval','num_cells_expressed')])
write.table(deg_output, file='monocle.DEG.xls', col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)
ordergene <- rownames(deg)[1:min(1000, nrow(deg))]
cds <- setOrderingFilter(cds, ordergene)

pdf('monocle.ordergenes.pdf'); plot_ordering_genes(cds); dev.off()
png('monocle.ordergenes.png', width=7, height=7.5, units='in', res=600); plot_ordering_genes(cds); dev.off()

print('ReduceDimension and orderCells:')
cds <- reduceDimension(cds, max_components = 2, reduction_method='DDRTree')
cds <- orderCells(cds)

chosen_root_state <- NA
if (!is.na(root_col) && !is.na(root_group)) {
  if (!(root_col %in% colnames(pData(cds)))) stop(paste0('root_col not found: ', root_col))
  root_idx <- pData(cds)[[root_col]] %in% root_group
  if (sum(root_idx) == 0) stop(paste0('No cells matched root_group in ', root_col, ': ', root_group))
  state_tab <- sort(table(pData(cds)[root_idx, 'State']), decreasing = TRUE)
  chosen_root_state <- names(state_tab)[1]
  print(paste0('Chosen root_state = ', chosen_root_state, ' based on ', root_col, ' == ', root_group))
  print(state_tab)
  cds <- orderCells(cds, root_state = chosen_root_state)
}

object_path <- file.path(object_dir, 'monocle2_cds.rds')
saveRDS(cds, object_path)
writeLines(as.character(chosen_root_state), 'root_state.txt')

plotcol = c('sample','group','State', clustercolname)
if (!is.null(args$type)) plotcol = c(plotcol, typeuse)
if ('cell_type' %in% colnames(pbmc@meta.data)) plotcol = c(plotcol, 'cell_type')
plotcol = unique(plotcol)
for (i in plotcol) {
  if (i %in% colnames(cds@phenoData@data)) monocleplot(cds, colnames=i, ptsize=ptsize)
}

p = plot_cell_trajectory(cds, color_by='Pseudotime', size=ptsize, show_backbone=TRUE)
ggsave('monocle.trajectory.pseudotime.pdf', p, width=7, height=7.5)
ggsave('monocle.trajectory.pseudotime.png', p, width=7, height=7.5, dpi=600)
p = plot_complex_cell_trajectory(cds, x=1, y=2, color_by=clustercolname)
ggsave('monocle.trajectory.tree.pdf', p, width=7, height=7.5)
ggsave('monocle.trajectory.tree.png', p, width=7, height=7.5, dpi=600)

ncols = ceiling(sqrt(length(unique(cds@phenoData@data[, clustercolname]))))
nrows = ceiling(length(unique(cds@phenoData@data[, clustercolname])) / ncols)
p = plot_complex_cell_trajectory(cds, x=1, y=2, color_by=clustercolname) +
  facet_wrap(vars(!!sym(clustercolname)), nrow=nrows)
ggsave('monocle.trajectory.tree.faceted.pdf', p, width=4*ncols, height=4.2*nrows)
ggsave('monocle.trajectory.tree.faceted.png', p, width=4*ncols, height=4.2*nrows, dpi=600)

meta_cols = c('sample', clustercolname, 'Size_Factor', 'num_genes_expressed', 'Pseudotime', 'State')
meta_cols = meta_cols[meta_cols %in% colnames(cds@phenoData@data)]
pseudotime_df = cds@phenoData@data[, meta_cols, drop = FALSE]
df = data.frame(Cells=rownames(pseudotime_df), pseudotime_df)
write.table(df, 'monocle.metadata.xls', sep='\t', quote=FALSE, row.names=FALSE)

if (BEAM_flag) {
  n_branch = length(cds@auxOrderingData[['DDRTree']]$branch_points)
  for (branch_point in seq_len(n_branch)) {
    outbeam = paste0(outpath, '/BEAM_analysis.Node', branch_point)
    if (!dir.exists(outbeam)) dir.create(outbeam)
    print(paste0('Analysis Beam Node ', branch_point))
    BEAM_res = BEAM(cds, branch_point = branch_point, cores = 4, progenitor_method = 'duplicate')
    BEAM_res = BEAM_res[, c('gene_short_name','pval','qval')]
    BEAM_res = BEAM_res[order(BEAM_res$qval), ]
    BEAM_output = data.frame(gene = rownames(BEAM_res), BEAM_res, check.names = FALSE)
    write.csv(BEAM_output, paste0(outbeam, '/BEAM_res.csv'), row.names = FALSE)
    selectgenes = row.names(subset(BEAM_res, qval < 0.05))
    if (length(selectgenes) > 1) {
      p = plot_genes_branched_heatmap(cds[selectgenes[1:min(100, length(selectgenes))], ],
        branch_point = branch_point, num_clusters = 4, cores = 1,
        use_gene_short_name = TRUE, show_rownames = TRUE, return_heatmap = TRUE)
      ggsave(paste0(outbeam, '/BEAM_heatmap.png'), p$ph_res, width=7, height=8, dpi=600)
      ggsave(paste0(outbeam, '/BEAM_heatmap.pdf'), p$ph_res, width=7, height=8)
      p2 = plot_genes_branched_pseudotime(cds[selectgenes[1:min(20, length(selectgenes))], ],
        branch_point = branch_point, color_by = 'State', ncol = 4)
      ggsave(paste0(outbeam, '/BEAM_pseudotime_top20.png'), p2, width=15, height=12, dpi=600)
      ggsave(paste0(outbeam, '/BEAM_pseudotime_top20.pdf'), p2, width=15, height=12)
    }
  }
}

result_note <- c(
  '# Monocle2 拟时序分析结果',
  '',
  paste0('- 输入对象：', normalizePath(rdsfile, winslash = '/', mustWork = FALSE)),
  paste0('- 聚类列：', clustercolname),
  paste0('- root 条件：', ifelse(is.na(root_col), '未指定', paste0(root_col, ' = ', root_group))),
  paste0('- root state：', chosen_root_state),
  paste0('- 结果目录：', outpath),
  paste0('- Monocle2 CDS：', object_path),
  '- ordering genes 使用 dpFeature 差异检验筛选。',
  '- BEAM 仅在 --BEAM 开启时运行。'
)
writeLines(result_note, 'README_results_CN.md')
capture.output(sessionInfo(), file = 'sessionInfo.txt')
message('Monocle2 results: ', outpath)
message('Monocle2 object: ', object_path)

}
if (sys.nframe() == 0L) main()
