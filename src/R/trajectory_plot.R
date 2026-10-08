# 从已拟合轨迹绘图；只复用坐标和路径，不重新拟合或修改对象。
trajectory_npg <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
                    "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85")

# 复用原生 Slingshot 路径或 Monocle/PAGA 线段，可只显示指定 lineage。
add_trajectory_graph <- function(plot, result, lineage = NULL) {
  edges <- result$segments
  if (identical(result$method, "Slingshot")) {
    stopifnot(all(c("x", "y", "lineage") %in% names(edges)))
    if (!is.null(lineage)) edges <- edges[edges$lineage == paste0("Lineage", lineage), , drop = FALSE]
    plot + ggplot2::geom_path(data = edges, ggplot2::aes(x = x, y = y, group = lineage),
      inherit.aes = FALSE, colour = "black", linewidth = 0.5)
  } else {
    stopifnot(all(c("x", "y", "xend", "yend") %in% names(edges)))
    plot + ggplot2::geom_segment(data = edges, ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE, colour = "black", linewidth = 0.45)
  }
}

# 分类或 pseudotime 使用 NPG；缺失 pseudotime 保留为灰色点，不丢弃细胞。
plot_existing_trajectory <- function(result, color, palette = NULL, point_size = NULL,
                                     lineage = NULL, title = result$method) {
  cells <- result$cells
  stopifnot(all(c("x", "y", color) %in% names(cells)), nrow(cells) > 0,
            all(is.finite(cells$x)), all(is.finite(cells$y)))
  value <- cells[[color]]
  if (is.numeric(value)) {
    cells[[color]][!is.finite(value)] <- NA_real_
    cells <- cells[order(cells[[color]], na.last = FALSE), , drop = FALSE]
  }
  if (is.null(point_size)) point_size <- Seurat::AutoPointSize(data = cells, raster = FALSE)
  plot <- ggplot2::ggplot(cells, ggplot2::aes(x = x, y = y, colour = .data[[color]])) +
    ggplot2::geom_point(size = point_size) + ggplot2::coord_fixed() +
    ggplot2::labs(title = title, subtitle = NULL, x = "Trajectory 1", y = "Trajectory 2", colour = color) +
    ggplot2::theme_classic()
  if (is.numeric(value)) {
    plot <- plot + ggplot2::scale_colour_gradientn(
      colours = if (is.null(palette)) c("#FFF5F0", "#F39B7F", "#E64B35", "#7A0019") else palette,
      na.value = "#C8C8C8", guide = ggplot2::guide_colorbar(display = "rectangles"))
  } else {
    labels <- if (is.factor(value)) levels(droplevels(value)) else sort(unique(as.character(value)))
    if (is.null(palette)) {
      colors <- trajectory_npg
      if (length(labels) > length(colors)) colors <- c(colors,
        grDevices::colorRampPalette(trajectory_npg)(length(labels) - length(colors) + 2)[-c(1, length(labels) - length(colors) + 2)])
      palette <- stats::setNames(colors[seq_along(labels)], labels)
    }
    plot <- plot + ggplot2::scale_colour_manual(values = palette, na.value = "#C8C8C8")
  }
  add_trajectory_graph(plot, result, lineage)
}
