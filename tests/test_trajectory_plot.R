# 复用已有轨迹时保持全部细胞、坐标和路径；检查真实 ggplot 图层。
source("src/R/trajectory_plot.R")
result <- list(method = "Monocle2", cells = data.frame(x = 1:3, y = c(1, 3, 2),
  state = factor(c("a", "b", "a")), pseudotime = c(0, NA, 1)),
  segments = data.frame(x = 1, y = 1, xend = 3, yend = 2))
original <- result
categorical <- ggplot2::ggplot_build(plot_existing_trajectory(result, "state", point_size = 0.2))$data
continuous <- ggplot2::ggplot_build(plot_existing_trajectory(result, "pseudotime", point_size = 0.2))$data
stopifnot(identical(result, original), nrow(categorical[[1]]) == 3, nrow(continuous[[1]]) == 3,
          identical(categorical[[1]]$x, as.numeric(result$cells$x)),
          identical(categorical[[1]]$colour, trajectory_npg[c(1, 2, 1)]),
          "#C8C8C8" %in% continuous[[1]]$colour,
          identical(sort(continuous[[1]]$x), as.numeric(result$cells$x)))
result$method <- "Slingshot"
result$segments <- data.frame(x = 1:4, y = 1:4, lineage = rep(c("Lineage1", "Lineage2"), each = 2))
lineage <- ggplot2::ggplot_build(plot_existing_trajectory(result, "state", lineage = 1, point_size = 0.2))$data[[2]]
stopifnot(nrow(lineage) == 2)
cat("EXISTING_TRAJECTORY_PLOT_OK\n")
