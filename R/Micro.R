#' @import pbapply
#' @import kohonen
#' @import reshape2
#' @import circlize
#' @import ggplot2
#' @import patchwork
#' @import ggsci
#' @import ggpubr
#' @import permute
#' @import lattice
#' @importFrom flowCore read.flowSet sampleNames exprs flowFrame markernames polygonGate write.FCS
#' @importFrom vegan vegdist adonis2
#' @import dplyr
# @importFrom dplyr select mutate mutate_all arrange summarise mutate_at %>%
#' @import rstatix
#' @import forcats
#' @importFrom cowplot plot_grid
#' @importFrom caret rfe rfeControl rfFuncs trainControl twoClassSummary train
#' @importFrom ggh4x force_panelsizes
#' @import pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom linkET mantel_test qcorrplot geom_square geom_mark geom_couple correlate nice_curvature color_pal
#' @importFrom pROC roc ci.auc ggroc ci.se
#' @import tibble
#' @import viridis
#' @importFrom scales trans_format
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom stats aov as.dist cmdscale p.adjust reorder setNames
#' @importFrom utils combn setTxtProgressBar txtProgressBar write.csv write.table
#' @importFrom randomForest randomForest


#' @title Train a SOM and export sample cluster abundances
#'
#' @description
#' Pools a balanced random subset of gated events, trains a hexagonal self-organizing map, maps retained samples to its nodes, and writes counts and mapped FCS files.
#'
#' @param gated_fcs Named list of samples, each list containing name and a numeric event matrix data.
#' @param m Number of SOM nodes (default 2000)
#' @param n_hclust Reserved argument, default 300; currently unused by the implementation.
#' @param n_cells_sub Number of cells to subsample (default 300,000)
#' @param out_path Output directory path (default './')
#'
#'
#' @export
MBR_som <- function(gated_fcs = NULL, m = 2e3, n_hclust = 300,
                    n_cells_sub = 3e5, out_path = './') {

  # ================================================================
  # Preparation: Output Directory
  # ================================================================
  cat("==========================================\n")
  cat("Starting Analysis Pipeline\n")
  cat("==========================================\n")

  # Create output directory
  if (!dir.exists(out_path)) {
    dir.create(out_path, recursive = TRUE)
    cat("\u2713 Created output directory:", out_path, "\n")
  }

  # ================================================================
  # Stage 1: Sample Quality Control
  # ================================================================
  cat("\nStage 1: Sample Quality Control\n")
  cat("----------------------------------------\n")

  min_cells <- 2e5
  original_sample_count <- length(gated_fcs)

  # Check sample cell counts
  cell_counts <- sapply(gated_fcs, function(x) nrow(x$data))
  low_count_samples <- cell_counts < min_cells

  if (any(low_count_samples)) {
    dropped_samples <- names(gated_fcs)[low_count_samples]
    warning(paste0("Dropping low cell count samples: ",
                   paste(dropped_samples, collapse = ", "),
                   " (cell count < ", format(min_cells, scientific = FALSE), ")"))
    cat("\u26a0 Dropped", length(dropped_samples), "samples\n")
  }

  # Filter samples
  gated_fcs <- gated_fcs[!low_count_samples]
  cat("\u2713 Retained", length(gated_fcs), "/", original_sample_count, "samples\n")

  if (!length(gated_fcs)) stop("No samples meet the 200000-event minimum.")

  # ================================================================
  # Stage 2: Random Sampling
  # ================================================================
  cat("\nStage 2: Performing Random Sampling\n")
  cat("----------------------------------------\n")
  cat("Target sampling size:", format(n_cells_sub, scientific = FALSE), "cells\n")

  # Execute sampling
  sub_sample <- downsample(
    gated_fcs,
    n = n_cells_sub,
    samples = names(gated_fcs)
  )
  sub_sample <- list(name = "random sample", data = sub_sample)

  cat("\u2713 Sampling completed, obtained", nrow(sub_sample$data), "cells\n")

  # ================================================================
  # Stage 3: SOM Training
  # ================================================================
  cat("\nStage 3: Training Self-Organizing Map (SOM)\n")
  cat("----------------------------------------\n")
  cat("Calculating SOM grid size...\n")

  # Calculate SOM parameters
  cells_per_node <- nrow(sub_sample$data) / m
  cat("Expected cells per node:", round(cells_per_node, 2), "\n")

  # Train SOM
  sub_sample <- compute_som(sub_sample, n_cells = cells_per_node)

  cat("\u2713 SOM training completed\n")
  cat("SOM grid size:", dim(sub_sample$som$grid$pts)[1], "nodes\n")

  # ================================================================
  # Stage 4: Save SOM Information
  # ================================================================
  cat("\nStage 4: Saving SOM Codebook Information\n")
  cat("----------------------------------------\n")

  # Save SOM codebook
  cohonen_information <- as.data.frame(sub_sample[["som"]][["codes"]][[1]])

  som_file <- file.path(out_path, "cluster_information.csv")
  write.csv(cohonen_information, som_file)

  cat("\u2713 SOM codebook saved to:", som_file, "\n")
  cat("Codebook dimensions:", dim(cohonen_information), "\n")

  # ================================================================
  # Stage 5: Mapping to SOM
  # ================================================================
  cat("\nStage 5: Mapping All Samples to Trained SOM\n")
  cat("----------------------------------------\n")

  n_subset_large <- 1e10  # Use all cells

  cat("Mapping", length(gated_fcs), "samples to SOM...\n")

  # Map to trained SOM (with progress bar)
  SOM_fcs <- pbapply::pblapply(gated_fcs, map_som,
                               trained = sub_sample$som,
                               n_subset = n_subset_large)

  cat("\u2713 SOM mapping completed\n")

  # ================================================================
  # Stage 6: Cluster Assignment
  # ================================================================
  cat("\nStage 6: Assigning Cluster Labels\n")
  cat("----------------------------------------\n")

  # Create cluster number matrix
  cluster_number <- as.matrix(as.list(1:2025))
  colnames(cluster_number) <- as.character(2025)

  cat("Total clusters:", nrow(cluster_number), "\n")
  cat("Assigning cluster labels...\n")

  # Assign clusters (with progress bar)
  SOM_fcs <- pbapply::pblapply(SOM_fcs, assign_clusters, clusters = cluster_number)

  cat("\u2713 Cluster assignment completed\n")

  # ================================================================
  # Stage 7: Counting and Statistics
  # ================================================================
  cat("\nStage 7: Computing Cluster Statistics\n")
  cat("----------------------------------------\n")

  cat("Calculating cell counts for each cluster...\n")

  # Count observations (with progress bar)
  SOM_fcs <- pbapply::pblapply(SOM_fcs, count_observations,
                               clusters = colnames(cluster_number))

  # Get count tables
  count_tables <- get_counts(SOM_fcs)

  cat("\u2713 Counting completed\n")

  # ================================================================
  # Stage 8: Save Count Results
  # ================================================================
  cat("\nStage 8: Saving Count Results\n")
  cat("----------------------------------------\n")

  # Save count tables
  tmp_path <- file.path(out_path, "count_tables")
  if (!dir.exists(tmp_path)) dir.create(tmp_path, recursive = TRUE)

  count_file <- file.path(tmp_path, "count_table_SOM.save")
  save(count_tables, file = count_file, compress = TRUE)

  cat("\u2713 Raw count table saved to:", count_file, "\n")

  # Process and normalize counts
  raw_count_table <- as.data.frame(t(count_tables[[1]]))
  count_table_2025 <- raw_count_table / rowSums(raw_count_table)

  # Check for rare clusters
  rare_clusters <- sum(apply(count_table_2025 < 1e-4, 2, all))
  if (rare_clusters > 0) {
    warning(paste0(rare_clusters, ' clusters are below 0.01% of cells in ALL samples!'))
    cat("Found", rare_clusters, "clusters below 0.01% in every sample\n")
  } else {
    cat("\u2713 No clusters found below 0.01% in all samples (0 rare clusters)\n")
  }

  # Save normalized count table
  som_file <- file.path(out_path, "SOM.csv")
  write.csv(count_table_2025, som_file)

  cat("\u2713 Normalized count table saved to:", som_file, "\n")

  # ================================================================
  # Stage 9: Export FCS Files
  # ================================================================
  cat("\nStage 9: Exporting Mapped FCS Files\n")
  cat("----------------------------------------\n")

  # Create output directory
  dir_name <- file.path(out_path, "MappedFCS")
  if (!dir.exists(dir_name)) {
    dir.create(dir_name, recursive = TRUE)
    cat("\u2713 Created FCS output directory:", dir_name, "\n")
  }

  cat("Exporting", length(SOM_fcs), "FCS files...\n")

  # Export each sample's FCS file
  pb <- txtProgressBar(min = 0, max = length(SOM_fcs), style = 3)

  for (i in seq_along(SOM_fcs)) {
    # Prepare data
    fcs_data <- as.data.frame(SOM_fcs[[i]]$data)
    fcs_data$classes <- as.numeric(SOM_fcs[[i]]$classes)

    # Convert to matrix and create flow cytometry frame
    fcs_matrix <- as.matrix(fcs_data)
    colnames(fcs_matrix) <- names(fcs_data)
    current_frame <- flowFrame(fcs_matrix)

    # Generate filename and save
    file_name <- file.path(dir_name, SOM_fcs[[i]]$name)
    write.FCS(current_frame, file_name)

    # Update progress bar
    setTxtProgressBar(pb, i)
  }

  close(pb)
  cat("\n\u2713 FCS file export completed\n")

  cat("\n==========================================\n")
  cat("Analysis Pipeline Completed!\n")
  cat("==========================================\n")
  cat("Processed samples:", length(SOM_fcs), "\n")
  cat("Total clusters:", nrow(cluster_number), "\n")
  cat("Output directory:", out_path, "\n")
  cat("\nMain output files:\n")
  cat("- SOM codebook:", file.path(out_path, "cluster_information.csv"), "\n")
  cat("- Normalized counts:", file.path(out_path, "SOM.csv"), "\n")
  cat("- Raw counts:", file.path(tmp_path, "count_table_SOM.save"), "\n")
  cat("- Mapped FCS files:", dir_name, "\n")
  cat("==========================================\n")

  return(list(
    cohonen_information = cohonen_information,
    raw_count_table = raw_count_table,
    count_table = count_table_2025
  ))
}

#' Statistical test
#'
#' @description
#' Computes one group-comparison p-value per feature, optionally adjusts for multiple testing, and exports significant feature columns.
#'
#' @param data A numeric data frame with samples in rows and features in columns.
#' @param meta_data A data frame containing sample metadata. Must include the grouping column.
#' @param group_col A character string indicating the column in `meta_data` that defines group labels.
#' @param test_type Type of statistical test to apply. One of `"wilcox"`, `"kruskal"`, `"anova"`, `"t.test"`.
#' @param cutoff Numeric p-value cutoff for significance (default: 0.05).
#' @param correction Method for multiple testing correction. One of `"none"` (default), `"fdr"`, `"bonferroni"`, `"BH"`.
#' @param out_path Path to directory where result CSV files will be saved (default: current directory).
#'
#' @export
MBR_stat <- function(data = NULL, meta_data = NULL, group_col = NULL,
                     test_type = 'wilcox', cutoff = 0.05,
                     correction = 'none', out_path = './') {

  inputs <- .mbr_inputs(data, meta_data, group_col)
  data <- inputs$data
  meta_data <- inputs$meta_data

  cat(paste0(test_type, " test\n"))

  if (length(unlist(meta_data[group_col])) != nrow(data)) {
    stop("Error: The lengths of group_labels and number of samples are not equal.")
  }

  group_vec <- as.vector(unlist(meta_data[group_col]))

  if (!test_type %in% c('wilcox', 'kruskal', 'anova', 't.test', 'ttest')) {
    stop("test_type must be 'wilcox', 'kruskal', 'anova', 't.test' or 'ttest'.")
  }

  if (!correction %in% c('none', 'fdr', 'bonferroni', 'BH')) {
    stop("Unknown multiple-testing correction.")
  }

  group_vec <- factor(group_vec)

  if (nlevels(group_vec) < 2L ||
      (test_type %in% c('wilcox', 't.test', 'ttest') && nlevels(group_vec) != 2L)) {
    stop("The selected test requires two groups (or at least two for kruskal/anova).")
  }

  constant_features <- vapply(data, function(values) length(unique(values)) < 2L, logical(1))

  if (any(constant_features)) {
    warning(sprintf(
      "Skipped %d constant features with no variation; their p-values remain NA and they are excluded from significance selection.",
      sum(constant_features)
    ), call. = FALSE)
  }

  original_pvalues <- vapply(names(data), function(feature) {

    if (constant_features[[feature]]) return(NA_real_)

    t <- data[[feature]]

    tryCatch({

      value <- switch(test_type,
        wilcox = stats::wilcox.test(t ~ group_vec, exact = FALSE)$p.value,
        kruskal = stats::kruskal.test(t ~ group_vec)$p.value,
        anova = stats::anova(stats::aov(t ~ group_vec))$`Pr(>F)`[1],
        t.test = stats::t.test(t ~ group_vec)$p.value,
        ttest = stats::t.test(t ~ group_vec)$p.value
      )

      if (!is.finite(value)) {
        stop("Test returned a non-finite p-value.")
      }

      value

    }, error = function(e) {

      warning(sprintf("Feature '%s': %s", feature, conditionMessage(e)),
              call. = FALSE)

      NA_real_
    })

  }, numeric(1))

  pvalue_df <- data.frame(p.value = original_pvalues)

  if (correction == 'fdr') {
    pvalue_df$p.adj <- p.adjust(original_pvalues, method = 'fdr')

  } else if (correction == 'bonferroni') {
    pvalue_df$p.adj <- p.adjust(original_pvalues, method = 'bonferroni')

  } else if (correction == 'BH') {
    pvalue_df$p.adj <- p.adjust(original_pvalues, method = 'BH')

  } else {
    pvalue_df$p.adj <- original_pvalues
  }

  significant_data <- data[, !is.na(pvalue_df$p.adj) &
                             pvalue_df$p.adj < cutoff, drop = FALSE]

  write.csv(pvalue_df, file.path(out_path, "pvalue.csv"), row.names = TRUE)
  write.csv(significant_data, file.path(out_path, "significant.csv"), row.names = TRUE)

  cat("\tDONE.\n")

  return(list(
    pvalue = pvalue_df,
    significant = significant_data
  ))
}

#' Circular heatmap
#'
#' @description
#' Displays feature means by group and a separate track for the difference between the two group means.
#'
#' @param data A numeric data frame with samples in rows and features in columns.
#' @param meta_data A data frame containing sample metadata with grouping information.
#' @param group_col Character string specifying the column name in `meta_data` for group labels.
#' @param out_path Directory path to save the output PDF file (default: current directory).
#' @param width Width of the output PDF.
#' @param height Height of the output PDF.
#' @param point_colors Vector of two colors for points representing positive and negative mean differences (default: c("red", "blue")).
#' @param cell_colors Vector of three colors for the heatmap gradient (default: c("blue", "white", "red")).
#'
#'
#' @export
MBR_circle <- function(data = NULL, meta_data = NULL, group_col = NULL,
                     out_path = './', width = 8, height = 8,
                     point_colors = c("red", "blue"),
                     cell_colors = c("blue", "white", "red")) {
  inputs <- .mbr_inputs(data, meta_data, group_col)
  data <- inputs$data
  meta_data <- inputs$meta_data

  group_labels <- as.vector(unlist(meta_data[group_col]))
  mean_table <- t(sapply(data, function(t) tapply(t, group_labels, mean)))
  # Validate case collisions, then format only complete cluster IDs for display.
  cluster_ids <- .mbr_cluster_ids(rownames(mean_table), "Feature names")
  rownames(mean_table) <- sub('^V([0-9]+)$', 'Bin\\1', cluster_ids)
  min_val <- min(mean_table)
  max_val <- max(mean_table)
  mean_val <- (min_val + max_val) / 2
  col_fun = colorRamp2(c(min_val, mean_val, max_val), cell_colors)

  circos.clear()
  circos.par(gap.degree = 10)
  circos.heatmap(mean_table, col = col_fun,
                 rownames.side = 'outside')

  mean_diff <- as.vector(mean_table[,1] - mean_table[,2])
  suppressMessages({
    circos.track(ylim = range(mean_diff),
                 track.height = 0.1,
                 panel.fun = function(x, y) {
                   y = mean_diff[CELL_META$row_order]
                   circos.points(seq_along(y) - 0.5, y, cex = 1,
                                 col = ifelse(y > 0, point_colors[1], point_colors[2]))
                   circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "grey")
                   # add column label
                   circos.text(c(-0.8, -0.8),
                               c(1.8 * max(mean_diff) - 1.8 * min(mean_diff),
                                 3.2 * max(mean_diff) - 3.2 * min(mean_diff)),
                               rev(unique(group_labels)),
                               facing = 'downward',
                               cex = 0.5)
                 },
                 cell.padding = c(0.02, 0, 0.02, 0))
  })

  pdf(paste0(out_path, '/circle.pdf'), height = height, width = width)
  circos.clear()
  circos.par(gap.degree = 10)
  circos.heatmap(mean_table, col = col_fun,
                 rownames.side = 'outside')

  mean_diff <- as.vector(mean_table[,1] - mean_table[,2])
  suppressMessages({
    circos.track(ylim = range(mean_diff),
                 track.height = 0.1,
                 panel.fun = function(x, y) {
                   y = mean_diff[CELL_META$row_order]
                   circos.points(seq_along(y) - 0.5, y, cex = 1,
                                 col = ifelse(y > 0, point_colors[1], point_colors[2]))
                   circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "grey")
                   # add column label
                   circos.text(c(-0.8, -0.8),
                               c(1.8 * max(mean_diff) - 1.8 * min(mean_diff),
                                 3.2 * max(mean_diff) - 3.2 * min(mean_diff)),
                               rev(unique(group_labels)),
                               facing = 'downward',
                               cex = 0.5)
                 },
                 cell.padding = c(0.02, 0, 0.02, 0))
  })
  invisible(dev.off())
  cat("\tDONE.\n")
}

#' Violin plot
#'
#' @description
#' Draws a violin and box plot for one cluster and annotates a supplied raw or adjusted p-value.
#'
#' @param data A numeric data frame with samples in rows and lowercase v-prefixed cluster columns.
#' @param meta_data A data frame containing sample metadata, including grouping information.
#' @param group_col Character string specifying the column in `meta_data` with group labels.
#' @param pvalue_data A data frame or matrix containing p-values for each cluster.
#' @param p Character string specifying which p-value column to use from `pvalue_data`. Must be either `"p.value"` or `"p.adj"`. Default is `"p.value"`.
#' @param cluster Numeric value specifying which cluster to plot.
#' @param out_path Directory path to save the output PDF (default: current directory).
#' @param colors Group fill colors; default c("#FFADAD", "#DEDAF4").
#' @param width Width of the saved PDF.
#' @param height Height of the saved PDF.
#'
#' @export
MBR_violin <- function(data = NULL, meta_data = NULL, group_col = NULL,
                       pvalue_data = NULL, p = "p.value",
                       cluster = 1, out_path = './',
                       colors = c('#FFADAD', '#DEDAF4'),
                       width = 4, height = 4) {
  inputs <- .mbr_inputs(data, meta_data, group_col)
  data <- inputs$data
  meta_data <- inputs$meta_data



  # Match data and p-value labels independently; their prefix case may differ.
  cluster_id <- paste0("V", cluster)
  data_index <- .mbr_cluster_match(cluster_id, names(data), "Data columns")
  pvalue_index <- .mbr_cluster_match(cluster_id, rownames(pvalue_data), "P-value rows")
  cluster_values <- data[[data_index]]
  p_val_numeric <- as.numeric(pvalue_data[pvalue_index, p])
  p_val_label <- format(p_val_numeric, digits = 3, scientific = TRUE)

  plot_df <- data.frame(
    group = factor(meta_data[[group_col]]),
    abundance = cluster_values
  )

  max_val <- max(plot_df$abundance, na.rm = TRUE)
  annotation_y_pos <- max_val * 1.1

  p_plot <- ggplot(data = plot_df, aes(x = group, y = abundance)) +
    geom_violin(aes(fill = group), alpha = 0.5) +
    geom_boxplot(aes(fill = group), width = 0.1, outlier.shape = NA) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = colors) +
    labs(title = paste0('Cluster ', cluster),
         x = 'Group',
         y = 'Relative Abundance',
         fill = 'Group') +
    annotate("text", x = 1.5, y = annotation_y_pos,
             label = paste0(p, " = ", p_val_label),
             size = 4, color = "black")

  print(p_plot)

  file_name <- paste0('violin_cluster_', cluster, '_', p, '.pdf')
  pdf(file.path(out_path, file_name), width = width, height = height)
  print(p_plot)
  dev.off()

  cat("Done: ", file_name, "\n")
}

#' Ordinate Bray-Curtis dissimilarity
#'
#' @description
#' Calculates Bray-Curtis dissimilarities and display plot.
#'
#' @param data A numeric data frame or matrix of features (samples as rows).
#' @param out_path Output directory path to save plots and results (default: "./").
#' @param test Statistical test for group differences: one of "wilcox", "ttest", "kruskal", or "anova" (default: "wilcox").
#' @param meta_data A data frame of metadata corresponding to samples.
#' @param group_name The column name in `meta_data` to define groups for comparison.
#' @param colors Vector of colors for groups in plots (default: c('#E41A1C', '#377EB8')).
#' @param width Width of output plot.
#' @param height Height of output plot.
#'
#' @export
MBR_beta <- function(data, out_path = './', test = 'wilcox',
                   meta_data = NULL, group_name = NULL,
                   colors = c('#E41A1C', '#377EB8'),
                   width = 5, height = 5) {
  inputs <- .mbr_inputs(data, meta_data, group_name)
  data <- inputs$data
  meta_data <- inputs$meta_data


  if (!test %in% c('wilcox', 'ttest', 'kruskal', 'anova')) {
    stop("test should be one of 'wilcox', 'ttest', 'kruskal', or 'anova'")
  }

  df <- cbind(meta_data[group_name], data)
  colnames(df)[1] <- 'Health.State'
  rownames(df) <- NULL
  groups <- data.frame(sample = paste('sample', 1:nrow(df), sep = ''), group = df$Health.State)
  write.table(groups, paste0(out_path, '/.group.txt'), row.names = FALSE, sep = '\t', quote = FALSE)

  df$Health.State <- groups$sample
  rownames(df) <- df$Health.State
  dataT <- df[, -1]

  dist <- vegdist(dataT, method = "bray")
  adist <- as.dist(dist)

  sd <- groups
  rownames(sd) <- sd$sample
  sd$group <- factor(sd$group, levels = unique(sd$group))

  dist <- as.matrix(dist)[sd$sample, sd$sample]

  pcoa <- cmdscale(dist, k = 3, eig = TRUE)
  pc12 <- as.data.frame(pcoa$points[, c(1, 2)])
  colnames(pc12) <- c("pc_x", "pc_y")
  pc12$sample <- rownames(pc12)

  pc <- round(pcoa$eig / sum(pcoa$eig) * 100, 2)
  pc12 <- merge(pc12, sd, by = "sample")
  pc12$group <- factor(pc12$group, levels = levels(sd$group))

  cp <- combn(levels(pc12$group), 2)
  comp <- split(cp, col(cp))

  ADONIS <- suppressMessages(adonis2(adist ~ group, data = sd))
  TEST <- ADONIS$`Pr(>F)`[1]
  R2adonis <- round(ADONIS$R2[1], 3)

  sink(paste0(out_path, '/adonis.txt'))
  print(ADONIS)
  sink()

  p <- ggscatter(pc12, x = "pc_x", y = "pc_y",
                 color = "group", shape = "group", size = 3,
                 ellipse = TRUE, conf.int.level = 0.95,
                 mean.point = TRUE, star.plot = TRUE, star.plot.lty = 1) +
    geom_hline(yintercept = 0, color = '#B3B3B3') +
    geom_vline(xintercept = 0, color = '#B3B3B3') +
    scale_color_manual(values = colors) +
    ggtitle(bquote(bold("Bray Curtis") ~ "," ~ bold("Adonis R"^2) == .(R2adonis) ~ "," ~ bold("P") == .(TEST))) +
    theme(axis.title = element_blank(),
          legend.position = "top", legend.title = element_blank(),
          panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
          text = element_text(size = 12),
          plot.margin = unit(c(0, 0, 0, 0), 'cm'))

  ## test method
  if (test == "anova") {
    method_use <- "anova"
    stat_func <- function(formula, data) {
      aov_res <- aov(formula, data = data)
      summary(aov_res)[[1]][["Pr(>F)"]][1]
    }
    pval_y <- stat_func(pc_y ~ group, pc12)
    pval_x <- stat_func(pc_x ~ group, pc12)
  } else {
    method_map <- c(wilcox = "wilcox.test", ttest = "t.test", kruskal = "kruskal.test")
    method_use <- method_map[test]
    pval_y <- NULL
    pval_x <- NULL
  }

  pl <- ggboxplot(pc12, x = "group", y = "pc_y", fill = "group", palette = colors)

  if (test == "anova") {
    pl <- pl + stat_compare_means(method = "anova", label.y = max(pc12$pc_y) * 1.05)
  } else {
    pl <- pl + stat_compare_means(comparisons = comp, label = "p.signif", method = method_use)
  }

  pl <- pl +
    ylab(paste0("PCoA", 2, " (", pc[2], "%)")) +
    theme(panel.border = element_rect(color = "black", linewidth = 1.0, fill = NA),
          legend.position = "none",
          axis.text = element_text(face = 'bold'),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 12, angle = 60, hjust = 1, face = 'bold'),
          axis.title.y = element_text(size = 15, face = 'bold')) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'cm'))

  pt <- ggboxplot(pc12, x = "group", y = "pc_x", fill = "group", palette = colors) + coord_flip()

  if (test == "anova") {
    pt <- pt + stat_compare_means(method = "anova", label.x = 1.2)
  } else {
    pt <- pt + stat_compare_means(comparisons = comp, label = "p.signif", method = method_use)
  }

  pt <- pt +
    scale_x_discrete(limits = rev(levels(pc12$group))) +
    ylab(paste0("PCoA", 1, " (", pc[1], "%)")) +
    theme(panel.border = element_rect(color = "black", linewidth = 1.0, fill = NA),
          legend.position = "none",
          axis.text = element_text(size = 12, angle = 0, face = 'bold'),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 15, face = 'bold')) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'cm'))

  p0 <- ggplot() + theme(panel.background = element_blank(),
                         plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "lines"))

  p_final <- pl + p + p0 + pt +
    plot_layout(ncol = 2, nrow = 2, heights = c(4, 1), widths = c(1, 4)) +
    plot_annotation(theme = theme(plot.margin = margin()))

  ggsave(paste0(out_path, '/pcoa_', test, '.pdf'), plot = p_final, width = width, height = height)
  return(p_final)
}

#' Select features with recursive elimination
#'
#' @description
#' Uses cross-validated recursive feature elimination and produces feature importance and transformed abundance effect-size plots.
#'
#' @param data A data frame or matrix with features as columns and samples as rows.
#' @param out_path Output directory for saving results and plots (default: "./").
#' @param meta_data A data frame containing sample metadata.
#' @param group_name Column name in `meta_data` indicating the grouping variable (e.g., disease status).
#' @param nfolds_cv Number of cross-validation folds used in RFE (default: 5).
#' @param top_n_features Number of top features to display in the importance plot (default: 20).
#' @param rfe_size Largest requested subset size in 1:rfe_size (default 10); caret also evaluates the full predictor set.
#' @param ref_group Reference group name for effect size calculation (default: NULL).
#' @param colors A vector of two colors for group comparison plots (default: c('#E41A1C', '#377EB8')).
#' @param width Width of the output PDF plot.
#' @param height Height of the output PDF plot.
#'
#' @export
MBR_fs <- function(data = NULL, out_path = './',
                   meta_data = NULL,
                   group_name = NULL,
                   nfolds_cv = 5,
                   top_n_features = 20,
                   rfe_size = 10,
                   ref_group = NULL,
                   colors = c('#E41A1C', '#377EB8'),
                   width = 5, height = 5) {
  inputs <- .mbr_inputs(data, meta_data, group_name)
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("The randomForest backend is required for this analysis.")
  }
  data <- inputs$data
  meta_data <- inputs$meta_data


  set.seed(2025)

  df <- cbind(meta_data[group_name], data)
  colnames(df)[1] <- 'Health.State'
  rownames(df) <- NULL

  groups <- data.frame(sample = paste0('sample', seq_len(nrow(df))), group = df$Health.State)
  write.table(groups, file = file.path(out_path, '.group.txt'), row.names = FALSE, sep = '\t', quote = FALSE)

  df$Health.State <- groups$sample
  rownames(df) <- df$Health.State
  df <- df[, -1]
  df <- as.data.frame(t(df))
  df2 <- as.data.frame(t(df))
  df2$group <- factor(groups$group)

  control <- rfeControl(functions = rfFuncs, method = "cv", number = nfolds_cv)
  results <- rfe(df2[, 1:(ncol(df2) - 1)], df2[, ncol(df2)], sizes = 1:rfe_size, rfeControl = control)

  print(results)
  predictors(results)
  p1_out <- plot(results, type = c("g", "o"), main = 'Feature Selection', col = '#377EB8', lwd = 2)

  best_selection <- results$optVariables
  message("The number of selected features is ", length(best_selection))

  df2 <- df2[, c(best_selection, 'group')]
  write.table(df2, file = file.path(out_path, 'figure1.txt'), row.names = FALSE, sep = '\t', quote = FALSE)

  # Extract selected features from input data
  selected_features <- data[, best_selection, drop = FALSE]
  rownames(selected_features) <- rownames(data)

  varimp_all <- varImp(results)
  varimp_df <- varimp_all[best_selection, , drop = FALSE]

  top_n <- min(top_n_features, nrow(varimp_df))
  varimp_data <- data.frame(
    feature = rownames(varimp_df)[1:top_n],
    importance = varimp_df[1:top_n, 1]
  )
  top_features <- varimp_data$feature
  df2 <- df2[, c(top_features, 'group')]

  p2_out <- ggplot(tidyr::drop_na(varimp_data),
                   aes(x = reorder(feature, -importance), y = importance, fill = feature)) +
    geom_bar(stat = "identity") +
    labs(x = "Features", y = "Variable Importance") +
    geom_text(aes(label = round(importance, 2)), vjust = 1.6, color = "white", size = 4) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5)) +
    labs(title = "Variable Importance")

  esize <- df2 %>%
    reshape2::melt() %>%
    mutate(value = log(value + 1, 10)) %>%
    group_by(variable) %>%
    suppressMessages() %>%
    cohens_d(value ~ group, conf.level = 0.95, ci = TRUE, ref.group = ref_group) %>%
    arrange(desc(effsize))

  fig2 <- ggplot(esize, aes(y = reorder(variable, effsize), x = effsize, fill = variable)) +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.2) +
    geom_point(aes(x = effsize), size = 4, color = 'black') +
    geom_point(aes(x = effsize, color = ifelse(effsize > 0, 'positive', 'negative')), size = 3) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    scale_color_manual(values = c('positive' = colors[1], 'negative' = colors[2])) +
    labs(title = "Effect Size", x = "Cohen's d", y = '') +
    force_panelsizes(rows = 0.5, cols = 0.5)

  fig1 <- df2 %>%
    reshape2::melt() %>%
    mutate(value = log(value + 1, 10)) %>%
    ggplot(aes(x = value, y = fct_rev(factor(variable, esize$variable)), fill = group)) +
    geom_boxplot(position = 'dodge') +
    theme_minimal() +
    theme(legend.position = "top") +
    labs(x = "log10(abundance)", y = "Variable") +
    scale_fill_manual(values = colors) +
    force_panelsizes(rows = 0.5, cols = 0.5)

  fig3 <- merge(varimp_data, esize, by.x = 'feature', by.y = 'variable') %>%
    ggplot(aes(x = importance, y = fct_rev(factor(feature, esize$variable)))) +
    geom_bar(stat = "identity", color = 'black', aes(fill = ifelse(effsize > 0, 'positive', 'negative'))) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    scale_fill_manual(values = c('positive' = colors[1], 'negative' = colors[2])) +
    labs(title = "Variable Importance", y = '') +
    force_panelsizes(rows = 0.5, cols = 0.5)

  p3_out <- plot_grid(fig1, fig2, fig3, ncol = 3, align = 'h', axis = 't')

  pdf(file.path(out_path, 'feature_exploration.pdf'), width = width, height = height)
  print(p1_out)
  print(p2_out)
  print(p3_out)
  dev.off()

  show_plots <- function(plots) {
    for (plot in plots) {
      print(plot)
      readline(prompt = "Press [enter] to see the next plot")
    }
  }

  show_plots(list(p1_out, p2_out, p3_out))
  cat("\tDONE.\n")

  invisible(selected_features)
}

#' Show channel values for selected SOM phenotype
#'
#' @description
#' Display phenotype of the selected clusters in heatmap.
#'
#' @param data Data frame or matrix whose column names correspond to SOM cluster IDs (e.g., "V1", "V2", ...).
#' @param cohonen_information Data frame containing information per SOM cluster. Row names should match cluster IDs.
#' @param out_path Directory path where the heatmap PDF will be saved (default './').
#' @param color Color palette for the heatmap (default is reversed RdYlBu palette from RColorBrewer).
#' @param width Width of the output PDF.
#' @param height Height of the output PDF.
#' @param scale Scaling method for pheatmap: 'row', 'column', or 'none' (default 'row').
#' @param cluster_rows Logical, whether to cluster rows (default FALSE).
#' @param cluster_cols Logical, whether to cluster columns (default FALSE).
#' @param display_numbers Logical, whether to display values on heatmap cells (default TRUE).
#' @param boarder_color Color of the border around heatmap cells (default 'grey60').
#' @param legend (default TRUE).
#'
#'
#' @export
MBR_heatmap <- function(data = NULL, cohonen_information = NULL,
                      out_path = './',
                      color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                      width = 10, height = 10,
                      scale = 'row', cluster_rows = F, cluster_cols = F,
                      display_numbers = T, boarder_color = 'grey60',
                      legend = T) {
  # Preserve numeric abundance-column support, and match either prefix case
  # against the codebook without changing its original row names.
  cluster_ids <- .mbr_event_ids(colnames(data))
  cluster_index <- .mbr_cluster_match(cluster_ids, rownames(cohonen_information), "Codebook rows")
  df <- cohonen_information[cluster_index, , drop = FALSE]
  df_mat <- as.matrix(sapply(df, as.numeric))
  rownames(df_mat) <- rownames(df)
  p <- pheatmap(df, scale = scale, cluster_rows = cluster_rows,
           cluster_cols = cluster_cols, display_numbers = display_numbers,
           border_color = boarder_color, color = color, legend = legend)
  pdf(paste0(out_path, '/heatmap.pdf'), width = width, height = height)
  print(p)
  dev.off()
  print(p)
  cat("\tDONE.\n")
}


#' Visualize Mantel associations and feature correlations
#'
#' @description
#' Compares distance patterns of metadata blocks with feature data using linkET and overlays associations on a feature correlation heatmap.
#'
#' @param data Numeric data matrix or data frame (e.g., feature abundance data).
#' @param meta_data Metadata data frame containing clinical and demographic variables.
#' @param clinical_cols Character vector of column names in meta_data for clinical variables.
#' @param demographic_cols Character vector of column names in meta_data for demographic variables.
#' @param spec_select_names Named list defining labels for clinical and demographic variable groups (default list with A = "A", B = "B").
#' @param colors Color palette for the Pearson correlation heatmap (default RdBu palette with 11 colors).
#' @param out_path Directory path to save output PDF (default './').
#' @param width Width of output PDF.
#' @param height Height of output PDF.
#'
#' @export
MBR_mantel <- function(data = NULL, meta_data = NULL,
                       clinical_cols = NULL, demographic_cols = NULL,
                       spec_select_names = list(A = "A", B = "B"),
                       colors = RColorBrewer::brewer.pal(11, "RdBu"),
                       out_path = './', width = 8, height = 8) {

  inputs <- .mbr_inputs(data, meta_data, NULL)
  data <- inputs$data
  meta_data <- inputs$meta_data

  spec_select <- list()
  if (!is.null(spec_select_names$A)) {
    spec_select[[spec_select_names$A]] <- clinical_cols
  }
  if (!is.null(spec_select_names$B)) {
    spec_select[[spec_select_names$B]] <- demographic_cols
  }

  mantel <- suppressMessages(
    linkET::mantel_test(meta_data, data, spec_select = spec_select) %>%
      dplyr::mutate(
        rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                 labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
        pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                 labels = c("< 0.01", "0.01 - 0.05", ">= 0.05"))
      )
  )

  if (nrow(mantel) == 0) {
    stop("No Mantel associations were found.")
  }

  cor_tbl <- linkET::as_md_tbl(
    linkET::correlate(data),
    type = "lower",
    diag = FALSE
  )

  row_names <- attr(cor_tbl, "row_names")
  col_names <- attr(cor_tbl, "col_names")
  n <- length(col_names)

  mantel <- linkET::nice_curvature()(mantel, cor_tbl)

  plot_data <- as.data.frame(cor_tbl)
  plot_data$x <- match(plot_data$.colnames, col_names)
  plot_data$y <- match(plot_data$.rownames, rev(row_names))

  if (anyNA(plot_data$x) || anyNA(plot_data$y)) {
    stop("Unable to determine correlation plot coordinates.")
  }

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = x, y = y, fill = r)
  ) +
    ggplot2::geom_tile(
      colour = "grey35",
      linewidth = 0.25
    ) +
    linkET::geom_mark(
      mapping = ggplot2::aes(r = r, pvalue = p),
      sep = '\n',
      size = 1.8,
      sig_level = c(0.05, 0.01, 0.001),
      sig_thres = 0.05,
      color = "white"
    ) +
    ggplot2::scale_fill_gradientn(colours = colors) +
    ggplot2::scale_x_continuous(
      breaks = seq_len(n),
      labels = col_names
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq_len(n),
      labels = rev(row_names)
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )

  from <- mantel[[1]]
  to <- as.character(mantel[[2]])

  if (is.factor(from)) {
    groups <- levels(from)
    from <- as.character(from)
    groups <- groups[groups %in% from]
  } else {
    from <- as.character(from)
    groups <- unique(from[!is.na(from)])
  }

  m <- length(groups)

  if (m == 1) {
    x <- 0.5 + 0.82 * n
    y <- 0.5 + 0.7 * n
  } else if (m == 2) {
    x <- c(0.5 + 0.8 * n, 0.5 + 1.02 * n)
    y <- c(0.5 + 0.8 * n, 0.5 + 0.54 * n)
  } else {
    x <- seq(0.5 + 0.75 * n, 0.5 + 1.3 * n, length.out = m)
    y <- seq(0.5 + 0.9 * n, 0.5 + 0.3 * n, length.out = m)
  }

  names(x) <- groups
  names(y) <- groups

  row_names <- rev(row_names)
  xend <- setNames(n:1, row_names)
  yend <- setNames(seq_len(n), row_names)

  edge <- data.frame(
    x = unname(x[from]),
    y = unname(y[from]),
    xend = unname(xend[to]),
    yend = unname(yend[to]),
    curvature = mantel$.curvature,
    pd = mantel$pd,
    rd = mantel$rd
  )

  if (anyNA(edge[, c("x", "y", "xend", "yend", "curvature")])) {
    stop("Invalid coordinates in Mantel connections.")
  }

  for (curv in unique(edge$curvature)) {
    edge_part <- edge[edge$curvature == curv, , drop = FALSE]

    p <- p +
      ggplot2::geom_curve(
        data = edge_part,
        mapping = ggplot2::aes(
          x = x,
          y = y,
          xend = xend,
          yend = yend,
          colour = pd,
          linewidth = rd
        ),
        curvature = curv,
        ncp = 5,
        inherit.aes = FALSE
      )
  }

  nodes <- unique(data.frame(
    x = c(edge$x, edge$xend),
    y = c(edge$y, edge$yend)
  ))

  p <- p +
    ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(x = x, y = y),
      shape = 21,
      colour = "blue",
      fill = "red",
      size = 2,
      stroke = 0.5,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = data.frame(
        x = unname(x[groups]),
        y = unname(y[groups]),
        label = groups
      ),
      ggplot2::aes(x = x, y = y, label = label),
      hjust = 0,
      nudge_x = 0.5,
      size = 3.88,
      inherit.aes = FALSE
    ) +
    ggplot2::expand_limits(
      x = c(0.5, max(edge$x, na.rm = TRUE) + 0.2 * n),
      y = c(0.5, n + 0.5)
    ) +
    ggplot2::scale_linewidth_manual(values = c(0.5, 1, 2)) +
    ggplot2::scale_colour_manual(values = linkET::color_pal(3)) +
    ggplot2::guides(
      linewidth = ggplot2::guide_legend(
        title = "Mantel's r",
        override.aes = list(colour = "grey35"),
        order = 2
      ),
      colour = ggplot2::guide_legend(
        title = "Mantel's p",
        override.aes = list(linewidth = 3),
        order = 1
      ),
      fill = ggplot2::guide_colorbar(
        title = "Pearson's r",
        order = 3
      )
    )

  grDevices::pdf(
    file.path(out_path, "mantel.pdf"),
    width = width,
    height = height
  )

  print(p)
  grDevices::dev.off()

  print(p)
  cat("\tDONE.\n")
}



#' Evaluate random forest classification with ROC
#'
#' @description
#' Construct random forest classifier.
#'
#' @param data A data frame or matrix with features in columns and samples in rows.
#' @param meta_data A data frame containing metadata for the samples.
#' @param group_name The column name in `meta_data` indicating the grouping variable (e.g., disease vs control).
#' @param out_path Output directory for saving the ROC/confusion plot and group assignment file (default: "./").
#' @param reference_level The reference level of the group used for binary classification.
#' @param width Width of the output ROC/confusion plot PDF.
#' @param height Height of the output ROC/confusion plot PDF.
#' @param method Resampling method for training control in caret (default: "repeatedcv").
#' @param number Number of folds or resampling iterations (default: 5).
#' @param repeats Number of repeats for repeated cross-validation (default: 5).
#' @param test_data Optional numeric feature table for an independent test set, with the same feature names as data.
#' @param test_meta_data Metadata for test_data, with sample IDs and group_name.
#'
#'
#' @export
MBR_ml <- function(data = NULL, meta_data = NULL,
                   group_name = 'Group', out_path = './',
                   reference_level = 'A',
                   width = 5, height = 5,
                   method = "repeatedcv", number = 5, repeats = 5,
                   test_data = NULL, test_meta_data = NULL) {
  inputs <- .mbr_inputs(data, meta_data, group_name)
  .mbr_test_inputs(data, test_data, test_meta_data, group_name)
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("The randomForest backend is required for this analysis.")
  }
  data <- inputs$data
  meta_data <- inputs$meta_data


  set.seed(2025)

  # Create output directory if it doesn't exist
  if (!dir.exists(out_path)) {
    dir.create(out_path, recursive = TRUE)
    cat("Output directory created at:", out_path, "\n")
  }

  # Ensure data has no row names before processing
  rownames(data) <- NULL

  # Create groups dataframe with column named 'group' (fixed name)
  groups <- data.frame(
    sample = paste('sample', seq_len(nrow(data)), sep = ''),
    group = meta_data[[group_name]]
  )

  # Save group info
  group_file <- file.path(out_path, 'ml_group.txt')
  write.table(groups, group_file, row.names = FALSE, sep = '\t', quote = FALSE)
  cat("Group info saved to:", group_file, "\n")

  # Prepare data
  data <- cbind(data, Health.State = groups$sample)
  rownames(data) <- data$Health.State
  data <- data[, -ncol(data)]  # remove Health.State column
  df <- as.data.frame(t(data))
  train <- as.data.frame(t(df))
  train$group <- factor(groups$group)

  if (nlevels(train$group) != 2L || !reference_level %in% levels(train$group)) {
    stop("Supply exactly two classes and a reference_level present in the groups.")
  }
  if (any(make.names(levels(train$group)) != levels(train$group))) {
    stop("Class names must be valid R names for probability columns.")
  }
  train$group <- stats::relevel(train$group, ref = reference_level)
  if (!method %in% c("cv", "repeatedcv", "LOOCV")) {
    stop("Diagnostic predictions require method = 'cv', 'repeatedcv' or 'LOOCV'.")
  }

  # Train control
  fitControl <- caret::trainControl(
    method = method,
    number = number,
    repeats = if (method == "repeatedcv") repeats else NA,
    returnResamp = "final",
    classProbs = TRUE,
    savePredictions = "final",
    summaryFunction = caret::twoClassSummary
  )

  # Train RF
  rf <- caret::train(
    group ~ .,
    data = train,
    method = "rf",
    trControl = fitControl,
    metric = "ROC",
    verbose = FALSE
  )

  # Average repeated held-out predictions once per sample for the chosen tuning.
  diagnostic <- .mbr_evaluate(rf, reference_level, train, test_data, test_meta_data, group_name)
  rocs_train <- pROC::roc(response = diagnostic$truth,
                         predictor = diagnostic$probability,
                         levels = rev(levels(train$group)), direction = "<", quiet = TRUE)
  ci_auc_train <- suppressWarnings(round(as.numeric(pROC::ci.auc(rocs_train)), 3))
  ci_tb_train <- suppressWarnings(as.data.frame(pROC::ci.se(rocs_train)))
  ci_tb_train <- suppressWarnings(tibble::rownames_to_column(ci_tb_train, var = 'x'))
  ci_tb_train <- as.data.frame(sapply(ci_tb_train, as.numeric))
  names(ci_tb_train) <- c('x', 'low', 'mid', 'high')

  metrics_train <- caret::confusionMatrix(diagnostic$predicted, diagnostic$truth,
                                          positive = reference_level, mode = 'everything')


  g1_out <- pROC::ggroc(rocs_train, legacy.axes = TRUE) +
    ggplot2::coord_equal() +
    ggplot2::geom_ribbon(aes(x = 1 - x, ymin = low, ymax = high), data = ci_tb_train, alpha = 0.5, fill = 'lightblue') +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 'dashed', alpha = 0.7) +
    ggplot2::annotate("text", x = 0.5, y = 0.25, hjust = 0,
                      label = paste0('AUC: ', round(rocs_train$auc, 3), ' 95%CI: ',
                                     ci_auc_train[1], ' ~ ', ci_auc_train[3])) +
    ggplot2::annotate("text", x = 0.5, y = 0.20, hjust = 0,
                      label = paste0('Sensitivity: ', round(as.numeric(metrics_train$byClass[1]), 3))) +
    ggplot2::annotate("text", x = 0.5, y = 0.15, hjust = 0,
                      label = paste0('Specificity: ', round(as.numeric(metrics_train$byClass[2]), 3))) +
    ggplot2::annotate("text", x = 0.5, y = 0.10, hjust = 0,
                      label = paste0('F1: ', round(as.numeric(metrics_train$byClass[7]), 3))) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = if (is.null(test_data)) 'Cross-validation ROC' else 'Independent test ROC',
                  x = '1 - Specificity (false-positive rate)',
                  y = 'Sensitivity (true-positive rate)')


  # Save ROC plot
  plot_file <- file.path(
    out_path,
    if (is.null(test_data)) {
      "cross_validation_roc.pdf"
    } else {
      "independent_test_roc.pdf"
    }
  )

  pdf(plot_file, width = width, height = height)
  print(g1_out)
  dev.off()

  cat("ROC plot saved to:", plot_file, "\n")

  # Console plot
  print(g1_out)
  cat("\tDONE.\n")
  invisible(list(model = rf, predictions = diagnostic, roc = rocs_train,
                 confusion = metrics_train, plot = g1_out))
}

#' Evaluate random forest classification with a confusion heatmap
#'
#' @description
#' Confusion matrix for classification result.
#'
#' @param data A data frame or matrix with features as columns and samples as rows.
#' @param meta_data A data frame containing metadata for the samples.
#' @param group_name The column name in `meta_data` that indicates the grouping variable (e.g., case vs control).
#' @param out_path Directory path to save the confusion matrix plot and group assignment file (default: "./").
#' @param reference_level The reference level of the group used for binary classification (default: "A").
#' @param colors A vector of two colors used for the heatmap gradient (default: c('#00BFC4', '#F8766D')).
#' @param width Width of the output PDF plot.
#' @param height Height of the output PDF plot.
#' @param method Resampling method for training control in caret (default: "repeatedcv").
#' @param number Number of folds or resampling iterations (default: 5).
#' @param repeats Number of repeats for repeated cross-validation (default: 5).
#' @param test_data Optional numeric feature table for an independent test set, with the same feature names as data.
#' @param test_meta_data Metadata for test_data, with sample IDs and group_name. Supply both test arguments together.
#'
#'
#' @import caret
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
MBR_conf <- function(data = NULL, meta_data = NULL,
                     group_name = 'Group', out_path = './',
                     reference_level = 'A',
                     colors = c('#00BFC4', '#F8766D'),
                     width = 5, height = 5,
                     method = "repeatedcv", number = 5, repeats = 5,
                   test_data = NULL, test_meta_data = NULL) {
  inputs <- .mbr_inputs(data, meta_data, group_name)
  .mbr_test_inputs(data, test_data, test_meta_data, group_name)
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("The randomForest backend is required for this analysis.")
  }
  data <- inputs$data
  meta_data <- inputs$meta_data

  # Create output directory if it doesn't exist
  if (!dir.exists(out_path)) {
    dir.create(out_path, recursive = TRUE)
    cat("Output directory created at:", out_path, "\n")
  }

  # Generate sample names
  rownames(data) <- NULL
  groups <- data.frame(sample = paste0('sample', seq_len(nrow(data))),
                       group = meta_data[[group_name]])

  group_file <- file.path(out_path, 'ml_group.txt')
  write.table(groups, group_file, row.names = FALSE, sep = '\t', quote = FALSE)
  cat("Group info saved to:", group_file, "\n")

  # Prepare data
  data <- cbind(data, Health.State = groups$sample)
  rownames(data) <- data$Health.State
  data <- data[, -ncol(data)]
  df <- as.data.frame(t(data))
  train <- as.data.frame(t(df))
  train$group <- factor(groups$group)

  if (nlevels(train$group) != 2L || !reference_level %in% levels(train$group)) {
    stop("Supply exactly two classes and a reference_level present in the groups.")
  }
  if (any(make.names(levels(train$group)) != levels(train$group))) {
    stop("Class names must be valid R names for probability columns.")
  }
  train$group <- stats::relevel(train$group, ref = reference_level)
  if (!method %in% c("cv", "repeatedcv", "LOOCV")) {
    stop("Diagnostic predictions require method = 'cv', 'repeatedcv' or 'LOOCV'.")
  }

  # Train control
  fitControl <- caret::trainControl(
    method = method,
    number = number,
    repeats = if (method == "repeatedcv") repeats else NA,
    returnResamp = "final",
    classProbs = TRUE,
    savePredictions = "final",
    summaryFunction = caret::twoClassSummary
  )

  # Train model
  rf <- caret::train(
    group ~ .,
    data = train,
    method = "rf",
    trControl = fitControl,
    metric = "ROC",
    verbose = FALSE
  )

  # Confusion matrix
  diagnostic <- .mbr_evaluate(rf, reference_level, train, test_data, test_meta_data, group_name)
  pred <- diagnostic$predicted
  truth <- diagnostic$truth
  cm <- caret::confusionMatrix(pred, truth, positive = reference_level)
  cm_table <- as.data.frame(cm$table)

  # Plot

  g2_out <- ggplot(cm_table, aes(x = Reference, y = Prediction)) +
    geom_tile(aes(fill = log(Freq + 1))) +
    geom_text(aes(label = Freq)) +
    scale_fill_gradient2(low = colors[1], high = colors[2],
                         midpoint = mean(log(cm_table$Freq + 1))) +
    coord_equal() +
    theme_minimal() +
    labs(title = if (is.null(test_data)) 'Cross-validation confusion matrix' else 'Independent test confusion matrix',
         x = 'True Label',
         y = 'Predicted Label',
         fill = 'Log(Count)') +
    theme(legend.position = 'right')

  # Save confusion matrix plot
  plot_file <- file.path(
    out_path,
    if (is.null(test_data)) {
      "cross_validation_confusion_matrix.pdf"
    } else {
      "independent_test_confusion_matrix.pdf"
    }
  )

  pdf(file = plot_file, width = width, height = height)
  print(g2_out)
  dev.off()

  cat("Confusion matrix plot saved to:", plot_file, "\n")
  invisible(list(model = rf, predictions = diagnostic, confusion = cm))
}


#' Reclustering
#'
#' @description
#' Recluster existing cluster assignments into a specified target number of clusters.
#'
#' @param data A data frame containing numeric variables for clustering.
#' @param num_clusters Integer specifying the number of clusters to cut the dendrogram into.
#'
#'
#' @importFrom stats dist hclust cutree
#' @export
MBR_reclustering <- function(data, num_clusters) {
  if (missing(data) || !is.data.frame(data)) {
    stop("Please provide a valid data frame as 'data'.")
  }

  features <- data[, sapply(data, is.numeric)]
  if (ncol(features) == 0) {
    stop("No numeric columns found in the input data.")
  }

  features_scaled <- base::scale(features)
  dist_matrix <- dist(features_scaled, method = "euclidean")
  hc <- hclust(dist_matrix, method = "ward.D2")

  cluster_assignments <- cutree(hc, k = num_clusters)

  reclustered_information <- data
  reclustered_information$Cluster <- cluster_assignments

  cat("Reclustering complete. 'reclustered_information' created with dimensions:\n")
  print(dim(reclustered_information))

  invisible(reclustered_information)
}


#' Read FCS files from a directory
#'
#' @description
#' Reads matching files into a flowCore flowSet for event processing.
#'
#' @param rawdata_path Directory containing input FCS files.
#'
#'
#' @export
MBR_read <- function(rawdata_path) {
  # Check if the path exists
  if (!dir.exists(rawdata_path)) {
    stop("Directory does not exist: ", rawdata_path)
  }

  # List FCS files
  fcs_path <- list.files(rawdata_path, pattern = "fcs$", full.names = TRUE)

  if (length(fcs_path) == 0) {
    stop("No FCS files found in directory: ", rawdata_path)
  }

  # Read FCS files as flowSet
  fcs_files <- read.flowSet(files = fcs_path,
                            column.pattern = "\\*|Bits|Drop",
                            invert.pattern = TRUE,
                            alter.names = FALSE)

  return(fcs_files)
}

#' Transform and rename events from one FCS sample
#'
#' @description
#' Extracts one flowSet sample, transforms its channel values, and optionally renames columns.
#'
#' @param fcs_files A flowSet containing the samples to process.
#' @param file_index Index of the sample to process (default: 1).
#' @param transformation Vectorized function applied to event values.
#' @param column_mapping Optional named vector mapping original channel names to new names.
#'
#'
#' @export
MBR_process <- function(fcs_files,
                              file_index = 1,
                              transformation = function(x) 10^((4 * x) / 65000),
                              column_mapping = NULL) {

  # Check if file_index is valid
  if (file_index > length(fcs_files) || file_index < 1) {
    stop("Invalid file index. Must be between 1 and ", length(fcs_files))
  }

  # Get file name from flowSet
  file_name <- sampleNames(fcs_files)[file_index]

  # Extract data from specified file
  dat <- exprs(fcs_files@frames[[file_name]]) %>%
    as.data.frame()

  # Apply transformation to all columns except 'classes' if it exists
  if ("classes" %in% colnames(dat)) {
    dat <- dat %>%
      mutate_at(vars(-'classes'), transformation)
  } else {
    dat <- dat %>%
      mutate_all(transformation)
  }

  # Rename columns if mapping is provided
  if (!is.null(column_mapping)) {
    new_names <- colnames(dat)
    for (old_name in names(column_mapping)) {
      if (old_name %in% colnames(dat)) {
        new_names[which(colnames(dat) == old_name)] <- column_mapping[old_name]
      }
    }
    colnames(dat) <- new_names
  }

  return(dat)
}

#' Select and transform phenotype for plotting
#'
#' @description
#' Filters a phenotype table by row name, renames channels, and transforms all remaining columns.
#'
#' @param bins Numeric phenotype matrix or data frame with cluster row names.
#' @param selected_rows Character vector of row names to retain.
#' @param transformation Vectorized numeric transformation applied to every column.
#' @param column_mapping Optional named character vector mapping old channel names to new names.
#'
#'
#' @export
MBR_prepare <- function(bins,
                         selected_rows,
                         transformation = function(x) 10^((4 * x) / 65000),
                         column_mapping = NULL) {

  # Validate requests but retain the original input row order and spelling.
  selected_index <- .mbr_cluster_match(selected_rows, rownames(bins), "Bin rows")
  bins <- bins[seq_len(nrow(bins)) %in% selected_index, , drop = FALSE]

  # Convert to dataframe if not already
  bins <- as.data.frame(bins)

  # Rename columns if mapping is provided
  if (!is.null(column_mapping)) {
    new_names <- colnames(bins)
    for (old_name in names(column_mapping)) {
      if (old_name %in% colnames(bins)) {
        new_names[which(colnames(bins) == old_name)] <- column_mapping[old_name]
      }
    }
    colnames(bins) <- new_names
  }

  # Apply transformation
  bins <- bins %>%
    mutate_all(transformation)

  return(bins)
}

#' Density plot
#'
#' @description
#' Builds a two-channel density plot with optional channel labels.
#'
#' @param dat A data frame processed with MBR_process, containing the flow cytometry data.
#' @param bins A data frame processed with MBR_prepare, containing information about data bins, typically used for displaying cluster centroids. Defaults to `NULL`.
#' @param x Character string, the name of the column in `dat` to be used for the x-axis.
#' @param y Character string, the name of the column in `dat` to be used for the y-axis.
#' @param selected_rows Character vector, names of clusters of interest to highlight. These should correspond to values in the 'classes' column of `dat`. Defaults to `NULL`.
#' @param x_limits Numeric vector of length 2, specifying the lower and upper limits for the x-axis. Defaults to `c(0.9, 11000)`.
#' @param y_limits Numeric vector of length 2, specifying the lower and upper limits for the y-axis. Defaults to `c(0.9, 11000)`.
#' @param hex_bins Integer, the number of bins to use for the hexagonal binning. Defaults to 100.
#' @param point_color Character string, the color for highlighted points. Defaults to "grey20".
#' @param point_alpha Numeric, the transparency level for highlighted points (0 = transparent, 1 = opaque). Defaults to 0.3.
#' @param point_size Numeric, the size of highlighted points. Defaults to 0.5.
#' @param label_color Character string, the color for bin labels. Defaults to "white".
#' @param label_size Numeric, the size of bin labels. Defaults to 3.
#' @param theme_family Font family for the plot theme; default sans.
#'
#'
#' @export
MBR_flow_plot <- function(dat,
                             bins = NULL,
                             x,
                             y,
                             selected_rows = NULL,
                             x_limits = c(0.9, 11000),
                             y_limits = c(0.9, 11000),
                             hex_bins = 100,
                             point_color = "grey20",
                             point_alpha = 0.3,
                             point_size = 0.5,
                             label_color = "white",
                             label_size = 3,
                             theme_family = "sans") {

  # Create the base plot
  p <- dat %>%
    ggplot(aes(x = .data[[x]], y = .data[[y]])) +
    geom_hex(bins = hex_bins) +
    scale_fill_viridis(discrete = FALSE, trans = 'log') +
    scale_x_log10(breaks = c(1, 10, 100, 1000, 10000),
                  labels = trans_format("log10", scales::math_format(10^.x))) +
    scale_y_log10(breaks = c(1, 10, 100, 1000, 10000),
                  labels = trans_format("log10", scales::math_format(10^.x))) +
    coord_cartesian(xlim = x_limits, ylim = y_limits) +
    theme_bw() +
    theme(text = element_text(family = theme_family))

  # Repeated event classes are valid; only distinct labels must be unambiguous.
  if (!is.null(selected_rows)) {
    if (!"classes" %in% colnames(dat)) stop("dat must contain classes to select clusters.")
    event_ids <- .mbr_event_ids(dat[, "classes"])
    selected_index <- .mbr_cluster_match(selected_rows, unique(event_ids), "Event classes")
    selected_ids <- unique(event_ids)[selected_index]
    p <- p +
      geom_point(data = dat[event_ids %in% selected_ids, , drop = FALSE],
                 aes(x = .data[[x]], y = .data[[y]]),
                 color = point_color,
                 alpha = point_alpha,
                 size = point_size)
  }

  # Add bin labels if bins is provided
  if (!is.null(bins) && !is.null(selected_rows)) {
    # Align coordinates to requested labels, independent of bin row order/case.
    bin_index <- .mbr_cluster_match(selected_rows, rownames(bins), "Bin rows")
    bins <- bins[bin_index, , drop = FALSE]
    p <- p +
      annotate('text',
               x = bins[, x],
               y = bins[, y],
               label = selected_rows,
               color = label_color,
               size = label_size)
  }

  return(p)
}

#' Arrange multiple flow density plots
#'
#' @description
#' Calls MBR_flow_plot for each channel pair and arranges the resulting plots with ggpubr.
#'
#' @param dat A data processed with MBR_process.
#' @param bins A data processed with MBR_prepare
#' @param plot_params A list of channel.
#' @param selected_rows clusters of interest
#' @param ncol Integer, number of columns for arranging the plots in the grid. Defaults to 2.
#' @param nrow Integer, number of rows for arranging the plots in the grid. Defaults to 2.
#' @param common_legend Use a shared legend for all plots (default: `TRUE`).
#' @param ... Additional arguments forwarded to every MBR_flow_plot call.
#'
#'
#' @export
MBR_plot <- function(dat,
                                  bins = NULL,
                                  plot_params,
                                  selected_rows = NULL,
                                  ncol = 2,
                                  nrow = 2,
                                  common_legend = TRUE,
                                  ...) {

  # Create a list to store plots
  plots <- list()

  # Create each plot
  for (i in seq_along(plot_params)) {
    params <- plot_params[[i]]
    if (!("x" %in% names(params)) || !("y" %in% names(params))) {
      stop("Each plot_params item must contain 'x' and 'y' parameters")
    }

    plots[[i]] <- MBR_flow_plot(dat = dat,
                                   bins = bins,
                                   x = params$x,
                                   y = params$y,
                                   selected_rows = selected_rows,
                                   ...)
  }

  # Arrange plots in a grid
  combined_plot <- ggarrange(plotlist = plots,
                             ncol = ncol,
                             nrow = nrow,
                             common.legend = common_legend)

  return(combined_plot)
}

#' Export selected cluster events to a new FCS file
#'
#' @description
#' Filters processed events by cluster and rebuilds a flowFrame using original FCS parameter metadata.
#'
#' @param fcs_files A `flowSet` object read by `MBR_read()`, containing the original FCS data.
#' @param dat A data frame returned by `MBR_process()`, containing transformed and labeled events.
#' @param selected_rows Cluster labels (e.g., "V170", "V214") indicating the events to keep.
#' @param file_index An integer index indicating which FCS file to extract from `fcs_files`. Defaults to 1.
#' @param rawdata_path Parent directory in which FilteredFCS will be created.
#'
#'
#' @export
MBR_save <- function(fcs_files, dat, selected_rows, file_index = 1, rawdata_path) {
  # 1. Get the FlowFrame properly from flowSet using [[ ]]
  target_fcs_filename <- sampleNames(fcs_files)[file_index]
  original_flowframe <- fcs_files[[file_index]]

  # 2. Validate requested clusters before export, then preserve event order.
  event_ids <- .mbr_event_ids(dat[, "classes"])
  selected_index <- .mbr_cluster_match(selected_rows, unique(event_ids), "Event classes")
  selected_events_data <- dat[event_ids %in% unique(event_ids)[selected_index], , drop = FALSE]
  selected_events_data_for_fcs <- selected_events_data %>% dplyr::select(-classes)

  # 3. Extract and prepare parameter info
  original_params_df <- Biobase::pData(flowCore::parameters(original_flowframe))
  filtered_params_df <- original_params_df[match(colnames(selected_events_data_for_fcs), original_params_df$name), ]

  required_cols <- c("name", "desc", "range", "minRange", "maxRange")
  missing_cols <- setdiff(required_cols, colnames(filtered_params_df))
  if (length(missing_cols) > 0) {
    for (col in missing_cols) {
      filtered_params_df[[col]] <- NA
    }
    warning("Some required parameter columns are missing and have been filled with NA.")
  }

  filtered_params_df$name <- colnames(selected_events_data_for_fcs)
  rownames(filtered_params_df) <- colnames(selected_events_data_for_fcs)
  new_parameters_ADF <- Biobase::AnnotatedDataFrame(filtered_params_df)

  # 4. Rebuild description
  original_desc <- flowCore::keyword(original_flowframe)
  new_desc <- list()

  for (key in names(original_desc)) {
    if (!grepl("^\\$P[0-9]+", key)) {
      new_desc[[key]] <- original_desc[[key]]
    }
  }

  original_param_name_to_index <- setNames(rownames(original_params_df), original_params_df$name)

  for (i in 1:nrow(filtered_params_df)) {
    current_param_name <- filtered_params_df$name[i]
    original_p_index_str <- original_param_name_to_index[current_param_name]

    if (is.na(original_p_index_str)) {
      warning(paste("Original information for parameter", current_param_name, "was not found; setting defaults."))
      new_desc[[paste0("$P", i, "N")]] <- current_param_name
      new_desc[[paste0("$P", i, "S")]] <- current_param_name
      new_desc[[paste0("$P", i, "R")]] <- as.character(filtered_params_df$range[i])
      new_desc[[paste0("$P", i, "E")]] <- "0,0"
      new_desc[[paste0("$P", i, "F")]] <- "0"
      new_desc[[paste0("$P", i, "B")]] <- "16"
      next
    }

    for (key_orig in names(original_desc)) {
      if (grepl(paste0("^\\$P", original_p_index_str), key_orig)) {
        suffix <- sub(paste0("^\\$P", original_p_index_str), "", key_orig)
        new_key <- paste0("$P", i, suffix)
        new_desc[[new_key]] <- original_desc[[key_orig]]
      }
    }
  }

  new_desc[["$TOT"]] <- as.character(nrow(selected_events_data_for_fcs))
  new_desc[["$PAR"]] <- as.character(ncol(selected_events_data_for_fcs))

  # 5. Create new FlowFrame
  new_flowframe <- flowCore::flowFrame(
    exprs = as.matrix(selected_events_data_for_fcs),
    parameters = new_parameters_ADF,
    description = new_desc
  )

  # 6. Save to FilteredFCS directory
  output_dir <- file.path(rawdata_path, "FilteredFCS")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  output_filename <- file.path(output_dir, paste0(sub("\\.fcs$", "", target_fcs_filename), "_filtered.fcs"))
  write.FCS(new_flowframe, filename = output_filename)

  cat("Save complete: ", output_filename, "\n")
  invisible(output_filename)
}
