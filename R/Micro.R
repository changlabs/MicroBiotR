# Script guide: see R/README.md and docs/ for workflow, statistical interpretation,
# examples, side effects, and known implementation limitations.
# See docs/implementation-notes.md for behavior and remaining limitations.

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


#' @title Train a shared SOM and export sample cluster abundances
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
#' @details
#' Each named sample must contain name and a numeric data matrix with identical channels in identical order. The pooled training subset gives each sample an equal target contribution, capped by its available events. A SOM learns prototype vectors (codebook rows) representing similar multichannel events; mapping assigns each event to a best matching prototype. No explicit channel standardization is applied here. Channel units and transformations therefore affect the fit.
#'
#' The grid side is round(sqrt(m)); the actual node count is its square. Despite the n_hclust argument, this implementation does not perform hierarchical clustering. Node labels and count levels are fixed at 2025. The default m = 2000 produces a 45 by 45 grid, matching that mapping; other grid sizes can omit nodes or create unused levels. Samples with fewer than 200000 events are removed before sampling and training. If none remain, the function stops.
#'
#' Rare clusters are counted by feature column: a cluster is rare when its abundance is below 0.01% in every retained sample.
#'
#' Writes cluster_information.csv, SOM.csv, count_tables/count_table_SOM.save, and one FCS file per sample in MappedFCS. Assigns cohonen_information, raw_count_table, and count_table in the global environment. The codebook spelling is intentionally retained. Seeds are not set internally; call set.seed before use. Existing output files may be replaced.
#'
#' @return
#' No analytical object is returned; the last console call returns invisible NULL. Retrieve tables from the documented files or global variables.
#'
#' @examples
#' \dontrun{
#' # Large example: use real gated samples with at least 200000 events each.
#' # gated <- list(sample1 = list(name = "sample1.fcs", data = events1),
#' #               sample2 = list(name = "sample2.fcs", data = events2))
#' set.seed(42)
#' out <- tempfile("som-")
#' MBR_som(gated_fcs = gated, out_path = out)
#' stopifnot(file.exists(file.path(out, "SOM.csv")),
#'           all(abs(rowSums(count_table) - 1) < 1e-8))
#' }
#'
#' @export
MBR_som <- function(gated_fcs = NULL, m = 2e3, n_hclust = 300,
                    n_cells_sub = 3e5, out_path = './') {

  # ================================================================
  # Stage 1: Environment Setup and Validation
  # ================================================================
  cat("==========================================\n")
  cat("Starting Analysis Pipeline\n")
  cat("==========================================\n")

  # Create output directory
  if (!dir.exists(out_path)) {
    dir.create(out_path, recursive = TRUE)
    cat("✓ Created output directory:", out_path, "\n")
  }

  # ================================================================
  # Stage 5: Sample Quality Control
  # ================================================================
  cat("\nStep 4: Sample Quality Control\n")
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
    cat("⚠ Dropped", length(dropped_samples), "samples\n")
  }

  # Filter samples
  gated_fcs <- gated_fcs[!low_count_samples]
  cat("✓ Retained", length(gated_fcs), "/", original_sample_count, "samples\n")

  if (!length(gated_fcs)) stop("No samples meet the 200000-event minimum.")

  # ================================================================
  # Stage 2: Random Sampling
  # ================================================================
  cat("\nStep 1: Performing Random Sampling\n")
  cat("----------------------------------------\n")
  cat("Target sampling size:", format(n_cells_sub, scientific = FALSE), "cells\n")

  # Execute sampling
  sub_sample <- downsample(
    gated_fcs,
    n = n_cells_sub,
    samples = names(gated_fcs)
  )
  sub_sample <- list(name = "random sample", data = sub_sample)

  cat("✓ Sampling completed, obtained", nrow(sub_sample$data), "cells\n")

  # ================================================================
  # Stage 3: SOM Training
  # ================================================================
  cat("\nStep 2: Training Self-Organizing Map (SOM)\n")
  cat("----------------------------------------\n")
  cat("Calculating SOM grid size...\n")

  # Calculate SOM parameters
  cells_per_node <- nrow(sub_sample$data) / m
  cat("Expected cells per node:", round(cells_per_node, 2), "\n")

  # Train SOM
  sub_sample <- compute_som(sub_sample, n_cells = cells_per_node)

  cat("✓ SOM training completed\n")
  cat("SOM grid size:", dim(sub_sample$som$grid$pts)[1], "nodes\n")

  # ================================================================
  # Stage 4: Save SOM Information
  # ================================================================
  cat("\nStep 3: Saving SOM Codebook Information\n")
  cat("----------------------------------------\n")

  # Save SOM codebook
  cohonen_information <- as.data.frame(sub_sample[["som"]][["codes"]][[1]])
  assign("cohonen_information", cohonen_information, envir = .GlobalEnv)

  som_file <- file.path(out_path, "cluster_information.csv")
  write.csv(cohonen_information, som_file)

  cat("✓ SOM codebook saved to:", som_file, "\n")
  cat("Codebook dimensions:", dim(cohonen_information), "\n")

  # ================================================================
  # Stage 6: Mapping to SOM
  # ================================================================
  cat("\nStep 5: Mapping All Samples to Trained SOM\n")
  cat("----------------------------------------\n")

  n_subset_large <- 1e10  # Use all cells

  cat("Mapping", length(gated_fcs), "samples to SOM...\n")

  # Map to trained SOM (with progress bar)
  SOM_fcs <- pbapply::pblapply(gated_fcs, map_som,
                               trained = sub_sample$som,
                               n_subset = n_subset_large)

  cat("✓ SOM mapping completed\n")

  # ================================================================
  # Stage 7: Cluster Assignment
  # ================================================================
  cat("\nStep 6: Assigning Cluster Labels\n")
  cat("----------------------------------------\n")

  # Create cluster number matrix
  cluster_number <- as.matrix(as.list(1:2025))
  colnames(cluster_number) <- as.character(2025)

  cat("Total clusters:", ncol(cluster_number), "\n")
  cat("Assigning cluster labels...\n")

  # Assign clusters (with progress bar)
  SOM_fcs <- pbapply::pblapply(SOM_fcs, assign_clusters, clusters = cluster_number)

  cat("✓ Cluster assignment completed\n")

  # ================================================================
  # Stage 8: Counting and Statistics
  # ================================================================
  cat("\nStep 7: Computing Cluster Statistics\n")
  cat("----------------------------------------\n")

  cat("Calculating cell counts for each cluster...\n")

  # Count observations (with progress bar)
  SOM_fcs <- pbapply::pblapply(SOM_fcs, count_observations,
                               clusters = colnames(cluster_number))

  # Get count tables
  count_tables <- get_counts(SOM_fcs)

  cat("✓ Counting completed\n")

  # ================================================================
  # Stage 9: Save Count Results
  # ================================================================
  cat("\nStep 8: Saving Count Results\n")
  cat("----------------------------------------\n")

  # Save count tables
  tmp_path <- file.path(out_path, "count_tables")
  if (!dir.exists(tmp_path)) dir.create(tmp_path, recursive = TRUE)

  count_file <- file.path(tmp_path, "count_table_SOM.save")
  save(count_tables, file = count_file, compress = TRUE)

  cat("✓ Raw count table saved to:", count_file, "\n")

  # Process and normalize counts
  raw_count_table <- as.data.frame(t(count_tables[[1]]))
  count_table_2025 <- raw_count_table / rowSums(raw_count_table)

  # Check for rare clusters
  rare_clusters <- sum(apply(count_table_2025 < 1e-4, 2, all))
  if (rare_clusters > 0) {
    warning(paste0(rare_clusters, ' clusters are below 0.01% of cells in ALL samples!'))
    cat("Found", rare_clusters, "rare clusters present in all samples\n")
  } else {
  	cat("✓ No clusters found below 0.01% in all samples (0 rare clusters)\n")
  }

  # Save normalized count table
  som_file <- file.path(out_path, "SOM.csv")
  write.csv(count_table_2025, som_file)

  cat("✓ Normalized count table saved to:", som_file, "\n")

  # Save to global environment
  assign("raw_count_table", raw_count_table, envir = .GlobalEnv)
  assign("count_table", count_table_2025, envir = .GlobalEnv)

  # ================================================================
  # Stage 10: Export FCS Files
  # ================================================================
  cat("\nStep 9: Exporting Mapped FCS Files\n")
  cat("----------------------------------------\n")

  # Create output directory
  dir_name <- file.path(out_path, "MappedFCS")
  if (!dir.exists(dir_name)) {
    dir.create(dir_name, recursive = TRUE)
    cat("✓ Created FCS output directory:", dir_name, "\n")
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
  cat("\n✓ FCS file export completed\n")

  # ================================================================
  # Completion Summary
  # ================================================================
  cat("\n==========================================\n")
  cat("Analysis Pipeline Completed!\n")
  cat("==========================================\n")
  cat("Processed samples:", length(SOM_fcs), "\n")
  cat("Total clusters:", ncol(cluster_number), "\n")
  cat("Output directory:", out_path, "\n")
  cat("\nMain output files:\n")
  cat("- SOM codebook:", file.path(out_path, "cluster_information.csv"), "\n")
  cat("- Normalized counts:", file.path(out_path, "SOM.csv"), "\n")
  cat("- Raw counts:", file.path(tmp_path, "count_table_SOM.save"), "\n")
  cat("- Mapped FCS files:", dir_name, "\n")
  cat("==========================================\n")
}

#' Test each abundance feature between sample groups
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
#' @details
#' Numeric matrices and data frames are accepted. Sample row names are matched to metadata row names and reordered when needed; mismatched IDs stop the analysis. Without usable IDs, positional matching emits a warning.
#'
#' wilcox uses the unpaired, two-sided Wilcoxon rank-sum test for two groups. It compares rank distributions; interpreting it purely as a median comparison requires similarly shaped distributions. kruskal uses the Kruskal-Wallis rank test for two or more independent groups; a significant omnibus test does not identify the differing pairs. anova fits a one-way ANOVA testing equality of means, with independent errors and approximately normal residuals of similar variance. t.test (also ttest) uses the unpaired two-sided Welch test, allowing unequal variances, for two groups. None of these branches supports paired samples, covariate adjustment, or repeated measurements.
#'
#' BH and fdr both request Benjamini-Hochberg adjustment, controlling false discovery rate under its dependence assumptions. bonferroni controls family-wise error through a more conservative adjustment. none, the default, leaves p-values unadjusted. Selection uses p.adj < cutoff, strictly excluding equality. Feature-wise tests on proportions remain affected by compositional dependence.
#'
#' Wilcoxon tests use a normal approximation (exact = FALSE), including ties. Failed or non-finite feature tests emit a feature-specific warning and retain NA p-values; missing results are excluded from significance selection. Invalid test names stop before testing.
#'
#' @return
#' The significant feature data frame, returned invisibly by the final assign call; pvalue_data is also assigned globally.
#'
#' @references
#' \url{https://stat.ethz.ch/R-manual/R-devel/library/stats/html/wilcox.test.html}
#' \url{https://stat.ethz.ch/R-manual/R-devel/library/stats/html/t.test.html}
#' \url{https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html}
#'
#' @examples
#' set.seed(42)
#' data <- as.data.frame(matrix(runif(120, 0.01, 1), nrow = 20))
#' names(data) <- paste0("v", seq_len(ncol(data)))
#' data <- data / rowSums(data)
#' meta <- data.frame(Group = rep(c("A", "B"), each = 10))
#' out <- tempfile("MicroBiotR-")
#' dir.create(out)
#' MBR_stat(data, meta, "Group", correction = "BH", out_path = out)
#' stopifnot(nrow(pvalue_data) == ncol(data),
#'           all(pvalue_data$p.adj >= 0 & pvalue_data$p.adj <= 1),
#'           file.exists(file.path(out, "pvalue.csv")))
#'
#' @export
MBR_stat <- function(data = NULL, meta_data = NULL, group_col = NULL,
                     test_type = 'wilcox', cutoff = 0.05,
                     correction = 'none', out_path = './') {
  inputs <- .mbr_inputs(data, meta_data, group_col)
  data <- inputs$data
  meta_data <- inputs$meta_data


  cat(paste0(test_type, " test\n"))

  # Input validation
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
  # Warnings retain valid test results; failed tests remain explicitly missing.
  original_pvalues <- vapply(names(data), function(feature) {
    t <- data[[feature]]
    tryCatch({
      value <- switch(test_type,
        wilcox = stats::wilcox.test(t ~ group_vec, exact = FALSE)$p.value,
        kruskal = stats::kruskal.test(t ~ group_vec)$p.value,
        anova = stats::anova(stats::aov(t ~ group_vec))$`Pr(>F)`[1],
        t.test = stats::t.test(t ~ group_vec)$p.value,
        ttest = stats::t.test(t ~ group_vec)$p.value)
      if (!is.finite(value)) stop("Test returned a non-finite p-value.")
      value
    }, error = function(e) {
      warning(sprintf("Feature '%s': %s", feature, conditionMessage(e)), call. = FALSE)
      NA_real_
    })
  }, numeric(1))

  # Create a data frame for p-values
  pvalue_df <- data.frame(p.value = original_pvalues)

  # Apply correction if specified
  if (correction == 'fdr') {
    pvalue_df$p.adj <- p.adjust(original_pvalues, method = 'fdr')
  } else if (correction == 'bonferroni') {
    pvalue_df$p.adj <- p.adjust(original_pvalues, method = 'bonferroni')
  } else if (correction == 'BH') {
    pvalue_df$p.adj <- p.adjust(original_pvalues, method = 'BH')
  } else if (correction == 'none') {
    pvalue_df$p.adj <- original_pvalues
  } else {
    stop('Error: correction should be one of "none", "fdr", "bonferroni", "BH".')
  }


  # Filter significant data
  significant_data <- data[, !is.na(pvalue_df$p.adj) & pvalue_df$p.adj < cutoff, drop = FALSE]

  # Save CSV files
  write.csv(pvalue_df, file.path(out_path, "pvalue.csv"), row.names = TRUE)
  write.csv(significant_data, file.path(out_path, "significant.csv"), row.names = TRUE)

  cat("\tDONE.\n")

  # Assign to global environment
  assign("pvalue_data", pvalue_df, envir = .GlobalEnv)
  assign("significant_data", significant_data, envir = .GlobalEnv)
}

#' Draw group means in a circular heatmap
#'
#' @description
#' Displays feature means by group and a separate track for the difference between the first two group means.
#'
#' @param data A numeric data frame with samples in rows and features in columns.
#' @param meta_data A data frame containing sample metadata with grouping information.
#' @param group_col Character string specifying the column name in `meta_data` for group labels.
#' @param out_path Directory path to save the output PDF file (default: current directory).
#' @param width Width of the output PDF (default: 8).
#' @param height Height of the output PDF (default: 8).
#' @param point_colors Vector of two colors for points representing positive and negative mean differences (default: c("red", "blue")).
#' @param cell_colors Vector of three colors for the heatmap gradient (default: c("blue", "white", "red")).
#'
#' @details
#' Use a numeric data frame and row-aligned metadata. tapply orders the grouping levels; mean_table columns, rather than order of appearance in metadata, determine the subtraction first group minus second group. At least two groups are needed and the difference track is intended for two groups. Means are calculated without na.rm, so missing values propagate. The color midpoint is the midpoint of the overall minimum and maximum, not the overall mean. Constant means or differences can make color breaks or track limits degenerate.
#'
#' This is a descriptive plot, with no hypothesis test or uncertainty interval. Group means can hide sample heterogeneity. The function draws on the active graphics device, resets circlize state, then redraws into circle.pdf. out_path must already exist.
#'
#' @return
#' Invisible NULL from the final console call; the plot is drawn and saved.
#'
#' @examples
#' set.seed(42)
#' data <- as.data.frame(matrix(runif(120, 0.01, 1), nrow = 20))
#' names(data) <- paste0("v", seq_len(ncol(data)))
#' data <- data / rowSums(data)
#' meta <- data.frame(Group = rep(c("A", "B"), each = 10))
#' out <- tempfile("MicroBiotR-")
#' dir.create(out)
#' MBR_circle(data, meta, "Group", out_path = out)
#' stopifnot(file.exists(file.path(out, "circle.pdf")))
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
  rownames(mean_table) <- gsub('V', 'Bin', rownames(mean_table))
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

#' Plot one cluster abundance and an existing p-value
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
#' @param width Width of the saved PDF (default: 4).
#' @param height Height of the saved PDF (default: 4).
#'
#' @details
#' Numeric matrices and data frames are accepted and converted to a data frame for named column extraction. cluster = 1 selects the lowercase column v1 and the row v1 of pvalue_data. Metadata rows must match data order. p specifies an existing p-value column, normally p.value or p.adj. No test is performed by this function; annotation validity depends on the provenance and correction of the supplied table.
#'
#' A violin is a smoothed density estimate, and the overlaid box plot summarizes the median and interquartile range. The abundance label assumes data already contains proportions. Supply enough colors for all groups. The annotation is positioned at x = 1.5, appropriate for a two-group comparison. Writes violin_cluster_<cluster>_<p>.pdf to an existing directory and prints the plot.
#'
#' @return
#' Invisible NULL; the plot is printed and saved, rather than returned.
#'
#' @examples
#' set.seed(42)
#' data <- as.data.frame(matrix(runif(120, 0.01, 1), nrow = 20))
#' names(data) <- paste0("v", seq_len(ncol(data)))
#' data <- data / rowSums(data)
#' meta <- data.frame(Group = rep(c("A", "B"), each = 10))
#' out <- tempfile("MicroBiotR-")
#' dir.create(out)
#' MBR_stat(data, meta, "Group", correction = "BH", out_path = out)
#' MBR_violin(data, meta, "Group", pvalue_data, p = "p.adj",
#'            cluster = 1, out_path = out)
#' stopifnot(file.exists(file.path(out, "violin_cluster_1_p.adj.pdf")))
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



  cluster_name <- paste0("v", cluster)

  cluster_values <- data[[cluster_name]]

  p_val_numeric <- as.numeric(pvalue_data[cluster_name, p])
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

#' Ordinate Bray-Curtis dissimilarity and test group differences
#'
#' @description
#' Calculates Bray-Curtis dissimilarities, a principal coordinates display, and a PERMANOVA group test with coordinate-wise comparison plots.
#'
#' @param data A numeric data frame or matrix of features (samples as rows).
#' @param out_path Output directory path to save plots and results (default: "./").
#' @param test Statistical test for group differences: one of "wilcox", "ttest", "kruskal", or "anova" (default: "wilcox").
#' @param meta_data A data frame of metadata corresponding to samples.
#' @param group_name The column name in `meta_data` to define groups for comparison.
#' @param colors Vector of colors for groups in plots (default: c('#E41A1C', '#377EB8')).
#' @param width Width of output plot in inches (default: 5).
#' @param height Height of output plot in inches (default: 5).
#'
#' @details
#' Rows are samples and columns are nonnegative features. Metadata must be row-aligned. For samples i and j, Bray-Curtis is sum(abs(x_i - x_j)) / sum(x_i + x_j); it emphasizes abundance differences and ignores joint absences. Empty samples and missing values need to be handled before calling. Raw counts retain library-size effects; this function does not normalize input.
#'
#' Classical multidimensional scaling via cmdscale requests three axes and displays the first two. Use sufficiently many distinct samples to obtain three positive axes (at least four samples). Bray-Curtis need not be Euclidean, so negative eigenvalues can occur. The displayed percentages divide each eigenvalue by the sum of ALL eigenvalues; they are not necessarily conventional fractions of positive inertia. No correction for negative eigenvalues is requested.
#'
#' adonis2 tests group-associated variation in the full distance matrix using dependency-default unrestricted permutations. Its R-squared summarizes the fraction of sums of squares attributed to group. Group dispersion differences can also affect the result; examine dispersion separately. Exchangeability is required, and this interface has no strata or covariate arguments. Set a seed for reproducible permutations.
#'
#' The test argument changes only comparisons of plotted coordinate scores, not PERMANOVA. wilcox and ttest compare group pairs; kruskal is called through pairwise plot comparisons; anova gives an omnibus axis comparison. Axis tests do not replace a multivariate test. The wrapper specifies no multiple-testing correction for these annotations. Ellipses are visualization summaries, not a PERMANOVA confidence region. Writes .group.txt, adonis.txt, and pcoa_<test>.pdf to an existing out_path.
#'
#' @return
#' The combined patchwork plot, returned visibly.
#'
#' @references
#' \url{https://vegandevs.github.io/vegan/reference/vegdist.html}
#' \url{https://vegandevs.github.io/vegan/reference/adonis.html}
#' \url{https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cmdscale.html}
#'
#' @examples
#' set.seed(42)
#' data <- as.data.frame(matrix(runif(120, 0.01, 1), nrow = 20))
#' names(data) <- paste0("v", seq_len(ncol(data)))
#' data <- data / rowSums(data)
#' meta <- data.frame(Group = rep(c("A", "B"), each = 10))
#' out <- tempfile("MicroBiotR-")
#' dir.create(out)
#' p <- MBR_beta(data, out_path = out, meta_data = meta, group_name = "Group")
#' stopifnot(inherits(p, "patchwork"), file.exists(file.path(out, "adonis.txt")))
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
                 color = "group", shape = "group", linewidth = 3,
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

#' Select features with random forest recursive elimination
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
#' @param width Width of the output PDF plot in inches (default: 5).
#' @param height Height of the output PDF plot in inches (default: 5).
#'
#' @details
#' The wrapper sets the random seed to 2025. Row-aligned numeric features and a metadata grouping column are required. caret::rfe with rfFuncs evaluates requested subset sizes 1:rfe_size and, through caret defaults, the full predictor set using nfolds_cv-fold cross-validation; use rfe_size no greater than the feature count and enough observations per class for all folds. The selected subset is results$optVariables. top_n_features controls how many selected variables are displayed, not the size searched.
#'
#' Random forests combine trees fitted to bootstrap samples and random predictor subsets. RFE ranks predictors and compares candidate subsets using resampling. Importance can be unstable when predictors are correlated. Effect-size plots use log10(value + 1), despite an axis label saying log10(abundance). Cohen's d is a standardized difference in group means; its sign depends on group ordering and ref_group. The displayed 95 percent intervals and effects are calculated after selecting features on these same data, so they are exploratory rather than selection-adjusted inference. Use two groups for this workflow and keep feature selection inside an outer resampling loop when evaluating predictive performance.
#'
#' Writes .group.txt, figure1.txt (selected features plus group), and feature_exploration.pdf. Assigns MBR_selected_features using superassignment, normally into the global environment. Requires an existing output directory. Missing importance rows are removed with tidyr::drop_na. Three console plots are shown with readline prompts, so this workflow is interactive. The fitted RFE object is not returned.
#'
#' @return
#' Invisible NULL; selected feature data are available as MBR_selected_features and files.
#'
#' @references
#' \url{https://topepo.github.io/caret/recursive-feature-elimination.html}
#' \url{https://rpkgs.datanovia.com/rstatix/reference/cohens_d.html}
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' data <- as.data.frame(matrix(runif(120, 0.01, 1), nrow = 20))
#' names(data) <- paste0("v", seq_len(ncol(data)))
#' data <- data / rowSums(data)
#' meta <- data.frame(Group = rep(c("A", "B"), each = 10))
#' out <- tempfile("MicroBiotR-")
#' dir.create(out)
#' MBR_fs(data, out_path = out, meta_data = meta, group_name = "Group",
#'        nfolds_cv = 2, rfe_size = 3, top_n_features = 3, ref_group = "A")
#' stopifnot(nrow(MBR_selected_features) == nrow(data),
#'           all(names(MBR_selected_features) %in% names(data)))
#' }
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

  # Save selected features ONLY (no group) from input data to global environment
  MBR_selected_features <<- data[, best_selection, drop = FALSE]
  rownames(MBR_selected_features) <<- rownames(data)

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
}

#' Show channel values for selected SOM prototypes
#'
#' @description
#' Selects codebook rows using feature column names and renders a configurable pheatmap.
#'
#' @param data Data frame or matrix whose column names correspond to SOM cluster IDs (e.g., "V1", "V2", ...).
#' @param cohonen_information Data frame containing information per SOM cluster. Row names should match cluster IDs.
#' @param out_path Directory path where the heatmap PDF will be saved (default './').
#' @param color Color palette for the heatmap (default is reversed RdYlBu palette from RColorBrewer).
#' @param width Width of the output PDF in inches (default 10).
#' @param height Height of the output PDF in inches (default 10).
#' @param scale Scaling method for pheatmap: 'row', 'column', or 'none' (default 'row').
#' @param cluster_rows Logical, whether to cluster rows (default FALSE).
#' @param cluster_cols Logical, whether to cluster columns (default FALSE).
#' @param display_numbers Logical, whether to display values on heatmap cells (default TRUE).
#' @param boarder_color Color of the border around heatmap cells (default 'grey60').
#' @param legend Logical, whether to display the legend (default TRUE).
#'
#' @details
#' Only the column names of data select clusters; its abundance values are not plotted. Both v1 and V1 are converted to V1 for lookup, so cohonen_information must have matching UPPERCASE V-prefixed row names. MBR_som's codebook normally has numeric row names; adapt a copy before passing it here. Missing matches yield missing rows.
#'
#' The plotted matrix contains prototype channel values. Row scaling standardizes channels within each selected node; column scaling standardizes nodes within each channel; none preserves input units. These different choices change the interpretation of color. Constant rows or columns may fail when scaled. Optional clustering uses pheatmap defaults; it does not alter the package's original SOM assignments. A numeric conversion df_mat is computed but the original subset df is actually plotted. The argument boarder_color is intentionally spelled as in the implementation. Writes heatmap.pdf and draws the plot; out_path must exist.
#'
#' @return
#' Invisible NULL; the pheatmap is drawn and saved.
#'
#' @examples
#' set.seed(42)
#' data <- as.data.frame(matrix(runif(120, 0.01, 1), nrow = 20))
#' names(data) <- paste0("v", seq_len(ncol(data)))
#' data <- data / rowSums(data)
#' meta <- data.frame(Group = rep(c("A", "B"), each = 10))
#' out <- tempfile("MicroBiotR-")
#' dir.create(out)
#' codebook <- as.data.frame(matrix(seq_len(18), nrow = 6))
#' rownames(codebook) <- paste0("V", 1:6)
#' MBR_heatmap(data, codebook, out_path = out, scale = "none")
#' stopifnot(file.exists(file.path(out, "heatmap.pdf")))
#'
#' @export
MBR_heatmap <- function(data = NULL, cohonen_information = NULL,
                      out_path = './',
                      color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                      width = 10, height = 10,
                      scale = 'row', cluster_rows = F, cluster_cols = F,
                      display_numbers = T, boarder_color = 'grey60',
                      legend = T) {
  cluster_ids <- paste0('V', gsub('^[vV]', '', colnames(data)))
  rownames(cohonen_information) <- as.character(rownames(cohonen_information))
  df <- cohonen_information[cluster_ids, , drop = FALSE]
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
#' @param width Width of output PDF in inches (default 8).
#' @param height Height of output PDF in inches (default 8).
#'
#' @details
#' Samples must occupy matching rows of numeric data and metadata. clinical_cols and demographic_cols select metadata variable blocks; spec_select_names$A and $B supply their display names. NULL labels omit a block. Only appropriate numeric variables should enter distance calculations; arbitrary numeric encoding of categories gives arbitrary distances.
#'
#' A Mantel statistic correlates entries of two sample-distance matrices and uses permutations for significance. Distances and permutation settings are delegated to linkET::mantel_test defaults, not specified by this wrapper; consult the installed linkET help to identify them. Unequal units across metadata variables can dominate distance, so consider scaling a copy of metadata before use. Permutations assume exchangeable samples, and no blocked design or confounder adjustment is exposed. This is an association analysis, not evidence of causality.
#'
#' correlate(data) supplies the heatmap correlations (Pearson by the dependency default). Links bin Mantel r at 0.2 and 0.4 and p at 0.01 and 0.05; cut uses right-closed intervals, so boundary values enter the lower interval even though labels use less-than signs. Negative r values share the first category. The code specifies no multiple-testing adjustment. Writes mantel.pdf to an existing directory, prints the plot, and does not return the Mantel table.
#'
#' @return
#' Invisible NULL; the association figure is saved and printed.
#'
#' @references
#' \url{https://github.com/Hy4m/linkET/blob/master/R/mantel-test.R}
#'
#' @examples
#' set.seed(42)
#' data <- as.data.frame(matrix(runif(120, 0.01, 1), nrow = 20))
#' names(data) <- paste0("v", seq_len(ncol(data)))
#' data <- data / rowSums(data)
#' meta <- data.frame(Group = rep(c("A", "B"), each = 10))
#' out <- tempfile("MicroBiotR-")
#' dir.create(out)
#' # Include a strongly associated feature pair so the significance layer is exercised.
#' data$v2 <- data$v1 + seq_len(nrow(data)) * 1e-6
#' data <- data / rowSums(data)
#' clinical <- data.frame(age = seq_len(20), marker = runif(20), BMI = rnorm(20, 25, 2))
#' MBR_mantel(data, clinical, clinical_cols = c("age", "marker"),
#'            demographic_cols = "BMI", out_path = out)
#' stopifnot(file.exists(file.path(out, "mantel.pdf")))
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
    spec_select[[ spec_select_names$A ]] <- clinical_cols
  }
  if (!is.null(spec_select_names$B)) {
    spec_select[[ spec_select_names$B ]] <- demographic_cols
  }

  mantel <- suppressMessages(
    mantel_test(meta_data, data, spec_select = spec_select) %>%
    mutate(
      rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
               labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
      pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
               labels = c("< 0.01", "0.01 - 0.05", ">= 0.05"))
    )
  )

    p <- qcorrplot(correlate(data), type = "lower", diag = FALSE) +
    geom_square() +
    geom_mark(sep='\n', size = 1.8, sig_level = c(0.05, 0.01, 0.001),
              sig_thres = 0.05, color="white")+
    geom_couple(aes(colour = pd, size = rd),
                data = mantel,
                curvature = nice_curvature()) +
    scale_fill_gradientn(colours = colors) +
    scale_size_manual(values = c(0.5, 1, 2)) +
    scale_colour_manual(values = color_pal(3)) +
    guides(size = guide_legend(title = "Mantel's r",
                               override.aes = list(colour = "grey35"),
                               order = 2),
           colour = guide_legend(title = "Mantel's p",
                                 override.aes = list(size = 3),
                                 order = 1),
           fill = guide_colorbar(title = "Pearson's r", order = 3))

  #pdf
  pdf(paste0(out_path, '/mantel.pdf'), width = width, height = height)
  print(p)
  dev.off()

  print(p)
  cat("\tDONE.\n")
}

#' Evaluate random forest classification with ROC
#'
#' @description
#' Fits a caret random forest classifier and saves a diagnostic plot based on held-out resampling predictions, or an optional independent test set.
#'
#' @param data A data frame or matrix with features in columns and samples in rows.
#' @param meta_data A data frame containing metadata for the samples.
#' @param group_name The column name in `meta_data` indicating the grouping variable (e.g., disease vs control).
#' @param out_path Output directory for saving the ROC/confusion plot and group assignment file (default: "./").
#' @param reference_level The reference level of the group used for binary classification (default: "A").
#' @param width Width of the output ROC/confusion plot PDF (default: 5).
#' @param height Height of the output ROC/confusion plot PDF (default: 5).
#' @param method Resampling method for training control in caret (default: "repeatedcv").
#' @param number Number of folds or resampling iterations (default: 5).
#' @param repeats Number of repeats for repeated cross-validation (default: 5).
#' @param test_data Optional numeric feature table for an independent test set, with the same feature names as data.
#' @param test_meta_data Metadata for test_data, with sample IDs and group_name. Supply both test arguments together.
#'
#' @details
#' Numeric predictors and row-aligned metadata with exactly two classes are required. Class labels must be valid R names for caret probability columns, and reference_level must identify a present class. Resampling settings are passed to caret::trainControl; caret tunes a random forest using ROC as its selection metric. The wrapper requires the randomForest backend. number and repeats configure the selected method; repeats is relevant to repeatedcv.
#'
#' By default, plots use held-out caret predictions for the selected tuning parameters, averaged across repeats once per sample. These are cross-validation diagnostics, not independent validation after tuning or upstream feature selection. For independent evaluation, supply test_data and test_meta_data from samples unused in feature selection, preprocessing estimation and tuning. The fitted model, predictions and confusion metrics are returned invisibly.
#'
#' Writes ml_group.txt and the PDF in out_path; creates the directory recursively if absent. Existing files can be overwritten.
#'
#' reference_level is the positive class for both ROC probabilities and confusion metrics. ROC direction is fixed so larger probabilities indicate the positive class. No global warning options are modified. Resampling confidence intervals are descriptive and omit tuning and feature-selection uncertainty.
#'
#' @return
#' An invisible list containing the fitted model, sample-level diagnostic predictions and confusion metrics; MBR_ml also returns the ROC object and plot.
#'
#' @references
#' \url{https://search.r-project.org/CRAN/refmans/randomForest/html/randomForest.html}
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' data <- as.data.frame(matrix(runif(120, 0.01, 1), nrow = 20))
#' names(data) <- paste0("v", seq_len(ncol(data)))
#' data <- data / rowSums(data)
#' meta <- data.frame(Group = rep(c("A", "B"), each = 10))
#' out <- tempfile("MicroBiotR-")
#' dir.create(out)
#' MBR_ml(data, meta, group_name = "Group", out_path = out,
#'        reference_level = "A", number = 2, repeats = 1)
#' stopifnot(file.exists(file.path(out, "roc_confusion.pdf")))
#' }
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
  plot_file <- file.path(out_path, 'roc_confusion.pdf')
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
#' Fits a caret random forest classifier and saves a diagnostic plot based on held-out resampling predictions, or an optional independent test set.
#'
#' @param data A data frame or matrix with features as columns and samples as rows.
#' @param meta_data A data frame containing metadata for the samples.
#' @param group_name The column name in `meta_data` that indicates the grouping variable (e.g., case vs control).
#' @param out_path Directory path to save the confusion matrix plot and group assignment file (default: "./").
#' @param reference_level The reference level of the group used for binary classification (default: "A").
#' @param colors A vector of two colors used for the heatmap gradient (default: c('#00BFC4', '#F8766D')).
#' @param width Width of the output PDF plot in inches (default: 5).
#' @param height Height of the output PDF plot in inches (default: 5).
#' @param method Resampling method for training control in caret (default: "repeatedcv").
#' @param number Number of folds or resampling iterations (default: 5).
#' @param repeats Number of repeats for repeated cross-validation (default: 5).
#' @param test_data Optional numeric feature table for an independent test set, with the same feature names as data.
#' @param test_meta_data Metadata for test_data, with sample IDs and group_name. Supply both test arguments together.
#'
#' @details
#' Numeric predictors and row-aligned metadata with exactly two classes are required. Class labels must be valid R names for caret probability columns, and reference_level must identify a present class. Resampling settings are passed to caret::trainControl; caret tunes a random forest using ROC as its selection metric. The wrapper requires the randomForest backend. number and repeats configure the selected method; repeats is relevant to repeatedcv.
#'
#' By default, plots use held-out caret predictions for the selected tuning parameters, averaged across repeats once per sample. These are cross-validation diagnostics, not independent validation after tuning or upstream feature selection. For independent evaluation, supply test_data and test_meta_data from samples unused in feature selection, preprocessing estimation and tuning. The fitted model, predictions and confusion metrics are returned invisibly.
#'
#' Writes ml_group.txt and the PDF in out_path; creates the directory recursively if absent. Existing files can be overwritten.
#'
#' Unlike MBR_ml, this function does not set a random seed; call set.seed before use. The heatmap labels are event counts (samples here), while fill is natural log(Freq + 1), making small cells visible alongside large ones. It shows no uncertainty intervals or ROC curve.
#'
#' @return
#' An invisible list containing the fitted model, sample-level diagnostic predictions and confusion metrics; MBR_ml also returns the ROC object and plot.
#'
#' @references
#' \url{https://search.r-project.org/CRAN/refmans/randomForest/html/randomForest.html}
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' data <- as.data.frame(matrix(runif(120, 0.01, 1), nrow = 20))
#' names(data) <- paste0("v", seq_len(ncol(data)))
#' data <- data / rowSums(data)
#' meta <- data.frame(Group = rep(c("A", "B"), each = 10))
#' out <- tempfile("MicroBiotR-")
#' dir.create(out)
#' MBR_conf(data, meta, group_name = "Group", out_path = out,
#'        reference_level = "A", number = 2, repeats = 1)
#' stopifnot(file.exists(file.path(out, "confusion_matrix.pdf")))
#' }
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

  # Save
  plot_file <- file.path(out_path, 'confusion_matrix.pdf')
  pdf(file = plot_file, width = width, height = height)
  print(g2_out)
  dev.off()
  cat("Confusion matrix plot saved to:", plot_file, "\n")

  print(g2_out)
  cat("\tDONE.\n")
  invisible(list(model = rf, predictions = diagnostic, confusion = cm))
}


#' Recluster numeric features using Ward linkage
#'
#' @description
#' Standardizes numeric columns, performs Euclidean hierarchical clustering, and appends cluster membership to a copy of the input.
#'
#' @param data A data frame containing numeric variables for clustering.
#' @param num_clusters Integer specifying the number of clusters to cut the dendrogram into.
#'
#' @details
#' Non-numeric columns are preserved but excluded from distance calculations. All numeric columns, including numeric identifiers or previous labels, enter clustering; remove unwanted numeric fields from a copy beforehand. scale centers each feature and divides by its sample standard deviation, giving differently measured channels comparable weight. Constant columns, missing or infinite values cause invalid distances.
#'
#' Ward.D2 linkage merges groups using a criterion related to the increase in within-cluster sum of squares; with Euclidean distances it favors compact clusters. cutree chooses the requested number of groups, which must be between 1 and the number of rows. Labels are arbitrary identifiers, not ordered biological states. This descriptive partition has no p-value or estimate of the optimal number of clusters. A pre-existing Cluster column is overwritten in the returned copy. Assigns reclustered_information globally.
#'
#' @return
#' The input data frame with an added or replaced Cluster column, returned invisibly and assigned to reclustered_information.
#'
#' @references
#' \url{https://stat.ethz.ch/R-manual/R-devel/library/stats/html/hclust.html}
#'
#' @examples
#' x <- data.frame(channel1 = c(1, 2, 8, 9), channel2 = c(2, 1, 9, 8),
#'                 label = letters[1:4])
#' y <- MBR_reclustering(x, num_clusters = 2)
#' stopifnot(nrow(y) == nrow(x), length(unique(y$Cluster)) == 2,
#'           identical(y$label, x$label))
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

  assign("reclustered_information", reclustered_information, envir = .GlobalEnv)

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
#' @details
#' Lists non-recursive filenames matching fcs$ (case-sensitive, without requiring a dot before fcs). Uppercase .FCS files are not selected. Files are read with alter.names = FALSE and channels matching an asterisk, Bits, or Drop excluded. The wrapper leaves transformation and other reader settings at flowCore defaults; it does not reproduce Gating's explicit transformation = FALSE. It stops for a missing directory or no matching files. Reading is not gating or compensation; verify instrument channels and reader settings before analysis.
#'
#' @return
#' A flowCore flowSet with one flowFrame per selected file.
#'
#' @examples
#' \dontrun{
#' fcs <- MBR_read("path/to/fcs_directory")
#' stopifnot(inherits(fcs, "flowSet"), length(fcs) > 0)
#' }
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
#' @param fcs_files A flowCore flowSet, normally returned by MBR_read.
#' @param file_index One-based sample index.
#' @param transformation Vectorized numeric transformation function; defaults to 10^((4*x)/65000).
#' @param column_mapping Optional named character vector mapping original channel names to new names.
#'
#' @details
#' The default transformation is 10^((4*x)/65000), an exponential rescaling tied to the original instrument range. It is not a universal cytometry compensation, logarithm, or arcsinh transformation. Choose a transformation suitable for your acquisition scale. Every column except a column named exactly classes is transformed; cluster labels are retained as labels. The selected sample is found through sampleNames and frames.
#'
#' column_mapping is a named character vector whose names are old channel names and values are replacement names. Unmatched names are ignored. Renaming happens after transformation. The function does not gate events or normalize sample abundance. Preserve original channel names when planning MBR_save, whose parameter matching uses those names.
#'
#' @return
#' A data frame of transformed events; rows retain event order and classes is preserved when present.
#'
#' @examples
#' \donttest{
#' ff <- flowCore::flowFrame(matrix(c(10, 20, 30, 40, 1, 2), nrow = 2,
#'                                  dimnames = list(NULL, c("FSC", "SSC", "classes"))))
#' fs <- flowCore::flowSet(list(sample1 = ff))
#' x <- MBR_process(fs, transformation = identity, column_mapping = c(FSC = "scatter"))
#' stopifnot("scatter" %in% names(x), identical(as.numeric(x$classes), c(1, 2)))
#' }
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

#' Select and transform prototype rows for plotting
#'
#' @description
#' Filters a prototype table by row name, optionally renames channels, and transforms all remaining columns.
#'
#' @param bins Numeric prototype matrix or data frame with cluster row names.
#' @param selected_rows Character vector of row names to retain.
#' @param transformation Vectorized numeric transformation applied to every column.
#' @param column_mapping Optional named character vector mapping old channel names to new names.
#'
#' @details
#' selected_rows is matched against existing row names using membership, retaining the INPUT order rather than the order of selected_rows. Missing requested rows are silently omitted. Unlike MBR_process, every column is transformed, including any label column supplied accidentally. Keep only appropriate numeric channels. The default exponential transformation assumes the original 65000 instrument scale and should match the transformation used for plotted events. Renaming occurs before transformation; a named character vector maps old names to new names. For plot labels, explicitly align selected_rows to the resulting row order.
#'
#' @return
#' A data frame containing the selected, renamed and transformed rows.
#'
#' @examples
#' bins <- data.frame(FSC = c(10, 20, 30), SSC = c(40, 50, 60),
#'                    row.names = c("V1", "V2", "V3"))
#' x <- MBR_prepare(bins, c("V3", "V1"), transformation = identity,
#'                  column_mapping = c(FSC = "scatter"))
#' stopifnot(identical(rownames(x), c("V1", "V3")), "scatter" %in% names(x))
#'
#' @export
MBR_prepare <- function(bins,
                         selected_rows,
                         transformation = function(x) 10^((4 * x) / 65000),
                         column_mapping = NULL) {

  # Select rows
  bins <- bins[rownames(bins) %in% selected_rows, , drop = FALSE]

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

#' Plot event density with optional cluster highlighting
#'
#' @description
#' Builds a two-channel hexagonal density plot with logarithmic axes, highlighted events, and optional prototype labels.
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
#' @param theme_family Font family for the plot theme; default Helvetica.
#'
#' @details
#' dat must contain positive numeric x and y channels for log axes. geom_hex counts events per hexagon and the fill scale is logarithmic. hex_bins sets the spatial resolution, so colors summarize counts per bin rather than normalized probability. The hexbin backend must be installed. Events outside scale limits are removed for plot calculations.
#'
#' Highlight matching uses paste0("V", dat$classes), requiring uppercase labels such as V1 in selected_rows. Other functions may use lowercase v abundance columns; these conventions are not interchangeable without conversion. bins must have x and y columns and exactly one correctly ordered row per label in selected_rows. No statistical test is performed. This returns a plot without saving or printing it.
#'
#' @return
#' A ggplot object that can be printed or saved with ggplot2::ggsave.
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' dat <- data.frame(FSC = runif(100, 1, 1000), SSC = runif(100, 1, 1000),
#'                   classes = rep(1:2, 50))
#' p <- MBR_flow_plot(dat, x = "FSC", y = "SSC", selected_rows = "V1")
#' stopifnot(inherits(p, "ggplot"))
#' print(p)
#' }
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
                             theme_family = "Helvetica") {

  # Create the base plot
  p <- dat %>%
    ggplot(aes_string(x = x, y = y)) +
    geom_hex(bins = hex_bins) +
    scale_fill_viridis(discrete = FALSE, trans = 'log') +
    scale_x_log10(breaks = c(0, 10, 100, 1000, 10000),
                  labels = trans_format("log10", scales::math_format(10^.x)),
                  limits = x_limits) +
    scale_y_log10(breaks = c(1, 10, 100, 1000, 10000),
                  labels = trans_format("log10", scales::math_format(10^.x)),
                  limits = y_limits) +
    theme_bw() +
    theme(text = element_text(family = theme_family))

  # Add highlighted points if selected_rows is provided
  if (!is.null(selected_rows) && "classes" %in% colnames(dat)) {
    p <- p +
      geom_point(data = dat[paste0('V', dat[,"classes"]) %in% selected_rows, ],
                 aes_string(x = x, y = y),
                 color = point_color,
                 alpha = point_alpha,
                 size = point_size)
  }

  # Add bin labels if bins is provided
  if (!is.null(bins) && !is.null(selected_rows)) {
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
#' @param common_legend Logical, whether to use a common legend for all plots. Defaults to `TRUE`.
#' @param ... Additional arguments forwarded to every MBR_flow_plot call.
#'
#' @details
#' plot_params is a list of lists, each containing x and y channel names. Other entries in a pair are ignored; additional plot options must be passed through ... and apply to ALL panels. selected_rows, bins and dat follow the same uppercase cluster-label and row-order rules as MBR_flow_plot. ncol and nrow control the arrangement and common_legend requests a shared legend. Shared legends do not ensure comparable density scales across channel pairs; each panel computes its own hexagon counts. No file is saved.
#'
#' @return
#' The arranged ggpubr plot object.
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' dat <- data.frame(FSC = runif(100, 1, 1000), SSC = runif(100, 1, 1000),
#'                   DNA = runif(100, 1, 1000))
#' p <- MBR_plot(dat, plot_params = list(list(x = "FSC", y = "SSC"),
#'                                       list(x = "FSC", y = "DNA")), nrow = 1)
#' stopifnot(inherits(p, "ggplot"))
#' print(p)
#' }
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
#' @details
#' Selects events using uppercase V-prefixed cluster labels derived from dat$classes, then removes classes from exported channels. Exports the values in dat as supplied, which may already have been transformed by MBR_process; it does not restore raw instrument values. Channel names must match the original flowFrame for parameter metadata matching. Renamed channels can produce unmatched metadata or warnings, and inherited parameter ranges may be inconsistent with transformed values.
#'
#' Rebuilds parameter fields and FCS description keys, preserving non-parameter description entries and updating total event and channel counts. Missing required metadata columns are filled with NA with a warning. Biobase and flowCore accessors are qualified by namespace. No event order correspondence to raw data is checked. Test exports by rereading and checking event counts and channels before downstream use.
#'
#' Creates rawdata_path/FilteredFCS and writes <original_basename>_filtered.fcs. Existing output files can be overwritten. The public output-directory argument is rawdata_path, not output_dir.
#'
#' @return
#' The output filename, returned invisibly.
#'
#' @examples
#' \dontrun{
#' library(flowCore)
#' library(Biobase)
#' fcs <- MBR_read("path/to/mapped_fcs")
#' # Keep original names and values for a conservative export example.
#' dat <- MBR_process(fcs, transformation = identity)
#' selected <- c("V1", "V2")
#' out <- tempfile("filtered-")
#' filename <- MBR_save(fcs, dat, selected, rawdata_path = out)
#' check <- flowCore::read.FCS(filename, transformation = FALSE)
#' stopifnot(nrow(flowCore::exprs(check)) == sum(paste0("V", dat$classes) %in% selected))
#' }
#'
#' @export
MBR_save <- function(fcs_files, dat, selected_rows, file_index = 1, rawdata_path) {
  # 1. Get the FlowFrame properly from flowSet using [[ ]]
  target_fcs_filename <- sampleNames(fcs_files)[file_index]
  original_flowframe <- fcs_files[[file_index]]

  # 2. Filter selected events
  selected_events_data <- dat[paste0('V', dat[,"classes"]) %in% selected_rows, ]
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
  original_desc <- Biobase::description(original_flowframe)
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
