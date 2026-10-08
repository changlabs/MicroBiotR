# Script guide: see R/README.md and docs/ for workflow, statistical interpretation,
# examples, side effects, and known implementation limitations.
# Helpers support the analysis functions described in docs/.

# ====================================================================
# Helper Function Definitions
# ====================================================================

#' Pool a balanced random subset of sample events
#'
#' @description
#' Randomly samples up to an equal quota from each selected sample and row-binds the results.
#'
#' @param data Named list of sample lists, each containing a numeric data matrix with matching channels.
#' @param n Target total pooled event count, divided equally across selected samples; default 1000000.
#' @param samples Character vector of names in data to include; defaults to all sample names.
#'
#' @details
#' data is a named list of samples, each containing a data matrix with identical channels in the same order. n is a target for the combined number of events, not a per-sample size. The quota is n divided by the number of selected samples; smaller samples contribute all available rows and unused quotas are not redistributed. Fractional quotas are passed to sample and can be truncated. Sampling is without replacement and the function does not set a seed.
#'
#' The samples argument filters membership but retains input-list order. Select at least one existing nonempty sample. Supply matrices with at least two columns because row extraction does not use drop = FALSE. Pooling discards sample identity; use the original sample list for later sample-level counts. Balanced sampling reduces domination by deeply measured samples but does not remove biological or batch confounding.
#'
#' @return
#' A pooled event matrix, normally with at most n rows; no labels or sample identifiers are added.
#'
#' @examples
#' set.seed(42)
#' a <- matrix(seq_len(40), ncol = 2)
#' x <- list(A = list(name = "A", data = a), B = list(name = "B", data = a + 1))
#' y <- downsample(x, n = 12)
#' stopifnot(is.matrix(y), nrow(y) == 12, ncol(y) == 2)
#'
#' @export
downsample <- function(data, n = 1e6, samples = names(data)) {
  cat("Sampling from", length(samples), "samples...\n")
  
  # Filter specified samples
  tmp_data <- data[sapply(names(data), function(x) x %in% samples)]
  
  # Calculate sampling size per sample
  n_per_sample <- n / length(tmp_data)
  cat("Sampling size per sample:", round(n_per_sample), "\n")
  
  # Execute sampling
  tmp_sample <- lapply(tmp_data, function(x) {
    tmp_mat <- x[["data"]]
    actual_n <- min(nrow(tmp_mat), n_per_sample)
    rows <- sample(seq_len(nrow(tmp_mat)), size = actual_n)
    tmp_mat[rows, ]
  })
  
  # Combine all samples
  tmp_sample <- do.call(rbind, tmp_sample)
  return(tmp_sample)
}

#' Fit a hexagonal self-organizing map
#'
#' @description
#' Builds a square hexagonal kohonen grid and trains prototypes from a sample data matrix.
#'
#' @param data Sample list containing the fields described in Details.
#' @param n_cells Target average events per node used to size the square grid; default 100.
#'
#' @details
#' data is a list with numeric data matrix and optional sample metadata. The grid side is round(sqrt(nrow(data$data)/n_cells)); n_cells is a grid-sizing target, not a guaranteed occupancy. Ensure it yields a positive grid side. Training calls kohonen::som with keep.data = TRUE and otherwise dependency-default training settings. No explicit scaling, missing-value treatment, or seed is supplied. Channels with larger scales can dominate prototype distances.
#'
#' A SOM represents high-dimensional events with codebook vectors on a neighborhood-preserving grid. It is an unsupervised representation and does not imply biologically validated populations or significance tests. The input list's data field is removed from the returned copy, while training data remain stored inside som because keep.data is true.
#'
#' @return
#' A copy of data with som containing the trained kohonen object and with the top-level data field removed.
#'
#' @examples
#' set.seed(42)
#' x <- list(name = "training", data = matrix(runif(200), ncol = 2))
#' y <- compute_som(x, n_cells = 25)
#' stopifnot(inherits(y$som, "kohonen"), is.null(y$data), nrow(y$som$grid$pts) == 4)
#'
#' @export
compute_som <- function(data, n_cells = 100) {
  
  tmp_data <- data
  tmp_mat <- tmp_data$data
  
  # Calculate grid dimensions
  dimno <- round(sqrt(nrow(tmp_mat)) / sqrt(n_cells), 0)
  cat("SOM grid dimensions:", dimno, "x", dimno, "\n")
  
  # Create hexagonal grid
  grd <- somgrid(xdim = dimno, ydim = dimno, topo = "hexagonal")
  
  # Train SOM
  cat("Starting SOM training...\n")
  tmp_data[['som']] <- som(tmp_data$data, grid = grd, keep.data = TRUE)
  tmp_data$data <- NULL
  
  return(tmp_data)
}

#' Translate SOM node indices into cluster labels
#'
#' @description
#' Looks up each event node index in the first column of a node-to-cluster mapping.
#'
#' @param data Sample list containing the fields described in Details.
#' @param clusters Matrix or data frame whose first column maps node row positions to cluster labels.
#'
#' @details
#' Uses som$unit.classif when present, otherwise classes. A present SOM assignment takes precedence over existing classes. clusters is a matrix or data frame with one row per SOM node; row POSITION, not row name, controls lookup. The first column supplies cluster labels; further columns are ignored. Node indices are converted with as.numeric, so factors can be interpreted as level codes rather than displayed labels. Invalid lookup results may generate NA and a warning. No clustering algorithm is run here. This helper is internal: use MicroBiotR:::assign_clusters rather than a direct package export.
#'
#' @return
#' A copy of data whose classes field contains the looked-up cluster labels.
#'
#' @examples
#' x <- list(name = "sample", classes = c(1, 2, 1, 3))
#' y <- MicroBiotR:::assign_clusters(x, matrix(c(1, 1, 2), ncol = 1))
#' stopifnot(identical(as.numeric(y$classes), c(1, 1, 1, 2)))
#'
assign_clusters <- function(data, clusters) {
  tmp_data <- data
  
  # Get SOM node assignments
  if (exists("som", where = tmp_data) && !is.null(tmp_data$som$unit.classif)) {
    som_node_assignments <- tmp_data$som$unit.classif
  } else if (exists("classes", where = tmp_data) && !is.null(tmp_data$classes)) {
    som_node_assignments <- tmp_data$classes
  } else {
    stop("Data does not contain 'som$unit.classif' or 'classes' information for clustering.")
  }
  
  som_node_assignments_numeric <- as.numeric(som_node_assignments)
  
  # Validate cluster parameters
  if (!is.matrix(clusters) && !is.data.frame(clusters)) {
    stop("The 'clusters' parameter must be a matrix or data frame defining SOM node to cluster mapping.")
  }
  if (ncol(clusters) < 1) {
    stop("The 'clusters' matrix/data frame must have at least one column for cluster IDs.")
  }
  
  # Assign new cluster labels (with progress bar)
  new_cluster_assignments <- pbapply::pbsapply(som_node_assignments_numeric, function(x) {
    clusters[x, 1]
  })
  
  # Check for NA values
  if (any(is.na(new_cluster_assignments))) {
    warning("Some SOM node assignments resulted in NA (out of bounds for cluster definitions). Check SOM size vs. cluster mapping.")
  }
  
  tmp_data$classes <- new_cluster_assignments
  return(tmp_data)
}

#' Assign sample events to a trained SOM
#'
#' @description
#' Optionally samples events and maps them to the best matching nodes of a trained kohonen model.
#'
#' @param data Sample list containing the fields described in Details.
#' @param trained Trained kohonen object with compatible channel order and preprocessing.
#' @param n_subset Positive integer event limit; samples larger than this are randomly reduced.
#'
#' @details
#' data is a sample list with a numeric data matrix and normally a name. trained is the kohonen object itself, such as compute_som(...)$som, not the entire sample list. Match its channel number, ordering, units and preprocessing. When n_subset is smaller than the sample size, sample.int selects rows without replacement; otherwise all rows are retained. Use a positive integer limit and call set.seed for reproducible subsampling. Row extraction does not set drop = FALSE, so one-row or one-channel subsets can lose matrix shape.
#'
#' The trained codebook is fixed during mapping. classes contains SOM node IDs rather than biologically curated cluster labels; use assign_clusters for a separate metacluster mapping if appropriate.
#'
#' @return
#' A copy of data with data containing the retained events and classes containing their SOM node indices.
#'
#' @examples
#' set.seed(42)
#' x <- list(name = "sample", data = matrix(runif(200), ncol = 2))
#' trained <- compute_som(x, n_cells = 25)$som
#' y <- map_som(x, trained, n_subset = 20)
#' stopifnot(nrow(y$data) == 20, length(y$classes) == 20,
#'           all(y$classes %in% seq_len(nrow(trained$grid$pts))))
#'
#' @export
map_som <- function(data, trained, n_subset) {
  tmp_data <- data
  
  # Subsample (if needed)
  if (n_subset < nrow(tmp_data$data)) {
    tmp_mat <- tmp_data$data[sample.int(nrow(tmp_data$data), size = n_subset), ]
  } else {
    tmp_mat <- tmp_data$data
  }
  
  # Map to trained SOM
  tmp_som <- kohonen::map(trained, tmp_mat)
  tmp_data[["data"]] <- tmp_mat
  tmp_data[["classes"]] <- tmp_som$unit.classif
  
  return(tmp_data)
}

#' Count sample events over a full cluster index range
#'
#' @description
#' Counts classes into bins 1 through the requested maximum and stores a one-column count matrix.
#'
#' @param data Sample list containing the fields described in Details.
#' @param clusters A single positive integer or numeric string giving both the resolution name and maximum cluster ID.
#'
#' @details
#' clusters[1] is converted to a positive number and used as BOTH the count-resolution name and the maximum cluster ID. Supply a single positive integer or numeric string, not a vector of cluster labels. Classes outside 1:max_cluster_id, including missing values, do not contribute to the table. Empty valid bins receive zero counts. No warning identifies discarded classes, so verify that the sum equals the intended number of events. data$name supplies the output column name. Existing counts are replaced, not extended. Counts describe mapped events and are not normalized abundances.
#'
#' @return
#' A copy of data with counts holding one cluster-by-one-sample matrix, named by the requested resolution.
#'
#' @examples
#' x <- list(name = "sample", classes = c(1, 1, 3))
#' y <- count_observations(x, clusters = "4")
#' stopifnot(identical(as.numeric(y$counts[["4"]]), c(2, 0, 1, 0)),
#'           sum(y$counts[["4"]]) == length(x$classes))
#'
#' @export
count_observations <- function(data, clusters) {
  tmp_data <- data
  tmp_classes <- data$classes
  
  # Get cluster resolution name
  cluster_resolution_name <- clusters[1]
  max_cluster_id <- as.numeric(cluster_resolution_name)
  
  if (is.na(max_cluster_id) || max_cluster_id <= 0) {
    stop("The 'clusters' parameter must be a single string that can be converted to a positive number (e.g., '1024').")
  }
  
  # Ensure all possible cluster IDs are represented
  all_possible_cluster_levels <- as.character(1:max_cluster_id)
  
  # Count each cluster
  tmp_mat <- table(factor(tmp_classes, levels = all_possible_cluster_levels))
  
  # Convert to matrix
  tmp_mat <- as.matrix(tmp_mat)
  colnames(tmp_mat) <- data$name
  
  # Store count matrix
  tmp_counts <- list()
  tmp_counts[[cluster_resolution_name]] <- tmp_mat
  
  tmp_data$counts <- tmp_counts
  return(tmp_data)
}

#' Combine per-sample cluster count matrices
#'
#' @description
#' Collects count resolutions across samples, aligns cluster rows, and combines count columns.
#'
#' @param data Named list of sample lists, each containing a numeric data matrix with matching channels.
#'
#' @details
#' data is a named list of sample lists with counts produced by count_observations. For each resolution, samples lacking it are omitted; they do not get a zero column. Cluster IDs are the numeric row names of per-sample count matrices, sorted numerically. Missing cluster rows in a included sample are filled with zeros. Each input count matrix is treated as a single column. Result row names are lowercase v-prefixed IDs and column names come from outer sample-list names, not data$name. This establishes the abundance-column convention used by MBR_violin after transposing. No normalization or hypothesis test is performed.
#'
#' @return
#' A named list of numeric cluster-by-sample matrices, one per available count resolution.
#'
#' @examples
#' a <- count_observations(list(name = "A", classes = c(1, 1, 3)), "4")
#' b <- count_observations(list(name = "B", classes = c(2, 3)), "4")
#' y <- get_counts(list(A = a, B = b))
#' stopifnot(identical(dim(y[["4"]]), c(4L, 2L)),
#'           identical(rownames(y[["4"]]), paste0("v", 1:4)),
#'           identical(as.numeric(colSums(y[["4"]])), c(3, 2)))
#'
#' @export
get_counts <- function(data) {
  # Extract count list from each sample
  tmp_cl_counts <- lapply(data, "[[", "counts")
  
  # Get all unique cluster resolution names
  tmp_cl_names <- unique(unlist(lapply(tmp_cl_counts, names)))
  
  tmp_counts <- lapply(tmp_cl_names, function(n) {
    tmp_mat_list <- lapply(tmp_cl_counts, "[[", n)
    tmp_mat_list <- tmp_mat_list[!sapply(tmp_mat_list, is.null)]
    
    if (length(tmp_mat_list) == 0) {
      warning(paste0("No count data found for cluster resolution: ", n))
      return(NULL)
    }
    
    # Get all cluster IDs
    all_cluster_ids <- unique(unlist(lapply(tmp_mat_list, rownames)))
    all_cluster_ids <- sort(as.numeric(all_cluster_ids))
    
    # Standardize count lists
    standardized_counts_list <- lapply(tmp_mat_list, function(mat) {
      standard_vec <- rep(0, length(all_cluster_ids))
      names(standard_vec) <- all_cluster_ids
      if (length(rownames(mat)) > 0) {
        standard_vec[rownames(mat)] <- mat[, 1]
      }
      return(standard_vec)
    })
    
    # Combine matrices
    combined_mat <- do.call(cbind, standardized_counts_list)
    
    # Add "v" prefix to row names
    rownames(combined_mat) <- paste0("v", rownames(combined_mat))
    colnames(combined_mat) <- names(tmp_mat_list)
    
    return(combined_mat)
  })
  
  names(tmp_counts) <- tmp_cl_names
  return(tmp_counts)
}