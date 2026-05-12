# ====================================================================
# Helper Function Definitions
# ====================================================================

#' Data Downsampling Function
#' @param data Input data list
#' @param n Total sampling size
#' @param samples Sample names to be sampled
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

#' SOM Computation Function
#' @param data Input data
#' @param n_cells Expected cells per node
#' @export
compute_som <- function(data, n_cells = 100) {
  require(kohonen)
  
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

#' Cluster Label Assignment Function
#' @param data Input data
#' @param clusters Cluster definition matrix
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

#' Map to SOM Function
#' @param data Input data
#' @param trained Trained SOM
#' @param n_subset Subset size
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

#' Count Observations Function
#' @param data Input data
#' @param clusters Cluster names
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

#' Get Count Tables Function
#' @param data Input data list
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