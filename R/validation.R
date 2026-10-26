# Internal validation shared by sample-level analyses.
.mbr_inputs <- function(data, meta_data, group = NULL) {
  if (!(is.data.frame(data) || is.matrix(data)) || !nrow(data) || !ncol(data)) {
    stop("data must be a nonempty numeric matrix or data frame.")
  }
  ids <- rownames(data)
  if (!is.null(ids) && (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))) {
    stop("Sample IDs must be nonmissing, nonempty and unique.")
  }
  data <- as.data.frame(data)
  if (!all(vapply(data, is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(data)))) {
    stop("All features must be numeric and finite; handle missing values before analysis.")
  }
  if (!is.data.frame(meta_data) || nrow(meta_data) != nrow(data)) {
    stop("meta_data must be a data frame with one row per sample.")
  }
  if (is.null(ids) || identical(ids, as.character(seq_len(nrow(data)))) ||
      identical(rownames(meta_data), as.character(seq_len(nrow(meta_data))))) {
    warning("Sample IDs are unavailable: using positional metadata alignment. Set row names in both tables to verify alignment.", call. = FALSE)
  } else {
    meta_ids <- rownames(meta_data)
    if (anyDuplicated(ids) || anyDuplicated(meta_ids) || !setequal(ids, meta_ids)) {
      stop("Feature and metadata sample IDs must be unique and match exactly.")
    }
    meta_data <- meta_data[match(ids, meta_ids), , drop = FALSE]
  }
  if (!is.null(group)) {
    if (length(group) != 1L || !group %in% names(meta_data) ||
        anyNA(meta_data[[group]])) {
      stop("Supply a valid grouping column without missing values.")
    }
  }
  list(data = data, meta_data = meta_data)
}

# Average only held-out probabilities, giving each sample one diagnostic row.
.mbr_predictions <- function(model, positive) {
  predictions <- model$pred
  if (is.null(predictions) || !positive %in% names(predictions) ||
      any(!is.finite(predictions[[positive]]))) {
    stop("The model did not produce finite held-out class probabilities.")
  }
  probability <- stats::aggregate(predictions[[positive]],
                                  list(rowIndex = predictions$rowIndex), mean)
  names(probability)[2] <- "probability"
  rows <- match(probability$rowIndex, predictions$rowIndex)
  truth <- predictions$obs[rows]
  negative <- setdiff(levels(truth), positive)
  predicted <- factor(ifelse(probability$probability >= 0.5, positive, negative),
                      levels = levels(truth))
  data.frame(rowIndex = probability$rowIndex, probability = probability$probability,
             truth = truth, predicted = predicted)
}

.mbr_evaluate <- function(model, positive, train, test_data, test_meta_data, group) {
  if (is.null(test_data) && is.null(test_meta_data)) {
    return(.mbr_predictions(model, positive))
  }
  if (is.null(test_data) || is.null(test_meta_data)) {
    stop("Supply both test_data and test_meta_data.")
  }
  inputs <- .mbr_inputs(test_data, test_meta_data, group)
  features <- setdiff(names(train), "group")
  if (!setequal(names(inputs$data), features)) {
    stop("Test and training feature names must match exactly.")
  }
  truth <- factor(inputs$meta_data[[group]], levels = levels(train$group))
  if (anyNA(truth) || length(unique(truth)) != 2L) {
    stop("Test metadata must contain both training classes and no unknown classes.")
  }
  probability <- stats::predict(model, inputs$data[, features, drop = FALSE], type = "prob")[[positive]]
  negative <- setdiff(levels(truth), positive)
  data.frame(rowIndex = seq_along(truth), probability = probability, truth = truth,
             predicted = factor(ifelse(probability >= 0.5, positive, negative), levels = levels(truth)))
}

# Validate test inputs before fitting, including accidental sample reuse.
.mbr_test_inputs <- function(data, test_data, test_meta_data, group) {
  if (is.null(test_data) && is.null(test_meta_data)) return(invisible(NULL))
  if (is.null(test_data) || is.null(test_meta_data)) {
    stop("Supply both test_data and test_meta_data.")
  }
  inputs <- .mbr_inputs(test_data, test_meta_data, group)
  if (!setequal(colnames(data), names(inputs$data))) {
    stop("Test and training feature names must match exactly.")
  }
  train_ids <- rownames(data)
  test_ids <- rownames(test_data)
  usable <- function(ids, n) !is.null(ids) && !identical(ids, as.character(seq_len(n)))
  if (usable(train_ids, nrow(data)) && usable(test_ids, nrow(test_data)) &&
      length(intersect(train_ids, test_ids))) {
    stop("Independent test sample IDs must not overlap training sample IDs.")
  }
  invisible(NULL)
}

# Canonical cluster IDs used ONLY for matching, never to rename caller data.
# A cluster label is exactly v/V followed by digits: v12 and V12 become V12.
# Other feature names retain their spelling and remain case-sensitive. Leading
# zeros are preserved, so V01 and V1 remain distinct. Numeric SOM row names are
# not converted here; event classes are handled separately below.
# Set unique = FALSE for repeated event labels, where repeats are expected.
.mbr_cluster_ids <- function(labels, context = "Cluster labels", unique = TRUE) {
  labels <- as.character(labels)
  if (anyNA(labels) || any(!nzchar(labels))) {
    stop(context, " must be nonmissing and nonempty.", call. = FALSE)
  }
  ids <- sub("^[vV]([0-9]+)$", "V\\1", labels)
  if (unique && anyDuplicated(ids)) {
    stop(context, " contain duplicate or case-ambiguous IDs: ",
         paste(base::unique(ids[duplicated(ids)]), collapse = ", "), call. = FALSE)
  }
  ids
}

# Return indices in requested order after normalizing BOTH sides. This supports
# cross-table matches (e.g. V1 data and v1 p-value rows) without modifying either
# table. Ambiguous available/requested names fail rather than selecting the first
# match; missing IDs fail rather than yielding NA rows or an empty event export.
.mbr_cluster_match <- function(requested, available, context = "Available clusters") {
  wanted <- .mbr_cluster_ids(requested, "Requested clusters")
  supplied <- .mbr_cluster_ids(available, context)
  indices <- match(wanted, supplied)
  if (anyNA(indices)) {
    stop(context, " missing requested clusters: ",
         paste(requested[is.na(indices)], collapse = ", "), call. = FALSE)
  }
  indices
}

# SOM event classes may be numeric (1), numeric strings ("1"), or prefixed
# labels ("v1"/"V1"). Add a prefix only to wholly numeric strings. Validate
# distinct raw labels first so v1/V1 collisions are caught, while repeated
# occurrences of the same class remain valid and preserve event order.
.mbr_event_ids <- function(classes) {
  labels <- as.character(classes)
  distinct <- base::unique(labels)
  distinct <- ifelse(grepl("^[0-9]+$", distinct), paste0("V", distinct), distinct)
  .mbr_cluster_ids(distinct, "Event classes")
  labels <- ifelse(grepl("^[0-9]+$", labels), paste0("V", labels), labels)
  .mbr_cluster_ids(labels, "Event classes", unique = FALSE)
}
