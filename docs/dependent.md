# R/dependent.R walkthrough

This file defines support functions called by the SOM pipeline. Sourcing it defines functions and does not run training. The sequence is downsample → compute_som → map_som → assign_clusters → count_observations → get_counts. Sample-level counting remains separate from pooled model training so the final experimental unit is a sample rather than an individual event.

## `downsample` — Pool a balanced random subset of sample events

Randomly samples up to an equal quota from each selected sample and row-binds the results.

data is a named list of samples, each containing a data matrix with identical channels in the same order. n is a target for the combined number of events, not a per-sample size. The quota is n divided by the number of selected samples; smaller samples contribute all available rows and unused quotas are not redistributed. Fractional quotas are passed to sample and can be truncated. Sampling is without replacement and the function does not set a seed.

The samples argument filters membership but retains input-list order. Select at least one existing nonempty sample. Supply matrices with at least two columns because row extraction does not use drop = FALSE. Pooling discards sample identity; use the original sample list for later sample-level counts. Balanced sampling reduces domination by deeply measured samples but does not remove biological or batch confounding.

**Result:** A pooled event matrix, normally with at most n rows; no labels or sample identifiers are added.

**Help:** `?downsample` after installing the package.

Example and explicit checks:

```r
set.seed(42)
a <- matrix(seq_len(40), ncol = 2)
x <- list(A = list(name = "A", data = a), B = list(name = "B", data = a + 1))
y <- downsample(x, n = 12)
stopifnot(is.matrix(y), nrow(y) == 12, ncol(y) == 2)
```

## `compute_som` — Fit a hexagonal self-organizing map

Builds a square hexagonal kohonen grid and trains prototypes from a sample data matrix.

data is a list with numeric data matrix and optional sample metadata. The grid side is round(sqrt(nrow(data$data)/n_cells)); n_cells is a grid-sizing target, not a guaranteed occupancy. Ensure it yields a positive grid side. Training calls kohonen::som with keep.data = TRUE and otherwise dependency-default training settings. No explicit scaling, missing-value treatment, or seed is supplied. Channels with larger scales can dominate prototype distances.

A SOM represents high-dimensional events with codebook vectors on a neighborhood-preserving grid. It is an unsupervised representation and does not imply biologically validated populations or significance tests. The input list's data field is removed from the returned copy, while training data remain stored inside som because keep.data is true.

**Result:** A copy of data with som containing the trained kohonen object and with the top-level data field removed.

**Help:** `?compute_som` after installing the package.

Example and explicit checks:

```r
set.seed(42)
x <- list(name = "training", data = matrix(runif(200), ncol = 2))
y <- compute_som(x, n_cells = 25)
stopifnot(inherits(y$som, "kohonen"), is.null(y$data), nrow(y$som$grid$pts) == 4)
```

## `assign_clusters` — Translate SOM node indices into cluster labels

Looks up each event node index in the first column of a node-to-cluster mapping.

Uses som$unit.classif when present, otherwise classes. A present SOM assignment takes precedence over existing classes. clusters is a matrix or data frame with one row per SOM node; row POSITION, not row name, controls lookup. The first column supplies cluster labels; further columns are ignored. Node indices are converted with as.numeric, so factors can be interpreted as level codes rather than displayed labels. Invalid lookup results may generate NA and a warning. No clustering algorithm is run here. This helper is internal: use MicroBiotR:::assign_clusters rather than a direct package export.

**Result:** A copy of data whose classes field contains the looked-up cluster labels.

**Help:** `?assign_clusters` after installing the package.

Example and explicit checks:

```r
x <- list(name = "sample", classes = c(1, 2, 1, 3))
y <- MicroBiotR:::assign_clusters(x, matrix(c(1, 1, 2), ncol = 1))
stopifnot(identical(as.numeric(y$classes), c(1, 1, 1, 2)))
```

## `map_som` — Assign sample events to a trained SOM

Optionally samples events and maps them to the best matching nodes of a trained kohonen model.

data is a sample list with a numeric data matrix and normally a name. trained is the kohonen object itself, such as compute_som(...)$som, not the entire sample list. Match its channel number, ordering, units and preprocessing. When n_subset is smaller than the sample size, sample.int selects rows without replacement; otherwise all rows are retained. Use a positive integer limit and call set.seed for reproducible subsampling. Row extraction does not set drop = FALSE, so one-row or one-channel subsets can lose matrix shape.

The trained codebook is fixed during mapping. classes contains SOM node IDs rather than biologically curated cluster labels; use assign_clusters for a separate metacluster mapping if appropriate.

**Result:** A copy of data with data containing the retained events and classes containing their SOM node indices.

**Help:** `?map_som` after installing the package.

Example and explicit checks:

```r
set.seed(42)
x <- list(name = "sample", data = matrix(runif(200), ncol = 2))
trained <- compute_som(x, n_cells = 25)$som
y <- map_som(x, trained, n_subset = 20)
stopifnot(nrow(y$data) == 20, length(y$classes) == 20,
          all(y$classes %in% seq_len(nrow(trained$grid$pts))))
```

## `count_observations` — Count sample events over a full cluster index range

Counts classes into bins 1 through the requested maximum and stores a one-column count matrix.

clusters[1] is converted to a positive number and used as BOTH the count-resolution name and the maximum cluster ID. Supply a single positive integer or numeric string, not a vector of cluster labels. Classes outside 1:max_cluster_id, including missing values, do not contribute to the table. Empty valid bins receive zero counts. No warning identifies discarded classes, so verify that the sum equals the intended number of events. data$name supplies the output column name. Existing counts are replaced, not extended. Counts describe mapped events and are not normalized abundances.

**Result:** A copy of data with counts holding one cluster-by-one-sample matrix, named by the requested resolution.

**Help:** `?count_observations` after installing the package.

Example and explicit checks:

```r
x <- list(name = "sample", classes = c(1, 1, 3))
y <- count_observations(x, clusters = "4")
stopifnot(identical(as.numeric(y$counts[["4"]]), c(2, 0, 1, 0)),
          sum(y$counts[["4"]]) == length(x$classes))
```

## `get_counts` — Combine per-sample cluster count matrices

Collects count resolutions across samples, aligns cluster rows, and combines count columns.

data is a named list of sample lists with counts produced by count_observations. For each resolution, samples lacking it are omitted; they do not get a zero column. Cluster IDs are the numeric row names of per-sample count matrices, sorted numerically. Missing cluster rows in a included sample are filled with zeros. Each input count matrix is treated as a single column. Result row names are lowercase v-prefixed IDs and column names come from outer sample-list names, not data$name. This establishes the abundance-column convention used by MBR_violin after transposing. No normalization or hypothesis test is performed.

**Result:** A named list of numeric cluster-by-sample matrices, one per available count resolution.

**Help:** `?get_counts` after installing the package.

Example and explicit checks:

```r
a <- count_observations(list(name = "A", classes = c(1, 1, 3)), "4")
b <- count_observations(list(name = "B", classes = c(2, 3)), "4")
y <- get_counts(list(A = a, B = b))
stopifnot(identical(dim(y[["4"]]), c(4L, 2L)),
          identical(rownames(y[["4"]]), paste0("v", 1:4)),
          identical(as.numeric(colSums(y[["4"]])), c(3, 2)))
```

