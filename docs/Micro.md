# R/Micro.R walkthrough

This file defines the public analysis and FCS utilities. Sourcing it defines functions; their analyses run only when called. Package namespace import comments are used during package generation, not library calls when sourcing. Prefer loading an installed package to sourcing the script alone.

The main workflow is gated events → shared SOM → per-sample counts → row-normalized proportions → group tests, descriptive plots and optional predictive exploration. The FCS read/process/prepare/plot/save functions form a separate inspection and export workflow.

## `MBR_som` — Train a shared SOM and export sample cluster abundances

Pools a balanced random subset of gated events, trains a hexagonal self-organizing map, maps retained samples to its nodes, and writes counts and mapped FCS files.

Each named sample must contain name and a numeric data matrix with identical channels in identical order. The pooled training subset gives each sample an equal target contribution, capped by its available events. A SOM learns prototype vectors (codebook rows) representing similar multichannel events; mapping assigns each event to a best matching prototype. No explicit channel standardization is applied here. Channel units and transformations therefore affect the fit.

The grid side is round(sqrt(m)); the actual node count is its square. Despite the n_hclust argument, this implementation does not perform hierarchical clustering. Node labels and count levels are fixed at 2025. The default m = 2000 produces a 45 by 45 grid, matching that mapping; other grid sizes can omit nodes or create unused levels. Samples with fewer than 200000 events are removed before sampling and training. If none remain, the function stops.

Rare clusters are counted by feature column: a cluster is rare when its abundance is below 0.01% in every retained sample.

Writes cluster_information.csv, SOM.csv, count_tables/count_table_SOM.save, and one FCS file per sample in MappedFCS. Assigns cohonen_information, raw_count_table, and count_table in the global environment. The codebook spelling is intentionally retained. Seeds are not set internally; call set.seed before use. Existing output files may be replaced.

**Result:** No analytical object is returned; the last console call returns invisible NULL. Retrieve tables from the documented files or global variables.

```r
# Large example: use real gated samples with at least 200000 events each.
# gated <- list(sample1 = list(name = "sample1.fcs", data = events1),
#               sample2 = list(name = "sample2.fcs", data = events2))
set.seed(42)
out <- tempfile("som-")
MBR_som(gated_fcs = gated, out_path = out)
stopifnot(file.exists(file.path(out, "SOM.csv")),
          all(abs(rowSums(count_table) - 1) < 1e-8))
```

## `MBR_stat` — Test each abundance feature between sample groups

Computes one group-comparison p-value per feature, optionally adjusts for multiple testing, and exports significant feature columns.

Numeric matrices and data frames are accepted. Sample row names are matched to metadata row names and reordered when needed; mismatched IDs stop the analysis. Without usable IDs, positional matching emits a warning.

wilcox uses the unpaired, two-sided Wilcoxon rank-sum test for two groups. It compares rank distributions; interpreting it purely as a median comparison requires similarly shaped distributions. kruskal uses the Kruskal-Wallis rank test for two or more independent groups; a significant omnibus test does not identify the differing pairs. anova fits a one-way ANOVA testing equality of means, with independent errors and approximately normal residuals of similar variance. t.test (also ttest) uses the unpaired two-sided Welch test, allowing unequal variances, for two groups. None of these branches supports paired samples, covariate adjustment, or repeated measurements.

BH and fdr both request Benjamini-Hochberg adjustment, controlling false discovery rate under its dependence assumptions. bonferroni controls family-wise error through a more conservative adjustment. none, the default, leaves p-values unadjusted. Selection uses p.adj < cutoff, strictly excluding equality. Feature-wise tests on proportions remain affected by compositional dependence.

Wilcoxon tests use a normal approximation (exact = FALSE), including ties. Failed or non-finite feature tests emit a feature-specific warning and retain NA p-values; missing results are excluded from significance selection. Invalid test names stop before testing.

**Result:** The significant feature data frame, returned invisibly by the final assign call; pvalue_data is also assigned globally.

```r
set.seed(42)
data <- as.data.frame(matrix(runif(120, 0.01, 1), nrow = 20))
names(data) <- paste0("v", seq_len(ncol(data)))
data <- data / rowSums(data)
meta <- data.frame(Group = rep(c("A", "B"), each = 10))
out <- tempfile("MicroBiotR-")
dir.create(out)
MBR_stat(data, meta, "Group", correction = "BH", out_path = out)
stopifnot(nrow(pvalue_data) == ncol(data),
          all(pvalue_data$p.adj >= 0 & pvalue_data$p.adj <= 1),
          file.exists(file.path(out, "pvalue.csv")))
```

## `MBR_circle` — Draw group means in a circular heatmap

Displays feature means by group and a separate track for the difference between the first two group means.

Use a numeric data frame and row-aligned metadata. tapply orders the grouping levels; mean_table columns, rather than order of appearance in metadata, determine the subtraction first group minus second group. At least two groups are needed and the difference track is intended for two groups. Means are calculated without na.rm, so missing values propagate. The color midpoint is the midpoint of the overall minimum and maximum, not the overall mean. Constant means or differences can make color breaks or track limits degenerate.

This is a descriptive plot, with no hypothesis test or uncertainty interval. Group means can hide sample heterogeneity. The function draws on the active graphics device, resets circlize state, then redraws into circle.pdf. out_path must already exist.

**Result:** Invisible NULL from the final console call; the plot is drawn and saved.

```r
set.seed(42)
data <- as.data.frame(matrix(runif(120, 0.01, 1), nrow = 20))
names(data) <- paste0("v", seq_len(ncol(data)))
data <- data / rowSums(data)
meta <- data.frame(Group = rep(c("A", "B"), each = 10))
out <- tempfile("MicroBiotR-")
dir.create(out)
MBR_circle(data, meta, "Group", out_path = out)
stopifnot(file.exists(file.path(out, "circle.pdf")))
```

## `MBR_violin` — Plot one cluster abundance and an existing p-value

Draws a violin and box plot for one cluster and annotates a supplied raw or adjusted p-value.

Numeric matrices and data frames are accepted and converted to a data frame for named column extraction. cluster = 1 selects the lowercase column v1 and the row v1 of pvalue_data. Metadata rows must match data order. p specifies an existing p-value column, normally p.value or p.adj. No test is performed by this function; annotation validity depends on the provenance and correction of the supplied table.

A violin is a smoothed density estimate, and the overlaid box plot summarizes the median and interquartile range. The abundance label assumes data already contains proportions. Supply enough colors for all groups. The annotation is positioned at x = 1.5, appropriate for a two-group comparison. Writes violin_cluster_<cluster>_<p>.pdf to an existing directory and prints the plot.

**Result:** Invisible NULL; the plot is printed and saved, rather than returned.

```r
set.seed(42)
data <- as.data.frame(matrix(runif(120, 0.01, 1), nrow = 20))
names(data) <- paste0("v", seq_len(ncol(data)))
data <- data / rowSums(data)
meta <- data.frame(Group = rep(c("A", "B"), each = 10))
out <- tempfile("MicroBiotR-")
dir.create(out)
MBR_stat(data, meta, "Group", correction = "BH", out_path = out)
MBR_violin(data, meta, "Group", pvalue_data, p = "p.adj",
           cluster = 1, out_path = out)
stopifnot(file.exists(file.path(out, "violin_cluster_1_p.adj.pdf")))
```

## `MBR_beta` — Ordinate Bray-Curtis dissimilarity and test group differences

Calculates Bray-Curtis dissimilarities, a principal coordinates display, and a PERMANOVA group test with coordinate-wise comparison plots.

Rows are samples and columns are nonnegative features. Metadata must be row-aligned. For samples i and j, Bray-Curtis is sum(abs(x_i - x_j)) / sum(x_i + x_j); it emphasizes abundance differences and ignores joint absences. Empty samples and missing values need to be handled before calling. Raw counts retain library-size effects; this function does not normalize input.

Classical multidimensional scaling via cmdscale requests three axes and displays the first two. Use sufficiently many distinct samples to obtain three positive axes (at least four samples). Bray-Curtis need not be Euclidean, so negative eigenvalues can occur. The displayed percentages divide each eigenvalue by the sum of ALL eigenvalues; they are not necessarily conventional fractions of positive inertia. No correction for negative eigenvalues is requested.

adonis2 tests group-associated variation in the full distance matrix using dependency-default unrestricted permutations. Its R-squared summarizes the fraction of sums of squares attributed to group. Group dispersion differences can also affect the result; examine dispersion separately. Exchangeability is required, and this interface has no strata or covariate arguments. Set a seed for reproducible permutations.

The test argument changes only comparisons of plotted coordinate scores, not PERMANOVA. wilcox and ttest compare group pairs; kruskal is called through pairwise plot comparisons; anova gives an omnibus axis comparison. Axis tests do not replace a multivariate test. The wrapper specifies no multiple-testing correction for these annotations. Ellipses are visualization summaries, not a PERMANOVA confidence region. Writes .group.txt, adonis.txt, and pcoa_<test>.pdf to an existing out_path.

**Result:** The combined patchwork plot, returned visibly.

```r
set.seed(42)
data <- as.data.frame(matrix(runif(120, 0.01, 1), nrow = 20))
names(data) <- paste0("v", seq_len(ncol(data)))
data <- data / rowSums(data)
meta <- data.frame(Group = rep(c("A", "B"), each = 10))
out <- tempfile("MicroBiotR-")
dir.create(out)
p <- MBR_beta(data, out_path = out, meta_data = meta, group_name = "Group")
stopifnot(inherits(p, "patchwork"), file.exists(file.path(out, "adonis.txt")))
```

## `MBR_fs` — Select features with random forest recursive elimination

Uses cross-validated recursive feature elimination and produces feature importance and transformed abundance effect-size plots.

The wrapper sets the random seed to 2025. Row-aligned numeric features and a metadata grouping column are required. caret::rfe with rfFuncs evaluates requested subset sizes 1:rfe_size and, through caret defaults, the full predictor set using nfolds_cv-fold cross-validation; use rfe_size no greater than the feature count and enough observations per class for all folds. The selected subset is results$optVariables. top_n_features controls how many selected variables are displayed, not the size searched.

Random forests combine trees fitted to bootstrap samples and random predictor subsets. RFE ranks predictors and compares candidate subsets using resampling. Importance can be unstable when predictors are correlated. Effect-size plots use log10(value + 1), despite an axis label saying log10(abundance). Cohen's d is a standardized difference in group means; its sign depends on group ordering and ref_group. The displayed 95 percent intervals and effects are calculated after selecting features on these same data, so they are exploratory rather than selection-adjusted inference. Use two groups for this workflow and keep feature selection inside an outer resampling loop when evaluating predictive performance.

Writes .group.txt, figure1.txt (selected features plus group), and feature_exploration.pdf. Assigns MBR_selected_features using superassignment, normally into the global environment. Requires an existing output directory. Missing importance rows are removed with tidyr::drop_na. Three console plots are shown with readline prompts, so this workflow is interactive. The fitted RFE object is not returned.

**Result:** Invisible NULL; selected feature data are available as MBR_selected_features and files.

```r
set.seed(42)
data <- as.data.frame(matrix(runif(120, 0.01, 1), nrow = 20))
names(data) <- paste0("v", seq_len(ncol(data)))
data <- data / rowSums(data)
meta <- data.frame(Group = rep(c("A", "B"), each = 10))
out <- tempfile("MicroBiotR-")
dir.create(out)
MBR_fs(data, out_path = out, meta_data = meta, group_name = "Group",
       nfolds_cv = 2, rfe_size = 3, top_n_features = 3, ref_group = "A")
stopifnot(nrow(MBR_selected_features) == nrow(data),
          all(names(MBR_selected_features) %in% names(data)))
```

## `MBR_heatmap` — Show channel values for selected SOM prototypes

Selects codebook rows using feature column names and renders a configurable pheatmap.

Only the column names of data select clusters; its abundance values are not plotted. Both v1 and V1 are converted to V1 for lookup, so cohonen_information must have matching UPPERCASE V-prefixed row names. MBR_som's codebook normally has numeric row names; adapt a copy before passing it here. Missing matches yield missing rows.

The plotted matrix contains prototype channel values. Row scaling standardizes channels within each selected node; column scaling standardizes nodes within each channel; none preserves input units. These different choices change the interpretation of color. Constant rows or columns may fail when scaled. Optional clustering uses pheatmap defaults; it does not alter the package's original SOM assignments. A numeric conversion df_mat is computed but the original subset df is actually plotted. The argument boarder_color is intentionally spelled as in the implementation. Writes heatmap.pdf and draws the plot; out_path must exist.

**Result:** Invisible NULL; the pheatmap is drawn and saved.

```r
set.seed(42)
data <- as.data.frame(matrix(runif(120, 0.01, 1), nrow = 20))
names(data) <- paste0("v", seq_len(ncol(data)))
data <- data / rowSums(data)
meta <- data.frame(Group = rep(c("A", "B"), each = 10))
out <- tempfile("MicroBiotR-")
dir.create(out)
codebook <- as.data.frame(matrix(seq_len(18), nrow = 6))
rownames(codebook) <- paste0("V", 1:6)
MBR_heatmap(data, codebook, out_path = out, scale = "none")
stopifnot(file.exists(file.path(out, "heatmap.pdf")))
```

## `MBR_mantel` — Visualize Mantel associations and feature correlations

Compares distance patterns of metadata blocks with feature data using linkET and overlays associations on a feature correlation heatmap.

Samples must occupy matching rows of numeric data and metadata. clinical_cols and demographic_cols select metadata variable blocks; spec_select_names$A and $B supply their display names. NULL labels omit a block. Only appropriate numeric variables should enter distance calculations; arbitrary numeric encoding of categories gives arbitrary distances.

A Mantel statistic correlates entries of two sample-distance matrices and uses permutations for significance. Distances and permutation settings are delegated to linkET::mantel_test defaults, not specified by this wrapper; consult the installed linkET help to identify them. Unequal units across metadata variables can dominate distance, so consider scaling a copy of metadata before use. Permutations assume exchangeable samples, and no blocked design or confounder adjustment is exposed. This is an association analysis, not evidence of causality.

correlate(data) supplies the heatmap correlations (Pearson by the dependency default). Links bin Mantel r at 0.2 and 0.4 and p at 0.01 and 0.05; cut uses right-closed intervals, so boundary values enter the lower interval even though labels use less-than signs. Negative r values share the first category. The code specifies no multiple-testing adjustment. Writes mantel.pdf to an existing directory, prints the plot, and does not return the Mantel table.

**Result:** Invisible NULL; the association figure is saved and printed.

```r
set.seed(42)
data <- as.data.frame(matrix(runif(120, 0.01, 1), nrow = 20))
names(data) <- paste0("v", seq_len(ncol(data)))
data <- data / rowSums(data)
meta <- data.frame(Group = rep(c("A", "B"), each = 10))
out <- tempfile("MicroBiotR-")
dir.create(out)
# Include a strongly associated feature pair so the significance layer is exercised.
data$v2 <- data$v1 + seq_len(nrow(data)) * 1e-6
data <- data / rowSums(data)
clinical <- data.frame(age = seq_len(20), marker = runif(20), BMI = rnorm(20, 25, 2))
MBR_mantel(data, clinical, clinical_cols = c("age", "marker"),
           demographic_cols = "BMI", out_path = out)
stopifnot(file.exists(file.path(out, "mantel.pdf")))
```

## `MBR_ml` — Evaluate random forest classification with ROC

Fits a caret random forest classifier and saves a diagnostic plot based on held-out resampling predictions, or an optional independent test set.

Numeric predictors and row-aligned metadata with exactly two classes are required. Class labels must be valid R names for caret probability columns, and reference_level must identify a present class. Resampling settings are passed to caret::trainControl; caret tunes a random forest using ROC as its selection metric. The wrapper requires the randomForest backend. number and repeats configure the selected method; repeats is relevant to repeatedcv.

By default, plots use held-out caret predictions for the selected tuning parameters, averaged across repeats once per sample. These are cross-validation diagnostics, not independent validation after tuning or upstream feature selection. For independent evaluation, supply test_data and test_meta_data from samples unused in feature selection, preprocessing estimation and tuning. The fitted model, predictions and confusion metrics are returned invisibly.

Writes ml_group.txt and the PDF in out_path; creates the directory recursively if absent. Existing files can be overwritten.

reference_level is the positive class for both ROC probabilities and confusion metrics. ROC direction is fixed so larger probabilities indicate the positive class. No global warning options are modified. Resampling confidence intervals are descriptive and omit tuning and feature-selection uncertainty.

**Result:** An invisible list containing the fitted model, sample-level diagnostic predictions and confusion metrics; MBR_ml also returns the ROC object and plot.

```r
set.seed(42)
data <- as.data.frame(matrix(runif(120, 0.01, 1), nrow = 20))
names(data) <- paste0("v", seq_len(ncol(data)))
data <- data / rowSums(data)
meta <- data.frame(Group = rep(c("A", "B"), each = 10))
out <- tempfile("MicroBiotR-")
dir.create(out)
MBR_ml(data, meta, group_name = "Group", out_path = out,
       reference_level = "A", number = 2, repeats = 1)
stopifnot(file.exists(file.path(out, "roc_confusion.pdf")))
```

## `MBR_conf` — Evaluate random forest classification with a confusion heatmap

Fits a caret random forest classifier and saves a diagnostic plot based on held-out resampling predictions, or an optional independent test set.

Numeric predictors and row-aligned metadata with exactly two classes are required. Class labels must be valid R names for caret probability columns, and reference_level must identify a present class. Resampling settings are passed to caret::trainControl; caret tunes a random forest using ROC as its selection metric. The wrapper requires the randomForest backend. number and repeats configure the selected method; repeats is relevant to repeatedcv.

By default, plots use held-out caret predictions for the selected tuning parameters, averaged across repeats once per sample. These are cross-validation diagnostics, not independent validation after tuning or upstream feature selection. For independent evaluation, supply test_data and test_meta_data from samples unused in feature selection, preprocessing estimation and tuning. The fitted model, predictions and confusion metrics are returned invisibly.

Writes ml_group.txt and the PDF in out_path; creates the directory recursively if absent. Existing files can be overwritten.

Unlike MBR_ml, this function does not set a random seed; call set.seed before use. The heatmap labels are event counts (samples here), while fill is natural log(Freq + 1), making small cells visible alongside large ones. It shows no uncertainty intervals or ROC curve.

**Result:** An invisible list containing the fitted model, sample-level diagnostic predictions and confusion metrics; MBR_ml also returns the ROC object and plot.

```r
set.seed(42)
data <- as.data.frame(matrix(runif(120, 0.01, 1), nrow = 20))
names(data) <- paste0("v", seq_len(ncol(data)))
data <- data / rowSums(data)
meta <- data.frame(Group = rep(c("A", "B"), each = 10))
out <- tempfile("MicroBiotR-")
dir.create(out)
MBR_conf(data, meta, group_name = "Group", out_path = out,
       reference_level = "A", number = 2, repeats = 1)
stopifnot(file.exists(file.path(out, "confusion_matrix.pdf")))
```

## `MBR_reclustering` — Recluster numeric features using Ward linkage

Standardizes numeric columns, performs Euclidean hierarchical clustering, and appends cluster membership to a copy of the input.

Non-numeric columns are preserved but excluded from distance calculations. All numeric columns, including numeric identifiers or previous labels, enter clustering; remove unwanted numeric fields from a copy beforehand. scale centers each feature and divides by its sample standard deviation, giving differently measured channels comparable weight. Constant columns, missing or infinite values cause invalid distances.

Ward.D2 linkage merges groups using a criterion related to the increase in within-cluster sum of squares; with Euclidean distances it favors compact clusters. cutree chooses the requested number of groups, which must be between 1 and the number of rows. Labels are arbitrary identifiers, not ordered biological states. This descriptive partition has no p-value or estimate of the optimal number of clusters. A pre-existing Cluster column is overwritten in the returned copy. Assigns reclustered_information globally.

**Result:** The input data frame with an added or replaced Cluster column, returned invisibly and assigned to reclustered_information.

```r
x <- data.frame(channel1 = c(1, 2, 8, 9), channel2 = c(2, 1, 9, 8),
                label = letters[1:4])
y <- MBR_reclustering(x, num_clusters = 2)
stopifnot(nrow(y) == nrow(x), length(unique(y$Cluster)) == 2,
          identical(y$label, x$label))
```

## `MBR_read` — Read FCS files from a directory

Reads matching files into a flowCore flowSet for event processing.

Lists non-recursive filenames matching fcs$ (case-sensitive, without requiring a dot before fcs). Uppercase .FCS files are not selected. Files are read with alter.names = FALSE and channels matching an asterisk, Bits, or Drop excluded. The wrapper leaves transformation and other reader settings at flowCore defaults; it does not reproduce Gating's explicit transformation = FALSE. It stops for a missing directory or no matching files. Reading is not gating or compensation; verify instrument channels and reader settings before analysis.

**Result:** A flowCore flowSet with one flowFrame per selected file.

```r
fcs <- MBR_read("path/to/fcs_directory")
stopifnot(inherits(fcs, "flowSet"), length(fcs) > 0)
```

## `MBR_process` — Transform and rename events from one FCS sample

Extracts one flowSet sample, transforms its channel values, and optionally renames columns.

The default transformation is 10^((4*x)/65000), an exponential rescaling tied to the original instrument range. It is not a universal cytometry compensation, logarithm, or arcsinh transformation. Choose a transformation suitable for your acquisition scale. Every column except a column named exactly classes is transformed; cluster labels are retained as labels. The selected sample is found through sampleNames and frames.

column_mapping is a named character vector whose names are old channel names and values are replacement names. Unmatched names are ignored. Renaming happens after transformation. The function does not gate events or normalize sample abundance. Preserve original channel names when planning MBR_save, whose parameter matching uses those names.

**Result:** A data frame of transformed events; rows retain event order and classes is preserved when present.

```r
ff <- flowCore::flowFrame(matrix(c(10, 20, 30, 40, 1, 2), nrow = 2,
                                 dimnames = list(NULL, c("FSC", "SSC", "classes"))))
fs <- flowCore::flowSet(list(sample1 = ff))
x <- MBR_process(fs, transformation = identity, column_mapping = c(FSC = "scatter"))
stopifnot("scatter" %in% names(x), identical(as.numeric(x$classes), c(1, 2)))
```

## `MBR_prepare` — Select and transform prototype rows for plotting

Filters a prototype table by row name, optionally renames channels, and transforms all remaining columns.

selected_rows is matched against existing row names using membership, retaining the INPUT order rather than the order of selected_rows. Missing requested rows are silently omitted. Unlike MBR_process, every column is transformed, including any label column supplied accidentally. Keep only appropriate numeric channels. The default exponential transformation assumes the original 65000 instrument scale and should match the transformation used for plotted events. Renaming occurs before transformation; a named character vector maps old names to new names. For plot labels, explicitly align selected_rows to the resulting row order.

**Result:** A data frame containing the selected, renamed and transformed rows.

```r
bins <- data.frame(FSC = c(10, 20, 30), SSC = c(40, 50, 60),
                   row.names = c("V1", "V2", "V3"))
x <- MBR_prepare(bins, c("V3", "V1"), transformation = identity,
                 column_mapping = c(FSC = "scatter"))
stopifnot(identical(rownames(x), c("V1", "V3")), "scatter" %in% names(x))
```

## `MBR_flow_plot` — Plot event density with optional cluster highlighting

Builds a two-channel hexagonal density plot with logarithmic axes, highlighted events, and optional prototype labels.

dat must contain positive numeric x and y channels for log axes. geom_hex counts events per hexagon and the fill scale is logarithmic. hex_bins sets the spatial resolution, so colors summarize counts per bin rather than normalized probability. The hexbin backend must be installed. Events outside scale limits are removed for plot calculations.

Highlight matching uses paste0("V", dat$classes), requiring uppercase labels such as V1 in selected_rows. Other functions may use lowercase v abundance columns; these conventions are not interchangeable without conversion. bins must have x and y columns and exactly one correctly ordered row per label in selected_rows. No statistical test is performed. This returns a plot without saving or printing it.

**Result:** A ggplot object that can be printed or saved with ggplot2::ggsave.

```r
set.seed(42)
dat <- data.frame(FSC = runif(100, 1, 1000), SSC = runif(100, 1, 1000),
                  classes = rep(1:2, 50))
p <- MBR_flow_plot(dat, x = "FSC", y = "SSC", selected_rows = "V1")
stopifnot(inherits(p, "ggplot"))
print(p)
```

## `MBR_plot` — Arrange multiple flow density plots

Calls MBR_flow_plot for each channel pair and arranges the resulting plots with ggpubr.

plot_params is a list of lists, each containing x and y channel names. Other entries in a pair are ignored; additional plot options must be passed through ... and apply to ALL panels. selected_rows, bins and dat follow the same uppercase cluster-label and row-order rules as MBR_flow_plot. ncol and nrow control the arrangement and common_legend requests a shared legend. Shared legends do not ensure comparable density scales across channel pairs; each panel computes its own hexagon counts. No file is saved.

**Result:** The arranged ggpubr plot object.

```r
set.seed(42)
dat <- data.frame(FSC = runif(100, 1, 1000), SSC = runif(100, 1, 1000),
                  DNA = runif(100, 1, 1000))
p <- MBR_plot(dat, plot_params = list(list(x = "FSC", y = "SSC"),
                                      list(x = "FSC", y = "DNA")), nrow = 1)
stopifnot(inherits(p, "ggplot"))
print(p)
```

## `MBR_save` — Export selected cluster events to a new FCS file

Filters processed events by cluster and rebuilds a flowFrame using original FCS parameter metadata.

Selects events using uppercase V-prefixed cluster labels derived from dat$classes, then removes classes from exported channels. Exports the values in dat as supplied, which may already have been transformed by MBR_process; it does not restore raw instrument values. Channel names must match the original flowFrame for parameter metadata matching. Renamed channels can produce unmatched metadata or warnings, and inherited parameter ranges may be inconsistent with transformed values.

Rebuilds parameter fields and FCS description keys, preserving non-parameter description entries and updating total event and channel counts. Missing required metadata columns are filled with NA with a warning. Biobase and flowCore accessors are qualified by namespace. No event order correspondence to raw data is checked. Test exports by rereading and checking event counts and channels before downstream use.

Creates rawdata_path/FilteredFCS and writes <original_basename>_filtered.fcs. Existing output files can be overwritten. The public output-directory argument is rawdata_path, not output_dir.

**Result:** The output filename, returned invisibly.

```r
library(flowCore)
library(Biobase)
fcs <- MBR_read("path/to/mapped_fcs")
# Keep original names and values for a conservative export example.
dat <- MBR_process(fcs, transformation = identity)
selected <- c("V1", "V2")
out <- tempfile("filtered-")
filename <- MBR_save(fcs, dat, selected, rawdata_path = out)
check <- flowCore::read.FCS(filename, transformation = FALSE)
stopifnot(nrow(flowCore::exprs(check)) == sum(paste0("V", dat$classes) %in% selected))
```

## Local functions

`stat_func(formula, data)` exists only in the ANOVA branch of `MBR_beta`. It fits `aov`, extracts the first F-test p-value, and is called for the first and second PCoA axes; the resulting `pval_x` and `pval_y` variables are not used in the final plot, which adds its own `stat_compare_means` layer.

`show_plots(plots)` exists only inside `MBR_fs`. It prints each plot and calls `readline` after every one, including the last. It returns the loop result invisibly and does not save plots; the enclosing function writes the PDF first. These local functions have no standalone `?` API. Their documentation belongs to the enclosing functions.
