# Statistical interpretation

## Experimental unit and proportions

Individual flow events train the SOM, but downstream group comparisons use one abundance row per biological sample. Events from the same sample are not independent biological replicates. Functions match and reorder metadata by sample row names. IDs must match exactly; tables without usable IDs trigger a positional-alignment warning.

Relative abundances sum to one, so features are constrained and correlated. A change in one population can change another population's proportion without changing its absolute cell count. These scripts do not implement a specialized compositional model, absolute-load correction, batch correction or covariate adjustment. Document gating, transformations, sample exclusions and normalization before interpreting results.

## SOM and hierarchical partitions

A SOM learns multichannel prototype vectors on a grid, then maps each event to a matching prototype. Neighborhood structure helps represent similar signal profiles, while occupancy gives a sample-level feature table. The code does not standardize channels before training, so channel scale affects distance. Prototype IDs are descriptive labels, not confirmed taxa. See the [kohonen package manual](https://cran.r-project.org/web/packages/kohonen/kohonen.pdf) and installed `?kohonen::som` for dependency defaults.

Reclustering is a separate operation: numeric columns are centered and divided by their standard deviations, Euclidean distances are computed, and Ward.D2 builds a dendrogram. A cut at a user-chosen number of groups gives a partition; the function does not test that this number is optimal. See [R's hierarchical clustering reference](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/hclust.html).

## Feature-wise tests in MBR_stat

* **Wilcoxon:** two independent groups, two-sided rank-sum test. It is not automatically a test of medians under differently shaped distributions. The wrapper uses the normal approximation (`exact = FALSE`) so tied observations retain a valid result. Failed or non-finite tests are reported as NA with a feature-specific warning. See [R's Wilcoxon reference](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/wilcox.test.html).
* **Welch t-test:** two independent groups, equality of means with unequal variances allowed. Small samples need suitable distributional conditions. The wrapper does not request a paired test. See [R's t-test reference](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/t.test.html).
* **Kruskal-Wallis:** omnibus rank comparison across groups; a significant result does not say which pairs differ. See [R's Kruskal-Wallis reference](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/kruskal.test.html).
* **One-way ANOVA:** omnibus mean comparison using independent errors, approximately normal residuals and common variance. The wrapper does not perform post-hoc tests or model diagnostics. See [R's aov reference](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/aov.html).

For many features, unadjusted significance can generate many false positives. `BH` and `fdr` both select Benjamini-Hochberg false discovery rate adjustment; `bonferroni` selects family-wise adjustment. `none` is the wrapper default. See [R's adjustment reference](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html). A selected feature is not an estimated effect size: this function exports only p-values and feature values. No paired or repeated-measures design is modeled.

## Bray-Curtis, PCoA and PERMANOVA

Bray-Curtis compares nonnegative sample abundance profiles; it ignores joint absence and depends on input normalization. `MBR_beta` leaves data unchanged before computing distances. See [vegan's distance reference](https://vegandevs.github.io/vegan/reference/vegdist.html).

PCoA embeds distances into coordinate axes. Bray-Curtis can generate negative eigenvalues; the wrapper uses `cmdscale` without correction and divides by the sum of all eigenvalues for axis percentages. Treat those percentages accordingly. See [R's cmdscale reference](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cmdscale.html).

PERMANOVA tests a group term in the full distance matrix using permutations. Differences in within-group dispersion can also influence significance, and unrestricted permutations require exchangeable samples. Consider dispersion assessment when interpreting results. See [vegan's adonis2 reference](https://vegandevs.github.io/vegan/reference/adonis.html). This wrapper does not accept restrictions for paired, longitudinal or blocked data. Coordinate-wise plot annotations ask separate questions and have no explicitly configured multiplicity correction.

## Feature selection, effect sizes and predictive evaluation

`MBR_fs` uses random forest RFE with cross-validation to choose a subset. Selected feature effects are then computed on the same samples using `log10(value + 1)` and Cohen's d. Selection can bias subsequent interpretation; use an outer validation design if predictive performance is the objective. The operation and subset-size control follow [caret's RFE documentation](https://topepo.github.io/caret/recursive-feature-elimination.html). Cohen's d is a standardized mean difference; its sign depends on group/reference ordering, and it is not a probability or a fold change. See [rstatix's effect-size reference](https://rpkgs.datanovia.com/rstatix/reference/cohens_d.html).

`MBR_ml` and `MBR_conf` tune through caret and use held-out probabilities for the chosen tuning parameters, averaged across repeats once per sample. `reference_level` consistently identifies the positive class. These resampling diagnostics can remain optimistic after tuning or prior feature selection. For independent evaluation, supply `test_data` and `test_meta_data` containing samples excluded from feature selection, preprocessing estimation and tuning. Test features must match training features; overlapping sample IDs are rejected. The fitted model and diagnostic results are returned invisibly. Confidence intervals on resampling ROC curves are descriptive and omit selection uncertainty. See [caret's training reference](https://topepo.github.io/caret/model-training-and-tuning.html).

## Mantel and correlation summaries

The Mantel test compares sample-distance patterns; its permutation test requires exchangeable sample labels. `MBR_mantel` delegates distance and permutation choices to linkET defaults. The [linkET implementation](https://github.com/Hy4m/linkET/blob/master/R/mantel-test.R) and installed help determine the actual choices for a dependency version. Unscaled metadata units, categorical encodings and confounding can change the meaning of the resulting distances. Mantel associations and the Pearson feature heatmap are different quantities. Color/line categories are descriptive bins, and the wrapper supplies no explicit multiplicity correction.

## Descriptive figures

Circular group means show averages; violins show smoothed distributions; heatmaps show codebook channel values; flow hexagons show counts in two-dimensional bins. These are not interchangeable measures of abundance or uncertainty. A violin annotation imports an existing p-value and does not perform its own test. Scaled codebook colors express standardized values rather than absolute signal intensity.
