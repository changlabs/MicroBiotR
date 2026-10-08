# Observed implementation limitations

The current implementation addresses the requested statistical, alignment, dependency and tutorial issues. The remaining limitations below should be considered when choosing an analysis.

| Area | Existing behavior and practical consequence |
| --- | --- |
| SOM resolution | Grid size follows `m`, but mapping/count levels remain fixed at 2025; `n_hclust` is unused. Keep the default grid when reproducing this workflow and inspect assignments. |
| Cluster names | Counts use `v` prefixes; flow selection uses `V`; heatmap codebook lookup forces `V`; default SOM codebook rows may be numeric. Convert and align copies explicitly. |
| Empty/small subsets | Several helpers subset without `drop = FALSE`; a single row or channel may lose matrix dimensions. |
| Count completeness | `count_observations` drops IDs outside its supplied range; `get_counts` omits samples lacking a resolution. Check sums and column names. |
| PCoA | Three requested axes require adequate distinct samples; negative Bray-Curtis eigenvalues are not corrected and reported percentages include all eigenvalues in their denominator. |
| Feature selection | RFE and subsequent effect exploration use the same dataset; effect intervals are not adjusted for selection. The function pauses with `readline` after plots. |
| FCS export | Channel matching expects original names, and the exported values are the supplied processed values, not reconstructed raw values. Metadata and parameter-key reconstruction warrant round-trip inspection on real files. |
| Gating workflow | Fixed raw-unit gates, at least six files, absolute `/R_out` outputs and private openCyto `.boundary` dependency. Plot lists are printed rather than each plot explicitly; check that PDFs contain intended figures. Several helper arguments are unused. |
| Mantel rendering | A random fixture with no marked feature correlations failed in a linkET plot layer locally; the documented correlated fixture renders. Investigate this edge case separately. |
| Output folders | Most plotting/testing functions require an existing out_path. SOM/classifier/export functions create their own documented directories. Common filenames can be overwritten. |
| Session effects | SOM, tests, RFE and reclustering assign output objects globally. Feature selection and MBR_ml reset the seed to 2025. circlize state is reset by MBR_circle. |

## Requested issue fixes

| Issue | Behavior |
| --- | --- |
| 1 | Invalid test names stop; tied Wilcoxon data use a normal approximation; failed tests retain NA with warnings. |
| 2 | Sample IDs align metadata to feature rows; mismatched IDs stop. Missing IDs produce a positional-alignment warning. |
| 4 | ROC/confusion use the same positive class and held-out predictions. Optional independent test tables bypass resampling diagnostics; tuning/selection must exclude those samples. |
| 5 | pROC, viridis, tidyr, randomForest and Biobase are declared; relevant accessors use explicit namespaces. |
| 6 | Low-event samples are removed before pooling and SOM fitting; an empty retained set stops. |
| 8 | Rare clusters are counted across feature columns. |
| 9 | Tutorial processing/export use sample 1 consistently; sequencing analyses use aligned metadata copies. |
| 10 | Shared validation converts numeric matrices to data frames and rejects nonnumeric or non-finite features. |

Issues 3 (fixed SOM count levels and unused n_hclust) and 7 (complex FCS metadata/export reconstruction) remain outside this repair scope.
