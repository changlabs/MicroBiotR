# MicroBiotR cluster-label case consistency report

> Historical snapshot from 12 October 2026. Its test counts and limitations describe that revision. For current results and review instructions, see [the tutorial test report](tutorial-local-test-report.md).

Date: 12 October 2026 <br>
Branch: `codex/package-development` <br>
Implementation commit: [`2cc4477`](https://github.com/changlabs/MicroBiotR/commit/2cc4477)

## Findings

Cluster labels were handled inconsistently across functions. Some lookups expected lowercase labels such as `v1`, while other functions expected uppercase labels such as `V1`. This could cause failed lookups, incorrect or missing p-value annotations, or empty event selections.

The original violin-plot issue was caused by constructing `paste0("v", cluster)` regardless of the input column's spelling. The p-value lookup used that same name, so differences in prefix case between the data and p-value table were also problematic.

Review identified related inconsistencies in heatmaps, bin preparation, flow plots, FCS export and circle-plot display labels.

## Matching policy

The fix compares complete cluster labels using a canonical uppercase prefix. For example, `v1` and `V1` both compare as `V1`. Input tables are not renamed or modified by this normalization.

Only labels matching `^[vV][0-9]+$` are normalized. Unrelated names remain case-sensitive, and leading zeros remain significant. Numeric event classes are supported by a separate event-label helper.

| Case | Example | Result |
| --- | --- | --- |
| Prefix case | `v1` versus `V1` | Equivalent for cluster matching. |
| Mixed case between tables | Data column `V1`, p-value row `v1` | Matched independently and successfully. |
| Missing request | Request `V9` when only `v1` and `V2` exist | Explicit error identifying the missing cluster. |
| Ambiguous names | One table contains both `v1` and `V1` | Explicit error rather than selecting the first match. |
| Repeated event classes | `v1, v1, V2` | Accepted; multiple events in one class are expected. |
| Ambiguous event encoding | Both `v1` and `V1`, or both `1` and `V1` | Explicit ambiguity error; use one spelling or encoding per class. |
| Unrelated names | `VitaminV` versus `vitaminV` | Distinct and case-sensitive. |
| Leading zeros | `V01` versus `V1` | Distinct IDs. |
| Numeric event classes | Classes `1, 2`; requests `v1, V2` | Matches successfully. |

## Helpers added and their locations

Three internal helpers were added to [`R/validation.R`](../R/validation.R), alongside the existing shared input-validation functions. They are not exported as public package APIs. The line numbers below refer to implementation commit `2cc4477`.

### `.mbr_cluster_ids()` — line 104

```r
.mbr_cluster_ids(labels, context = "Cluster labels", unique = TRUE)
```

Converts labels to character and normalizes only complete v/V-prefixed numeric cluster labels. It rejects missing or empty labels. By default, it also rejects duplicate normalized IDs, including ambiguous pairs such as `v1` and `V1`. The `context` argument identifies the relevant input in error messages.

```r
.mbr_cluster_ids(c("v1", "V2", "VitaminV", "v01"))
# c("V1", "V2", "VitaminV", "V01")

.mbr_cluster_ids(c("v1", "V1"))
# Error: duplicate or case-ambiguous IDs: V1
```

`unique = FALSE` permits repeated labels when normalizing individual events. The event helper validates the distinct class labels before using this option.

Used directly by `MBR_circle`, and internally by the other two helpers.

### `.mbr_cluster_match()` — line 121

```r
.mbr_cluster_match(requested, available, context = "Available clusters")
```

Normalizes both requested and available labels, validates their uniqueness, and returns matching integer indices in requested order. Missing requests raise an error instead of creating NA rows or silently empty selections. Available labels retain their original spelling in the source table.

```r
.mbr_cluster_match(c("V2", "v1"), c("V1", "v2"))
# c(2L, 1L)

.mbr_cluster_match("V9", c("v1", "V2"))
# Error: missing requested clusters: V9
```

Used by `MBR_violin`, `MBR_heatmap`, `MBR_prepare`, `MBR_flow_plot` and `MBR_save`. The `MBR_plot` wrapper inherits the behavior through `MBR_flow_plot`.

### `.mbr_event_ids()` — line 136

```r
.mbr_event_ids(classes)
```

Handles numeric event classes, numeric strings and prefixed labels. It adds `V` to wholly numeric strings, validates distinct class encodings for ambiguity, then normalizes the full event vector while retaining repeated occurrences and original event order.

```r
.mbr_event_ids(c(1, 1, 2))
# c("V1", "V1", "V2")

.mbr_event_ids(c("v1", "v1", "V2"))
# c("V1", "V1", "V2")

.mbr_event_ids(c("1", "V1"))
# Error: duplicate or case-ambiguous IDs: V1
```

Used by `MBR_flow_plot` and `MBR_save` for event selection. `MBR_heatmap` also uses it to preserve support for numeric abundance-column names.

## Function changes and examples

Public-function changes are in [`R/Micro.R`](../R/Micro.R).

| Function | Example | Corrected behavior |
| --- | --- | --- |
| `MBR_violin` | Data column `V1`, p-value row `v1`, `cluster = 1` | Selects the correct abundance column and p-value. The reverse combination also works. |
| `MBR_heatmap` | Abundance columns `v2, V1`; codebook rows `v1, V2` | Matches codebook rows in abundance-column order: cluster 2, then cluster 1. |
| `MBR_prepare` | Bin rows `v1, V2, v3`; requested labels `V3, V1` | Retains `v1, v3`, preserving original bin order and names. |
| `MBR_flow_plot` / `MBR_plot` | Event classes `1, 1, 2, 3`; requested labels `v2, V1` | Highlights events in clusters 1 and 2, preserving event order. |
| Flow-plot bin labels | Bin rows `v1, V2`; requested labels `v2, V1` | Aligns each label with its cluster's coordinates, regardless of bin-row order or prefix case. |
| `MBR_save` | Event classes `1, 1, 2`; request `v1` or `V1` | Exports the same two cluster-1 events with their values and order preserved. |
| `MBR_circle` | Feature names `v1, V2, VitaminV` | Displays `Bin1, Bin2, VitaminV`; only complete cluster labels are reformatted. |

Example calls using existing workflow objects:

```r
# abundances has column V1; pvalues has row v1.
MBR_violin(abundances, metadata, "Group", pvalues, cluster = 1,
           out_path = output_directory)

# Prefix case may differ between abundance columns and codebook rows.
MBR_heatmap(abundances, codebook, out_path = output_directory)

# Selection accepts either prefix case and preserves input bin-row order.
selected_bins <- MBR_prepare(bins, c("V3", "v1"), transformation = identity)

# Bin coordinates are matched by cluster ID before annotations are placed.
plot <- MBR_flow_plot(events, bins = bins, x = "FSC", y = "SSC",
                      selected_rows = c("v2", "V1"))

# Either spelling selects the same cluster-1 events for export.
MBR_save(fcs_files, events, selected_rows = "v1", file_index = 1,
         rawdata_path = output_directory)
```

These examples assume the usual valid inputs, metadata, channels, output directories and package dependencies.

## Scope and intentional behavior changes

Beyond accepting either prefix case, the changes deliberately make missing requests and ambiguous IDs errors. Workflows that previously depended on silent omission must supply valid, unambiguous cluster IDs.

Circle-plot display labels now handle lowercase and uppercase cluster IDs consistently. Flow-plot annotations now use cluster-ID matching to align coordinates rather than relying on supplied bin-row order.

The normalization is scoped to the functions listed above. It does not make every feature-name lookup in the package case-insensitive or change statistical calculations, SOM training or classifier algorithms. Heatmap codebook row names still need v/V-prefixed cluster labels; bare numeric codebook row names are not normalized by the cluster-matching helper.

## Documentation and validation

Detailed comments were added to the helpers and their call sites. Function help pages and repository guides were updated to explain the new behavior.

Validation completed successfully:

- Installation of the package into an isolated temporary R library.
- Existing regression suite: [`tests/regression.R`](../tests/regression.R).
- New checks: [`tests/cluster-labels.R`](../tests/cluster-labels.R), run against both current source and the installed package.
- Coverage of all four uppercase/lowercase combinations between violin data and p-value rows, mixed-case heatmap matching, bin selection, missing and ambiguous IDs, unrelated names, event repetitions and circle labels.
- Inspection of actual flow-plot layers to verify event selection and bin-coordinate alignment.
- FCS export and rereading to verify event values, channel names and event order for both requested prefix cases and numeric/prefixed classes.
- `git diff --check`.

The FCS checks emitted existing flowCore deprecation warnings. A full `R CMD check` was not run for this change.

## Delivery

The implementation, tests and supporting documentation were pushed to `codex/package-development` in [commit `2cc4477`](https://github.com/changlabs/MicroBiotR/commit/2cc4477). This report was created afterward and is not part of that commit. No merge to the main branch was performed as part of this work.
