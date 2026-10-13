# MicroBiotR development validation report

**Date:** 13 October 2026 <br>
**Package:** MicroBiotR 1.1.1 <br>
**Branch:** `codex/package-development` <br>
**Starting revision:** `2cc4477` (cluster-label consistency fixes)

## Outcome

The local development package installs and builds successfully. `R CMD check --no-manual --no-build-vignettes` reports **0 errors, 0 warnings, 0 notes — Status: OK**. All three regression suites pass. The revised README analyses and the new downloadable R Notebook run successfully on the public example datasets. GitHub publication remains pending the repository owner's permission.

## Repairs included

| Finding | Repair and verification |
| --- | --- |
| Heatmap requested 184 selected features when RFE retained 13 | Plot all `MBR_selected_features`; the revised example produces the 13-feature heatmap. |
| RNA metadata became a vector | Preserve data-frame shape with `drop = FALSE`; the RNA statistics and ROC complete on 98 samples. The same defensive subsetting is used for 16S metadata. |
| Hundreds of undefined p-value warnings | Skip constant features explicitly, keep `NA` p-values, and issue one count summary. Other failures still produce feature-specific warnings. Statistical output tables are unchanged. |
| Unsupported `linewidth` in the beta plot | Use ggscatter's supported point `size` argument. |
| Deprecated flow plotting and FCS metadata access | Replace `aes_string()` with modern data-mask mappings and use `flowCore::keyword()`. |
| Plot limits removed valid events/bins | Crop with Cartesian coordinates while retaining finite positive observations in density calculations. |
| Incorrect SOM cluster-count message | Report 2,025 mapping rows rather than the mapping table's one column. Internal stage labels and the rare-cluster message now describe their actual operations. |
| Package-check warnings and notes | Use portable Unicode escapes; correct standalone-gating help sections; expand DESCRIPTION; add namespace imports and intentional data-mask bindings; assign selected features explicitly to the global environment. |

## Repository documentation and link audit

The review covered R source and comments, installed help, README tutorial text, folder guides, standalone Gating explanations, regression instructions, historical review reports, and figure references. Corrections include outdated uppercase-only matching claims, numeric default codebook claims, manual annotation-order requirements, the incorrect 70-file prerequisite for sample 1, and references to feature-selection superassignment. Repetitive generated comments and duplicate installation commands were removed. Original figures remain clearly labeled as historical illustrations, rather than newly generated validation results.

Older review reports are explicitly labeled as dated snapshots and link to this current report. Their original revision-specific test counts remain historical evidence. All local Markdown/notebook file links resolve. **All 40 existing public documentation and dataset URLs returned HTTP 200 after redirects.** The SVG namespace and Git transport URL were excluded from the HTTP page check. See [the link audit](test-results/link-check.json). Successful URL checks establish availability, not that third-party documentation will remain unchanged.

## Downloadable tutorial

The [R Notebook source](../tutorials/MicroBiotR-tutorial.Rmd) contains all 15 README analysis steps. The [notebook guide](../tutorials/README.md) describes workflow selection and rendering. The [distributed HTML example](../tutorials/MicroBiotR-tutorial.html) contains actual 16S outputs and embedded plots.

The notebook supports `cytometry`, `16s`, `rna`, and `all` modes, separate output folders, optional Zenodo downloads with published MD5 verification, cached SOM reuse, and code-only rendering. It declares both `html_notebook` and `html_document` formats. Its Flatly theme, Tango highlighting, floating three-level contents, numbered sections, paged data frames, 16px prose, 14px code, and slate/orange/red/green alert styling match the supplied `1_data_integrity_check.Rmd` reference. The CSS is embedded in the notebook, with no reference-project stylesheet dependency. It never installs packages automatically. README links explain downloading, RStudio execution, and HTML rendering.

| Notebook validation | Result |
| --- | --- |
| All R chunks parse; all 15 steps are present | Pass |
| Default 16S execution and standalone HTML | Pass |
| RNA execution and standalone HTML | Pass |
| Cytometry downstream execution with validated cached SOM outputs | Pass |
| Combined execution with cached SOM outputs | Pass |
| Notebook preview and HTML with `execute = FALSE` | Pass |
| P-value tables from every executed workflow versus original baselines | Identical |
| HTML figures embedded rather than dependent on external image files | Pass |
| Browser inspection of distributed HTML | Four 16S plots loaded; no horizontal page overflow; computed typography and colors match the reference |
| Missing/corrupt example input | Rejected with a clear error; corrupt files are not replaced |
| Optional download of the two small 16S files | Pass; published MD5 checks match and timeout settings are restored |

Full SOM training was executed in the initial real-data test. Notebook cytometry/all rendering reused those outputs because the fitting and mapping algorithms were unchanged. This is not a claim that training was repeated inside the notebook.

## Real-data test scope

All six example files from [Zenodo record 20286378](https://zenodo.org/records/20286378) matched their published MD5 checksums.

| Workflow | Observed results |
| --- | --- |
| Cytometry | 88 retained samples; 299,992 training events; 2,025 × 7 codebook; 88 × 2,025 normalized count table; row sums equal 1; 88 mapped FCS files; 199 significant clusters; 13 RFE-selected clusters. Initial full SOM workflow took approximately 14 minutes 33 seconds. |
| 16S | 88 samples; 930 input taxa; 125 significant taxa; 42 RFE-selected taxa; ROC output generated. |
| RNA | 98 samples; 58,524 input genes; 3,839 significant genes; ROC output generated. |
| Selected-population export | All 31 events from V9/V10 exactly match source values and event order across the seven acquisition channels after rereading the exported FCS. |

## Tests and environment

- `tests/regression.R`: statistical ties/failures, sample IDs, matrices, SOM exclusion/diagnostics, real RFE/classifiers, held-out probabilities, independent evaluation, and FCS export.
- `tests/cluster-labels.R`: prefix-case combinations, missing/ambiguous labels, flow annotations, and export values/order.
- `tests/tutorial-warnings.R`: all five statistical test names with constant features, channel names containing spaces, preservation of out-of-view density observations, and export without deprecated metadata access.
- Help/example checks, package installation/build, and `git diff --check` pass.

Test environment: R 4.5.2 on Apple silicon macOS; ggplot2 4.0.3, caret 7.0.1, linkET 0.1.0, flowCore 2.20.0, Biobase 2.73.1, rmarkdown 2.31, knitr 1.51, and RStudio's bundled Pandoc. Dependencies were already installed. Installation from scratch on another machine and remote GitHub installation remain to be verified by collaborators. The [package-check log](test-results/R-CMD-check.log) is included with machine-specific directory paths replaced by placeholders.

## Remaining limitations

- linkET still emits upstream ggplot2 deprecations for `aes_()` and the line `size` aesthetic. Mantel plots complete; warnings are not suppressed.
- One constant-feature summary remains expected per 16S/RNA statistical call (406 taxa; 13,252 genes). Their p-values remain `NA`.
- Missing/nonpositive channel values remain invalid for logarithmic axes. The plot repair does not conceal invalid input.
- SOM count levels remain fixed at 2,025 and `n_hclust` remains unused. Keep the documented default resolution.
- Standalone raw-data gating, complex instrument-specific FCS metadata, and biological/predictive validity are outside these computational checks. See [implementation limitations](implementation-notes.md).
- Original README plots were not regenerated. The new downloadable HTML contains fresh 16S plots.

## Reproduce after downloading the review checkout

Install the local checkout after installing its dependencies, then run:

```sh
R CMD INSTALL .
Rscript tests/regression.R
Rscript tests/cluster-labels.R
Rscript tests/tutorial-warnings.R
R CMD build --no-build-vignettes .
R CMD check --no-manual --no-build-vignettes MicroBiotR_1.1.1.tar.gz
```

For a fresh HTML tutorial run, download the notebook, install the tutorial dependencies, and run from its folder:

```r
rmarkdown::render("MicroBiotR-tutorial.Rmd", output_format = "html_document",
                  params = list(workflow = "16s", download_data = TRUE,
                                output_dir = "results-review"))
```

Record the tested revision, `sessionInfo()`, commands, and unexpected results when sharing feedback. Package check does not render the notebook; notebook validation is a separate step. Tutorial data/results and local IDE histories are excluded from the prepared Git changes.
