# MicroBiotR branch review and local testing

> Historical snapshot from 8 October 2026. Its test counts and limitations describe that revision. For current results and review instructions, see [the tutorial test report](tutorial-local-test-report.md) and [the runnable notebook](../tutorials/MicroBiotR-tutorial.Rmd).

Review branch: [codex/package-development](https://github.com/changlabs/MicroBiotR/tree/codex/package-development) <br>
Changes: [compare with main](https://github.com/changlabs/MicroBiotR/compare/main...codex/package-development) <br>
Reviewed commit: [`f75ea5f`](https://github.com/changlabs/MicroBiotR/commit/f75ea5f) <br>
Validation date: 8 October 2026

These changes are on the review branch; they have not been merged into `main`.

## Quick route: RStudio with dependencies already installed

Start a fresh RStudio session and run this in the Console. Cloning the repository is unnecessary for running the tutorial.

```r
# Install the review branch without installing or upgrading dependencies.
devtools::install_github(
  "changlabs/MicroBiotR",
  ref = "codex/package-development",
  dependencies = FALSE,
  upgrade = "never"
)
# Load the installed review package.
library(MicroBiotR)
# Open help for the statistical and classification functions.
?MBR_stat
?MBR_ml
```

This assumes devtools is installed as well. [GitHub installer options](https://remotes.r-lib.org/reference/install_github.html) document branch selection and dependency handling. Installation replaces MicroBiotR in the active library; it does not reinstall the declared dependencies.

To run the tutorial:

1. In RStudio, create a new project in an empty directory using **File → New Project → New Directory → New Project**. Use it to hold the example data and generated outputs.
2. Download the matching data/metadata pair from [Zenodo](https://doi.org/10.5281/zenodo.20286378) into that project directory: `data_gated.save` + `meta.txt` for cytometry, `bac_data.txt` + `bac_meta.txt` for 16S, or `rnaibd.csv` + `meta_rna.txt` for RNA.
3. Open the [review-branch tutorial](https://github.com/changlabs/MicroBiotR/blob/codex/package-development/README.md), skip Installation, and run the code chunks in the Console or an R script. Run cytometry steps in order; the 16S and RNA sections are separate workflows. If starting those separately, first load the packages listed at the beginning of the tutorial, then run the preprocessing in that section.
4. Inspect plots, generated files and function help. Cluster IDs, feature counts and sample selections are study-specific and should be reviewed before adapting the tutorial to another dataset.

The package installation does not download the example dataset or the repository's regression scripts. For the regression suite and package checks, use the checkout route below.

## Download and install

Clone the review branch from a terminal:

```bash
# Download the review branch.
git clone --branch codex/package-development https://github.com/changlabs/MicroBiotR.git
# Enter the downloaded repository.
cd MicroBiotR
# Check which revision you downloaded.
git log -1 --oneline
```

Alternatively, open the branch link, select **Code → Download ZIP**, and extract it. ZIP downloads do not include Git history. GitHub access is required if the repository is private.

Open the downloaded folder in RStudio and set it as the working directory. In a fresh R session, run:

```r
# Install installation tools if they are missing.
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# Install the cytometry and metadata dependencies from Bioconductor.
BiocManager::install(c("Biobase", "flowCore"), ask = FALSE, update = FALSE)
# Install the dependency hosted on GitHub.
devtools::install_github("Hy4m/linkET", upgrade = "never")
# Install this local review checkout and its required dependencies.
devtools::install(".", dependencies = NA, upgrade = "never")
# Load the package and inspect its updated help.
library(MicroBiotR)
?MBR_stat
?MBR_ml
```

Installing uses the active R library and can replace an existing MicroBiotR installation there. Use a separate R library if you need to retain another installation. Development tools may be needed when dependencies compile from source. Do not use the README's unqualified GitHub installation command for this review: it installs the default branch rather than this checkout.

## Test locally

From a terminal in the downloaded repository, using the same R installation/library as above:

```bash
# Run the synthetic regression suite against the installed review package.
Rscript tests/regression.R
# Build the source package.
R CMD build .
# Check installation, help, examples and regression tests.
R CMD check --no-manual --no-build-vignettes MicroBiotR_1.1.1.tar.gz
```

The regression suite should finish with **Regression checks passed.** The build creates a source archive, and the check creates a `MicroBiotR.Rcheck` directory. Both are local outputs.

For the tutorial, download the relevant data and metadata pair from [Zenodo](https://doi.org/10.5281/zenodo.20286378), place them in your R working directory, and follow the [branch README](https://github.com/changlabs/MicroBiotR/blob/codex/package-development/README.md). The gated cytometry dataset is about 539 MB; 16S and RNA examples use separate files. Run analyses in a fresh output directory because several functions reuse output filenames.

## Changes made

Documentation now covers the package R scripts and standalone Gating workflow, including statistical assumptions, function help, examples and folder guides. The README has commented tutorial code, dataset/reference links, explanatory figure captions and a 2×2 workflow figure with tutorial step numbers.

The issue numbers below match the earlier review list.

| Issue | Problem | Fix |
| --- | --- | --- |
| 1 | Test warnings/errors became p = 1, hiding failures and valid tied-data results. | Wilcoxon uses a normal approximation; failed/non-finite tests produce a feature-specific warning and NA. Invalid test names stop, and missing results are excluded from significance selection. |
| 2 | Metadata alignment depended only on row order/count. | Sample IDs match and reorder metadata. Mismatched IDs stop; missing usable IDs trigger a positional-alignment warning. |
| 4 | ROC probabilities and confusion metrics used inconsistent positive classes; plots used final-model OOB results. | Both use `reference_level` as positive. Diagnostics use held-out CV probabilities, averaged once per sample across repeats. Optional `test_data` and `test_meta_data` support independent evaluation; fitted models and diagnostics are returned invisibly. |
| 5 | Required dependencies were undeclared or depended on attached packages. | Declared pROC, viridis, tidyr, randomForest and Biobase; qualified relevant dependency accessors and removed unnecessary package-attachment calls. |
| 6 | Excluded low-event samples influenced SOM training. | Samples below 200000 events are removed before pooling/training; all-excluded inputs stop. |
| 8 | The rare-cluster diagnostic counted sample rows. | It now counts feature columns below 0.01% in every retained sample. |
| 9 | Tutorial processing/export indices and sequencing metadata references were inconsistent. | Processing/export consistently use sample 1, export uses acquisition-scale events with original channel names, and sequencing analyses use matched metadata copies. |
| 10 | Some downstream functions accepted matrices but processed individual elements rather than feature columns. | Shared validation converts numeric matrices to data frames and rejects nonnumeric/non-finite feature inputs. |

## Remaining issues and limitations

Issues **3 and 7 were deliberately excluded** from the repairs:

| Issue | Remaining concern |
| --- | --- |
| 3 | SOM mapping/count levels remain fixed at 2025 even when `m` changes; `n_hclust` remains unused. Use the default grid for current testing. |
| 7 | Complex FCS metadata and parameter-key reconstruction still need real-file round-trip validation. Namespace access was repaired under issue 5, but the export algorithm was not changed. |

Other limitations to consider:

- CV diagnostics can remain optimistic after tuning or feature selection. Independent test samples must be excluded from feature selection, estimated preprocessing and tuning. Shared sample IDs are rejected, but repeated subjects with different IDs cannot be detected automatically.
- Feature-selection effect estimates and intervals use the same samples used for selection; they remain exploratory. RFE also pauses between displayed plots.
- Cluster naming differs across steps (`v`, `V`, and numeric prototype row names). Small subsets can lose dimensions in some helpers, and count aggregation can omit out-of-range clusters or samples without a requested resolution.
- Bray–Curtis PCoA can have negative eigenvalues; the implementation does not correct them. PERMANOVA and feature tests do not implement paired/blocked designs through these interfaces.
- Standalone Gating uses instrument-specific thresholds, requires at least six files, writes to `/R_out`, and depends on a private openCyto helper. Source it separately only after reviewing its paths and settings.
- A Mantel plotting fixture without marked correlations failed in a linkET layer; the correlated fixture renders. This edge case remains unresolved.
- Several functions assign objects globally, change RNG state or reuse output filenames.

## What passed and what needs checking again

On R 4.5.2, installation, source build, help syntax checks and regression tests passed. The full package check completed with **0 errors, 2 warnings and 3 notes**. The regression suite covers statistical failures/ties, IDs and matrices, real RFE/classifier execution, independent evaluation, mocked SOM filtering/rarity, and simple FCS export.

The remaining package-check findings are:

| Type | Finding |
| --- | --- |
| Warning 1 | Non-ASCII symbols in executable console messages in `R/Micro.R`. |
| Warning 2 | Gating helper usages are documented, but those functions are not installed package functions. |
| Note 1 | DESCRIPTION text should be a complete sentence. |
| Note 2 | Static analysis treats the randomForest backend as an unused import because caret invokes it dynamically. |
| Note 3 | Existing global plotting variables, global assignments and unqualified base/recommended-package functions need declarations or namespace qualification. |

Please check the following before merging:

1. **Fresh installation:** repeat installation on your OS and R version; verify GitHub/Bioconductor dependency availability. The earlier check used existing local dependencies, and remote repository-index queries were unavailable in its sandbox.
2. **Real tutorial data:** run cytometry, 16S and RNA examples; check sample IDs, retained samples, numerical outputs and plots. Existing README images are illustrations and were not regenerated from a full real-data run after the fixes.
3. **Full SOM:** run the default-grid pipeline on real gated events. Filtering and rarity were tested with mocked helpers, not a complete large-event fit. Confirm count totals, row normalization and sample identities.
4. **Real FCS round trips:** reread exports and compare selected event counts, channel names/values, parameter ranges and metadata with the matching source file. The passing synthetic export does not establish complex instrument compatibility.
5. **Predictive validation:** use an independent subject-level test set and ensure selection/preprocessing use training samples only. Verify probabilities, ROC direction and confusion counts for your chosen positive class.
6. **Statistical suitability:** confirm independence/replication, missing-value handling, normalization, multiple-testing choices and group dispersion where relevant. Passing software tests does not establish biological validity.
7. **Edge cases:** check tiny/single-column subsets, cluster-label conversions, excluded samples, and Mantel plots with weak or absent correlations. Review Gating paths/thresholds and verify its PDFs contain the expected plots.

Please report your OS, `sessionInfo()`, tested commit, commands, and any errors or unexpected outputs. Leave comments on a pull request or share them with the review group; no changes have been merged into `main`.
