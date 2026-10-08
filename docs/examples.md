# Run examples and check documentation

## Install the documented checkout

From the repository root, use a temporary library to keep your normal R library unchanged:

```sh
mkdir -p /tmp/microbiotr-doc-lib
R CMD INSTALL --library=/tmp/microbiotr-doc-lib .
```

This assumes the existing dependency requirements are installed, including the randomForest backend and Bioconductor dependencies. It installs the current source and help files; it does not publish anything to GitHub.

In R:

```r
.libPaths(c("/tmp/microbiotr-doc-lib", .libPaths()))
library(MicroBiotR)
?MBR_read
?MBR_stat
?compute_som
?Gating_workflow
example("count_observations", package = "MicroBiotR")
example("get_counts", package = "MicroBiotR")
example("MBR_reclustering", package = "MicroBiotR")
```

Every package function has an example with `stopifnot` checks where useful. Examples use `tempfile` output directories and synthetic values, not biological evidence. They may assign documented global objects, print figures, change RNG state and leave temporary output files. Execute them in a disposable R session if you want isolation.

## Example categories

* Ordinary examples use small generated objects, with no input-file prerequisite. Run with `example("MBR_stat", package = "MicroBiotR")` and inspect exported CSV files and global result tables.
* `\\donttest` examples use model training, FCS constructors or plot backends. Explicitly request them with `example("MBR_ml", package = "MicroBiotR", run.donttest = TRUE)`; this requires the relevant dependencies.
* `\\dontrun` examples require real FCS files, large gated data or manual interaction. Read and adapt the code before invoking it; placeholder objects and paths are intentionally not executable as-is. The complete SOM pipeline has fixed sample-size and mapping assumptions.

The [Micro.R](Micro.md), [dependent.R](dependent.md) and [Gating](Gating.md) guides repeat function-specific examples and checks. For FCS export, reread output with transformation disabled and compare event counts, channel names and values against the selected input. A file merely existing is not enough to establish a scientifically correct export.

## Documentation validation

Regenerate only R help (then inspect any incidental DESCRIPTION version metadata changes):

```r
roxygen2::roxygenize(".", roclets = "rd")
for (f in list.files("man", pattern = "\\.Rd$", full.names = TRUE)) {
  rd <- tools::parse_Rd(f)
  tools::checkRd(rd)
}
```

`Gating-workflow.Rd` is a manual topic; it is not generated from `R/` because its functions belong to a separate script. Do not run the full Gating script as a documentation check.

A source build verifies packaging; a full check also evaluates dependencies and existing code issues:

```sh
R CMD build --no-build-vignettes .
R CMD check --no-manual MicroBiotR_1.1.1.tar.gz
```

Keep build/check outputs outside the repository when possible. Inspect remaining package-check warnings described in the validation results.

## Regression checks

Run `Rscript tests/regression.R` after installation. These checks cover the requested behavior fixes with synthetic inputs; see [validation results](validation-results.md) for test scope and remaining limitations.
