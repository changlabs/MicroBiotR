# Package R source

R loads the top-level `.R` files in this directory when MicroBiotR is installed. Comments and help explain the analysis behavior; validation.R provides shared input and prediction checks.

| Script | Role | Detailed guide |
| --- | --- | --- |
| `Micro.R` | Public workflow: SOM, statistics, figures, feature selection, classification, FCS processing and export | [Micro.R walkthrough](../docs/Micro.md) |
| `globals.R` | Registers intentional column symbols used by plotting/data masks | [Validation guide](../docs/validation.md) |
| `validation.R` | Shared sample validation and held-out/test prediction checks | [Validation guide](../docs/validation.md) |
| `dependent.R` | Sampling, SOM training and mapping, label lookup, sample counts and count-table assembly | [dependent.R walkthrough](../docs/dependent.md) |

The separate [`Gating`](../Gating) workflow lives outside the package source and is **not executed or exported by the installed package**. See [its guide](../docs/Gating.md) before running it.

## Function help

After installing this checkout, use `library(MicroBiotR)` and `?MBR_stat` (or any other documented function). `?assign_clusters` opens the internal helper's help page; calling it requires `MicroBiotR:::assign_clusters`. Functions defined only inside another function are local closures and cannot be called as package APIs: `stat_func` implements the ANOVA branch of `MBR_beta`, and `show_plots` implements the interactive pause in `MBR_fs`.

Roxygen comments in these scripts are the source for `man/*.Rd`. Regenerate help only with `roxygen2::roxygenize(".", roclets = "rd")` when making documentation-only changes. When changing imports, update the roxygen import declarations and regenerate the namespace with the `namespace` roclet as well. Review generated changes against the intended scope. Current roxygen versions can update DESCRIPTION metadata, so inspect and restore that incidental change when the scope is documentation only.

## Input conventions

* Gated samples: named list of `list(name = "sample.fcs", data = numeric_event_matrix)`; columns must be identical and ordered consistently between samples.
* Count helper output: cluster rows (`v1`, `v2`, ...) and sample columns. Transpose before sample-level downstream analyses.
* Downstream abundances: numeric matrix or data frame with samples in rows and features in columns; sample row names align metadata. Missing IDs produce a positional-alignment warning.
* Cluster matching accepts `v1` and `V1` equivalently in violin, heatmap, bin preparation, flow plotting and export. Input names are preserved; unrelated feature names remain case-sensitive. Missing requests and duplicate or case-ambiguous IDs raise errors. Numeric event classes are supported; codebook rows should have v/V-prefixed cluster labels.

See [statistical interpretation](../docs/statistical-methods.md), [examples and checks](../docs/examples.md), and [known limitations](../docs/implementation-notes.md).
