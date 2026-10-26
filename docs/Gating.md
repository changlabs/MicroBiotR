# Standalone `Gating` R workflow

`Gating` is an R script without an `.R` extension. It is outside `R/` and is not executed by package installation. It is an instrument-specific analysis template with top-level commands, not a reusable package function. `source("Gating")` immediately loads libraries, reads FCS files, creates folders and writes output. Review its assumptions in a separate analysis copy before execution.

## Inputs and dependencies

The script attaches flowCore, flowWorkspace, openCyto, ggcyto, ggplot2 and ggthemes. Their requirements are additional to the package's declared imports; PDF hexagonal plotting also needs hexbin. Input defaults to `./` and only names ending in lowercase `.fcs` are selected. FCS reading explicitly disables transformation and preserves channel names, while excluding channels matching `*`, `Bits`, or `Drop` through a regular expression.

The following original channel-to-marker mapping must match the instrument: Hoechst Red → DNA, FITC → hIgA2, APC → hIgA1, Pe-TR → hIgG, and BV650 → hIgM. Scatter channels must be `FSC PAR` and `SSC`. Gate thresholds use raw instrument units, so compensation, acquisition range, gain changes and other instruments require independent assessment. Gates identify operational populations rather than confirming organism identity or viability.

## What each stage does and why

1. Creates `/R_out/data` and `/R_out/graphs` and reads files from the working directory. These absolute output paths need write permission; source only when that location is intended.
2. Sets marker labels and builds a flowWorkspace `GatingSet`, which maintains parent-child populations for each sample.
3. Sets seed 2025 and selects six distinct sample indices for preview. At least six input files are required; a smaller experiment fails here. A preview subset is for visual inspection, not an inferential sample design.
4. Creates a rectangular debris boundary on FSC PAR and SSC, both from 2000 to 63000, and adds `non-debris` below `root`. This removes low-scatter/background events and excludes events near the specified upper boundary. Recomputing applies the gate to every sample.
5. Draws and saves debris gate diagnostics. Inspect the actual boundary relative to event density across all samples rather than assuming six previews establish general validity.
6. Adds the DNA polygon under `non-debris`: FSC PAR 2000–65000 and Hoechst Red 18000–67000. The parent relationship intersects the DNA gate with the debris gate; it does not apply the two independently.
7. Recomputes and writes DNA and full gating figures.
8. Extracts events from the DNA population for seven channels, appends marker labels to selected channel names and retains each sample name. This produces the named sample structure needed by `MBR_som`.
9. Saves the object `fl_data_ig` in `FCM_gated/fl_data_ig.save` using `save`, not `saveRDS`. Use `load` into a separate environment to inspect it without overwriting existing objects.

Outputs are `/R_out/graphs/gating/Ig_gate_debris.pdf`, `Ig_gate_dna.pdf`, `gating_Ig.pdf`, and `/R_out/data/FCM_gated/fl_data_ig.save`. Existing files may be overwritten. The script creates many objects in the environment where it is sourced.

## Helper functions

These six functions belong to the standalone script. Installed help aliases such as `?gate_dna` explain them, but loading MicroBiotR does not define them. Evaluate only the relevant function definitions from a copy if you need isolated testing.

### `gate_debris(data)`

Calls the private `openCyto:::.boundary` helper with channels FSC PAR and SSC and the fixed limits above. Returns its boundary gate result for `gs_pop_add`. It does not infer thresholds from a fitted distribution. This depends on an internal openCyto API, so behavior can change between dependency versions. Validate the gate against raw units and sample-specific background. An isolated example requires a compatible flow object: `gate <- gate_debris(ig_gs)`; inspect `gate` before adding it.

### `plot_gate_debr(data, gate, path = tmp_path)`

Creates a scatter-channel hexagonal ggcyto display for the named population, adding its gate and population statistics. `data` is a flowWorkspace object accepted by ggcyto; `gate` is the subset name, such as `"root"`. `path` is unused, despite its global default. Returns the plot expression; on a caught error it prints a message and returns an empty ggplot titled `Error`. It does not write files itself. Example: `p <- plot_gate_debr(ig_gs[[1]], "root"); print(p)` after preparing compatible data. Gate statistics describe event fractions, not hypothesis-test results. The script calls `print(lapply(...))`; printing a list of plot objects is not equivalent to explicitly printing each plot, so inspect PDFs for populated pages.

### `gate_dna(data, channels = c("FSC PAR", "Hoechst Red"))`

Returns a flowCore polygonGate with four fixed corners. `data` is unused: thresholds are identical regardless of supplied sample. `channels` assigns the two coordinate names in order. A small isolated check needs no experiment once the function definition is available:

```r
g <- gate_dna(NULL)
stopifnot(inherits(g, "polygonGate"))
```

This checks the object type, not biological validity. Inspect the rectangle and event distribution before adding it under non-debris.

### `plot_gate_dna(data, gate, path = tmp_path)`

Plots FSC PAR against Hoechst Red for a named population with hexagon counts, gate and descriptive statistics. Input, unused path, error fallback and return conventions match `plot_gate_debr`. Example: `p <- plot_gate_dna(ig_gs[[1]], "non-debris"); print(p)`. PDF printing in the surrounding script uses the same list-printing pattern.

### `plot_gating(x, nodes = gs_get_pop_paths(ig_gs))`

Returns `autoplot(x, bins = 64, nrow = 1, ylim = "instrument", xlim = "instrument")`. `x` is a compatible gating hierarchy. `nodes` is unused; it neither selects nor filters populations. The default references the global ig_gs, but is normally not evaluated because the argument is unused. Example: `p <- plot_gating(ig_gs[[1]]); print(p)`. Again verify the enclosing PDF contains all intended panels.

### `get_sc_gated(data, node, channels = NULL)`

Extracts the sample name using pData, obtains events for the population through gh_pop_get_data and exprs, optionally subsets channels and appends `.marker` to names where marker labels exist. `data` is one gating hierarchy, `node` is a population such as `"DNA"`, and `channels` is a channel-name vector. Returns `list(name = sample_name, data = event_matrix)`. With channels = NULL all columns are retained but marker matching yields no selected labels. The subset lacks drop = FALSE, so use multiple channels and verify the matrix shape for sparse gates. The calling script selects seven channels.

```r
# After preparing ig_gs and evaluating this helper's definition:
x <- get_sc_gated(ig_gs[[1]], node = "DNA", channels = c("FSC PAR", "SSC"))
stopifnot(is.list(x), is.matrix(x$data), ncol(x$data) == 2,
          all(is.finite(x$data)))
```

## Check the saved output without rerunning gates

```r
e <- new.env(parent = emptyenv())
load("/R_out/data/FCM_gated/fl_data_ig.save", envir = e)
gated <- e$fl_data_ig
stopifnot(is.list(gated), length(gated) > 0,
          all(vapply(gated, function(x) is.matrix(x$data), logical(1))))
channels <- lapply(gated, function(x) colnames(x$data))
stopifnot(all(vapply(channels, identical, logical(1), channels[[1]])))
```

Also inspect sample names, event counts, retained fractions and channel distributions. The SOM pipeline later removes samples with fewer than 200000 events. No real FCS files are provided, so full gating examples require your own experiment and verified dependencies.
