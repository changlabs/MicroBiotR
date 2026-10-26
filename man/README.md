# Installed R help pages

The `.Rd` files provide help shown by `?<function_name>` after installing MicroBiotR. They include signatures, arguments, actual return behavior, side effects, methodological interpretation and examples.

Package-function pages are generated from roxygen comments in [`R/`](../R/README.md). `Gating-workflow.Rd` documents the standalone workflow and its six helper functions; those helpers are not installed APIs. Their help aliases can be looked up, but the functions exist only after evaluating the appropriate definitions from `Gating` in an analysis session. Do not source the full workflow merely to load helpers: it immediately reads data and writes outputs.

`\dontrun{}` examples need real files, large data, or manual interaction; `\donttest{}` examples are runnable but involve additional packages or model fits. Ordinary examples are small synthetic demonstrations. See [validation instructions](../docs/examples.md).

Regenerate package help from the repository root using `roxygen2::roxygenize(".", roclets = "rd")`. Preserve the manual Gating topic. Installing this checkout is necessary to see new help; an older installed package continues to show its old documentation.
