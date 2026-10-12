# Regression checks

`regression.R` runs with an installed MicroBiotR package and its declared dependencies. Run `Rscript tests/regression.R`, or use `R CMD check` on the built package. Synthetic inputs test statistical results, metadata matching, matrix support, excluded SOM samples, held-out probabilities and independent test evaluation. SOM filtering is tested with a mocked training helper so the test does not train a full map or alter the intentionally retained SOM-resolution behavior.

Run `Rscript tests/cluster-labels.R` from the repository root to check current source code for case-independent cluster matching, ambiguity/missing-ID errors, cross-table p-values, circle labels, and flow/export selection.
