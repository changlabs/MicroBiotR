# Shared analysis validation

`R/validation.R` defines internal helpers used by the sample-level functions. Numeric matrices are converted once to data frames so columns remain features in tests, group averages and named cluster extraction. Nonnumeric columns, non-finite values, empty inputs and invalid grouping columns stop before analysis.

`.mbr_inputs` preserves feature-row order and matches metadata by sample row name. Both ID sets must be unique and equal. Metadata is reordered to the features, avoiding accidental comparisons against shuffled labels. Default sequential row names are treated as unavailable IDs; these inputs retain positional behavior with a warning. Set meaningful row names on both tables to verify correspondence.

`.mbr_predictions` reads caret's saved held-out predictions for the chosen tuning parameters. Repeated cross-validation probabilities are averaged by original sample index, yielding one diagnostic row per sample. A probability of at least 0.5 classifies the sample as `reference_level`; that class is also positive in ROC estimation and confusion metrics. This ensures repeats do not multiply the confusion counts. These diagnostic estimates are conditional on the selected model and features, and are not independent validation of model selection.

`.mbr_test_inputs` checks optional independent test tables before fitting. It rejects missing paired arguments, mismatched features and overlap between meaningful training and test sample IDs. `.mbr_evaluate` uses predictions from the fitted model for independent test samples when supplied; otherwise it uses saved resampling probabilities. Both classes must occur in the test metadata. Matching IDs cannot detect shared subjects with different identifiers or preprocessing leakage, so the experimental split must precede feature selection and estimated transformations.

For an independent evaluation, choose features using training samples only, supply those feature columns in both tables, and call:

```r
# Train and tune on training samples; evaluate on separate test samples.
result <- MBR_ml(train_features, train_metadata, group_name = "Group",
                 reference_level = "CD", test_data = test_features,
                 test_meta_data = test_metadata, out_path = output_directory)
# Inspect per-sample probabilities and confusion metrics.
result$predictions
result$confusion
```

The regression checks in `tests/regression.R` cover tied and failed tests, matrix conversion, reordered and mismatched IDs, exclusion before SOM training, positive-class consistency, repeated predictions, and independent evaluation.
