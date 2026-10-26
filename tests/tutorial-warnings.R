# Run against the locally installed development package.
library(MicroBiotR)
collect_warnings <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart('muffleWarning')
  })
  list(value = value, warnings = warnings)
}
folder <- tempfile(); dir.create(folder)
x <- data.frame(variable = c(0, 0, 1, 2, 3, 4, 4, 5),
                zero = rep(0, 8), constant = rep(3, 8),
                row.names = paste0('s', 1:8))
meta <- data.frame(Group = rep(c('CD', 'HC'), each = 4), row.names = rownames(x))
for (test in c('wilcox', 'kruskal', 'anova', 't.test', 'ttest')) {
  result <- collect_warnings(MBR_stat(x, meta, 'Group', test_type = test, out_path = folder))
  stopifnot(length(result$warnings) == 1L,
            grepl('Skipped 2 constant features', result$warnings),
            all(is.na(pvalue_data[c('zero', 'constant'), 'p.value'])),
            is.finite(pvalue_data['variable', 'p.value']),
            !any(c('zero', 'constant') %in% names(significant_data)))
}
# Test modern mappings with acquisition names containing spaces, while retaining
# out-of-view positive events in hex-density calculations.
dat <- data.frame('FSC PAR' = c(1, 2, 5, 20), SSC = c(1, 2, 5, 20),
                  classes = c(1, 1, 2, 2), check.names = FALSE)
plot <- MBR_flow_plot(dat, x = 'FSC PAR', y = 'SSC', selected_rows = 'V1',
                      x_limits = c(.9, 10), y_limits = c(.9, 10))
result <- collect_warnings(ggplot2::ggplot_build(plot))
stopifnot(length(result$warnings) == 0L,
          sum(result$value$data[[1]]$count) == nrow(dat),
          nrow(result$value$data[[2]]) == 2L)
# Export without deprecated metadata access, and verify acquisition-scale values.
ff <- flowCore::flowFrame(as.matrix(dat))
fs <- flowCore::flowSet(list(example = ff))
result <- collect_warnings(MBR_save(fs, MBR_process(fs, transformation = identity),
                                    selected_rows = 'V1', rawdata_path = folder))
stopifnot(length(result$warnings) == 0L)
export <- flowCore::exprs(flowCore::read.FCS(file.path(folder, 'FilteredFCS/example_filtered.fcs'),
                                           transformation = FALSE))
stopifnot(identical(unname(colnames(export)), c('FSC PAR', 'SSC')),
          isTRUE(all.equal(as.vector(export), as.vector(as.matrix(dat[1:2, 1:2])))))
unlink(folder, recursive = TRUE)
cat('Tutorial warning checks passed.\n')
