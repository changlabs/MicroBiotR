library(MicroBiotR)
plot_file <- tempfile(fileext = ".pdf")
grDevices::pdf(plot_file)
inputs <- getFromNamespace('.mbr_inputs', 'MicroBiotR')
predictions <- getFromNamespace('.mbr_predictions', 'MicroBiotR')
expect_error <- function(expr) stopifnot(inherits(tryCatch(expr, error = identity), 'error'))
x <- data.frame(v1 = c(0, 0, 1, 1, 4, 4, 5, 5), v2 = rep(1, 8),
                row.names = paste0('s', 1:8))
meta <- data.frame(Group = rep(c('A', 'B'), each = 4), row.names = rownames(x))
stopifnot(identical(inputs(as.matrix(x), meta[8:1, , drop = FALSE], 'Group')$meta_data, meta))
bad <- meta; rownames(bad)[1] <- 'wrong'
expect_error(inputs(x, bad, 'Group'))
expect_error(inputs(transform(x, v1 = as.character(v1)), meta, 'Group'))
expect_error(inputs(x, meta, 'absent'))
folder <- tempfile(); dir.create(folder)
suppressWarnings(MBR_stat(as.matrix(x), meta[8:1, , drop = FALSE], 'Group', out_path = folder))
expected <- wilcox.test(x$v1 ~ meta$Group, exact = FALSE)$p.value
stopifnot(isTRUE(all.equal(pvalue_data['v1', 'p.value'], expected)),
          expected < 1, is.na(pvalue_data['v2', 'p.value']),
          !'v2' %in% names(significant_data))
expect_error(MBR_stat(x, meta, 'Group', test_type = 'invalid', out_path = folder))
# Repeats count each sample once; positive class can be the second alphabetic level.
fake <- list(pred = data.frame(rowIndex = c(1, 1, 2, 2),
                 obs = factor(c('A','A','B','B')), B = c(.1,.3,.8,1)))
p <- predictions(fake, 'B')
stopifnot(nrow(p) == 2, isTRUE(all.equal(p$probability, c(.2,.9))),
          identical(as.character(p$predicted), c('A','B')))
# Test the SOM exclusion before the expensive fit using isolated mocked helpers.
e <- new.env(parent = environment(MBR_som)); e$MBR_som <- MBR_som
environment(e$MBR_som) <- e
seen <- NULL
e$downsample <- function(data, ...) { seen <<- names(data); matrix(1, 4, 2) }
e$compute_som <- function(...) stop('training sentinel')
fcs <- list(keep = list(data = matrix(1, 200000, 2)),
            drop = list(data = matrix(99, 3, 2)))
suppressWarnings(expect_error(e$MBR_som(fcs, out_path = folder)))
stopifnot(identical(seen, 'keep'))
seen <- NULL
suppressWarnings(expect_error(e$MBR_som(fcs['drop'], out_path = folder)))
stopifnot(is.null(seen))
# Exercise the actual SOM diagnostic without fitting or exporting a full map.
e$compute_som <- function(...) list(som = list(codes = list(matrix(1, 2, 2)),
                                              grid = list(pts = matrix(1, 2, 2))))
e$map_som <- function(data, ...) data
e$assign_clusters <- function(data, ...) data
e$count_observations <- function(data, ...) data
e$get_counts <- function(...) list(matrix(c(0, 2, 0, 3), nrow = 2,
                              dimnames = list(c('v1','v2'), c('keep','other'))))
e$flowFrame <- function(...) stop('export sentinel')
messages <- character()
withCallingHandlers(expect_error(e$MBR_som(list(keep = fcs$keep, other = fcs$keep),
                                         out_path = folder)),
                    warning = function(w) {
                      messages <<- c(messages, conditionMessage(w))
                      invokeRestart('muffleWarning')
                    })
stopifnot(any(grepl('1 clusters are below', messages, fixed = TRUE)))
# Real RF resampling and independent evaluation, without attached tidyr/Biobase.
set.seed(16)
features <- data.frame(v1 = c(rnorm(15), rnorm(15, 3)), v2 = rnorm(30),
                       row.names = paste0('train', 1:30))
labels <- data.frame(Group = rep(c('A', 'B'), each = 15), row.names = rownames(features))
warn_before <- getOption('warn')
fit <- MBR_ml(as.matrix(features), labels[30:1, , drop = FALSE], reference_level = 'B',
              method = 'cv', number = 3, out_path = folder)
stopifnot(nrow(fit$predictions) == 30, fit$confusion$positive == 'B',
          identical(getOption('warn'), warn_before), fit$roc$levels[2] == 'B')
test <- features[1:10, ]; rownames(test) <- paste0('test', 1:10)
test_meta <- data.frame(Group = rep(c('A','B'), each = 5), row.names = rownames(test))
conf <- MBR_conf(features, labels, reference_level = 'B', method = 'cv', number = 3,
                 test_data = test, test_meta_data = test_meta, out_path = folder)
stopifnot(nrow(conf$predictions) == 10, sum(conf$confusion$table) == 10)
# Plotting and RFE accept matrix inputs without attaching tidyr.
MBR_circle(as.matrix(x), meta[8:1, , drop = FALSE], 'Group', out_path = folder)
MBR_violin(as.matrix(x), meta, 'Group', pvalue_data, cluster = 1, out_path = folder)
set.seed(42)
rfe_data <- matrix(runif(120, .01, 1), nrow = 20,
                   dimnames = list(paste0('r',1:20), paste0('v',1:6)))
rfe_meta <- data.frame(Group = rep(c('A','B'), each = 10), row.names = rownames(rfe_data))
MBR_fs(rfe_data, out_path = folder, meta_data = rfe_meta, group_name = 'Group',
       nfolds_cv = 2, rfe_size = 3, top_n_features = 3, ref_group = 'A')
stopifnot(nrow(MBR_selected_features) == 20, !'package:tidyr' %in% search())
# Verify acquisition-scale export works with no Biobase attached to the search path.
ff <- flowCore::flowFrame(matrix(c(1,2,3,4,1,2), nrow = 2,
                        dimnames = list(NULL, c('FSC','SSC','classes'))))
fs <- flowCore::flowSet(list(example = ff))
dat <- MBR_process(fs, transformation = identity)
MBR_save(fs, dat, selected_rows = 'V1', rawdata_path = folder)
# Invalid independent sample reuse stops before fitting.
expect_error(MBR_ml(features, labels, reference_level = 'B', test_data = features,
                    test_meta_data = labels, out_path = folder))
unlink(folder, recursive = TRUE)
grDevices::dev.off()
unlink(plot_file)
cat('Regression checks passed.\n')
