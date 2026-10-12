# Run from the repository root: Rscript tests/cluster-labels.R
# Source current code so these checks do not accidentally test an older install.
suppressPackageStartupMessages(library(dplyr))
if (file.exists('R/validation.R')) {
  source('R/validation.R')
  source('R/Micro.R')
} else {
  # R CMD check executes tests outside the source root: use its installed build.
  library(MicroBiotR)
  for (helper in c('.mbr_cluster_ids', '.mbr_cluster_match', '.mbr_event_ids')) {
    assign(helper, getFromNamespace(helper, 'MicroBiotR'))
  }
}
expect_error <- function(expr, pattern) {
  err <- tryCatch(expr, error = identity)
  stopifnot(inherits(err, 'error'), grepl(pattern, conditionMessage(err)))
}
stopifnot(identical(.mbr_cluster_ids(c('v1', 'V2', 'VitaminV', 'v01')),
                    c('V1', 'V2', 'VitaminV', 'V01')),
          identical(.mbr_cluster_match(c('V2', 'v1'), c('V1', 'v2')), c(2L, 1L)),
          identical(.mbr_event_ids(c(1, 1, 2)), c('V1', 'V1', 'V2')),
          identical(.mbr_event_ids(c('v1', 'v1', 'V2')), c('V1', 'V1', 'V2')))
expect_error(.mbr_cluster_ids(c('v1', 'V1')), 'ambiguous')
expect_error(.mbr_event_ids(c('v1', 'V1')), 'ambiguous')
expect_error(.mbr_event_ids(c('1', 'V1')), 'ambiguous')
expect_error(.mbr_cluster_match('V3', c('v1','v2')), 'missing')
expect_error(.mbr_cluster_match(c('v1','V1'), 'v1'), 'ambiguous')
expect_error(.mbr_cluster_match('vitaminV', 'VitaminV'), 'missing')
expect_error(.mbr_cluster_ids(c('v1', NA)), 'nonmissing')
# Preparation keeps input row order and names even with reversed mixed-case requests.
bins <- data.frame(x = 1:3, y = 4:6, row.names = c('v1', 'V2', 'v3'))
original <- bins
prepared <- MBR_prepare(bins, c('V3','V1'), transformation = identity)
stopifnot(identical(rownames(prepared), c('v1','v3')), identical(bins, original))
expect_error(MBR_prepare(bins, 'V9'), 'missing')
# Stop at the plotting boundary to inspect the actual selected input without
# requiring graphical dependencies; the production lookup path runs unchanged.
e <- new.env(parent = environment(MBR_violin))
e$MBR_violin <- MBR_violin; environment(e$MBR_violin) <- e
e$ggplot <- function(data, ...) {
  stopifnot(identical(data$abundance, 1:4))
  stop('plot boundary')
}
meta <- data.frame(Group=c('A','A','B','B'), row.names=paste0('s',1:4))
for (data_label in c('v1','V1')) for (p_label in c('v1','V1')) {
  dat <- setNames(data.frame(1:4, row.names=rownames(meta)), data_label)
  pv <- data.frame(p.value=.05, row.names=p_label)
  expect_error(e$MBR_violin(dat, meta, 'Group', pv), 'plot boundary')
}
expect_error(e$MBR_violin(dat, meta, 'Group', pv, cluster=9), 'missing')
pv_bad <- data.frame(p.value=c(.05,.1), row.names=c('v1','V1'))
expect_error(e$MBR_violin(dat, meta, 'Group', pv_bad), 'ambiguous')
# Check the p-value row index independently from the data column spelling/order.
e$format <- function(x, ...) { stopifnot(identical(x, .05)); stop('p-value verified') }
pv <- data.frame(p.value=c(.9,.05), row.names=c('v2','v1'))
expect_error(e$MBR_violin(dat, meta, 'Group', pv), 'p-value verified')
# Heatmap accepts either case on both sides and preserves the codebook row names.
h <- new.env(parent=environment(MBR_heatmap))
h$MBR_heatmap <- MBR_heatmap; environment(h$MBR_heatmap) <- h
h$pheatmap <- function(mat, ...) {
  stopifnot(identical(as.numeric(mat[,1]), c(20,10)))
  stop('heatmap boundary')
}
for (prefix in c('v','V')) {
  book <- data.frame(channel=c(10,20), row.names=paste0(prefix,1:2))
  abundance <- matrix(1,2,2,dimnames=list(NULL,c('v2','V1')))
  expect_error(h$MBR_heatmap(abundance, book), 'heatmap boundary')
}
expect_error(h$MBR_heatmap(abundance, book[1,,drop=FALSE]), 'missing')
# Circle labels change only full cluster IDs, never letters inside feature names.
circle <- new.env(parent=environment(MBR_circle))
circle$MBR_circle <- MBR_circle; environment(circle$MBR_circle) <- circle
circle$colorRamp2 <- function(...) NULL
circle$circos.clear <- circle$circos.par <- function(...) NULL
circle$circos.heatmap <- function(mat, ...) {
  stopifnot(identical(rownames(mat), c('Bin1','Bin2','VitaminV')))
  stop('circle boundary')
}
dat <- data.frame(v1=1:4,V2=2:5,VitaminV=3:6,row.names=rownames(meta))
expect_error(circle$MBR_circle(dat,meta,'Group'), 'circle boundary')
# Inspect actual ggplot layers to verify selected events and coordinate alignment.
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(viridis))
trans_format <- scales::trans_format
events <- data.frame(x=1:4, y=5:8, classes=c(1,1,2,3))
for (labels in list(c('v2','V1'), c('V2','v1'))) {
  plot <- suppressWarnings(MBR_flow_plot(events, bins=bins, x='x', y='y', selected_rows=labels))
  stopifnot(identical(plot$layers[[2]]$data$x, 1:3),
            identical(plot$layers[[3]]$data$x, c(2L,1L)))
}
expect_error(suppressWarnings(MBR_flow_plot(events, x='x', y='y', selected_rows='v9')), 'missing')
# Export real FCS files and reread them: each prefix case must retain the same
# event values, channel names and order, not merely the same number of rows.
suppressPackageStartupMessages(library(flowCore))
ff <- flowFrame(as.matrix(events[,c('x','y')]))
fs <- flowSet(list(example=ff))
folder <- tempfile('cluster-labels-'); dir.create(folder)
for (classes in list(events$classes, c('v1','v1','V2','v3'))) {
  events$classes <- classes
  for (selected in c('v1','V1')) {
    file <- MBR_save(fs, events, selected, rawdata_path=folder)
    exported <- exprs(read.FCS(file, transformation=FALSE))
    stopifnot(isTRUE(all.equal(unname(exported), unname(as.matrix(events[1:2,c('x','y')])), check.attributes=FALSE)),
              identical(unname(colnames(exported)), c('x','y')))
  }
}
expect_error(MBR_save(fs, events, 'V9', rawdata_path=folder), 'missing')
unlink(folder, recursive=TRUE)
cat('Cluster-label regression checks passed.\n')
