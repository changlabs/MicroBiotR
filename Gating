# ----------------------------
# Load libraries
# ----------------------------
library(flowCore)
library(flowWorkspace)
library(openCyto)
library(ggcyto)
library(ggplot2)
library(ggthemes)

# ----------------------------
# Define paths
# ----------------------------
data_path <- "/R_out/data"
graph_path <- "/R_out/graphs"
rawdata_path <- "./"

dir.create(data_path, recursive = TRUE, showWarnings = FALSE)
dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Load FCS data
# ----------------------------
fcs_files <- list.files(rawdata_path, pattern = "\\.fcs$", full.names = TRUE)
ig_fs <- read.flowSet(files = fcs_files, column.pattern = "\\*|Bits|Drop",
                      invert.pattern = TRUE, transformation = FALSE,
                      alter.names = FALSE)

tmp_st <- c("DNA", "hIgA2", "hIgA1", "hIgG", "hIgM")
names(tmp_st) <- c("Hoechst Red", "FITC", "APC", "Pe-TR", "BV650")
markernames(ig_fs) <- tmp_st

ig_gs <- GatingSet(ig_fs)

set.seed(2025)
tmp_sam <- sort(sample(seq_len(length(ig_gs)), 6, replace = FALSE))

# ----------------------------
# Debris gate
# ----------------------------
p <- ggcyto(data = gs_pop_get_data(ig_gs[tmp_sam]),
            mapping = aes(x = "FSC PAR", y = "SSC")) +
  geom_hex(bins = 128)

gate_debris <- function(data) {
  openCyto:::.boundary(data, channels = c("FSC PAR", "SSC"),
                       min = c(2000, 2000), max = c(6.3e4, 6.3e4))
}

p + geom_gate(gate_debris(ig_gs[tmp_sam]))

gs_pop_add(ig_gs, gate_debris(ig_gs), name = "non-debris", parent = "root")
recompute(ig_gs)

# ----------------------------
# Plot debris gate
# ----------------------------
tmp_path <- file.path(graph_path, "gating/")
dir.create(tmp_path, recursive = TRUE, showWarnings = FALSE)

pdf(file = file.path(tmp_path, "Ig_gate_debris.pdf"), onefile = TRUE,
    bg = "transparent", width = 7, height = 7)

plot_gate_debr <- function(data, gate, path = tmp_path) {
  tryCatch({
    ggcyto(data, subset = gate,
           mapping = aes(x =  "FSC PAR", y = "SSC")) +
      geom_hex(bins = 128) +
      geom_gate() +
      geom_stats()
  }, error = function(e) {
    message("⚠️ Error in plot_gate_debr: ", conditionMessage(e))
    ggplot() + ggtitle("Error")
  })
}

print(lapply(ig_gs, plot_gate_debr, gate = "root"))
dev.off()

# ----------------------------
# DNA gate
# ----------------------------
p <- ggcyto(data = ig_gs[tmp_sam],
            subset = "non-debris",
            mapping = aes(x = "FSC PAR", y = "DNA")) +
  geom_hex(bins = 128)

gate_dna <- function(data, channels = c("FSC PAR", "Hoechst Red")) {
  tmp_mat <- matrix(c(
    2.0e3, 6.5e4, 6.5e4, 2.0e3,   # FSC PAR
    1.8e4, 1.8e4, 6.7e4, 6.7e4    # DNA
  ), ncol = 2)
  colnames(tmp_mat) <- channels
  polygonGate(tmp_mat)
}

p + geom_gate(gate_dna(ig_gs))

gs_pop_add(ig_gs, gate_dna(ig_gs), name = "DNA", parent = "non-debris")
recompute(ig_gs)

pdf(file = file.path(tmp_path, "Ig_gate_dna.pdf"), onefile = TRUE,
    bg = "transparent", width = 7, height = 7)

plot_gate_dna <- function(data, gate, path = tmp_path) {
  tryCatch({
    ggcyto(data, subset = gate,
           mapping = aes(x =  "FSC PAR", y = "Hoechst Red")) +
      geom_hex(bins = 128) +
      geom_gate() +
      geom_stats()
  }, error = function(e) {
    message("⚠️ Error in plot_gate_dna: ", conditionMessage(e))
    ggplot() + ggtitle("Error")
  })
}

print(lapply(ig_gs, plot_gate_dna, gate = "non-debris"))
dev.off()

# ----------------------------
# Autoplot: show all gating
# ----------------------------
plot_gating <- function(x, nodes = gs_get_pop_paths(ig_gs)) {
  autoplot(x, bins = 64, nrow = 1, ylim = "instrument", xlim = "instrument")
}

pdf(file = file.path(tmp_path, "gating_Ig.pdf"), width = 3,
    height = 9, bg = "transparent", onefile = TRUE)
print(lapply(ig_gs, plot_gating))
dev.off()

# ----------------------------
# Extract gated single-cell data
# ----------------------------
get_sc_gated <- function(data, node, channels = NULL) {
  tmp_name <- pData(data[[1]])$name
  tmp_mat <- exprs(gh_pop_get_data(data, y = node))
  if (!is.null(channels)) {
    tmp_mat <- tmp_mat[, channels]
  }
  tmp_channels <- markernames(data)[names(markernames(data)) %in% channels]
  colnames(tmp_mat) <- sapply(colnames(tmp_mat), function(x) {
    ifelse(!is.na(tmp_channels[x]), paste0(x, ".", tmp_channels[x]), x)
  })
  list(name = tmp_name, data = tmp_mat)
}

fl_data_ig <- lapply(ig_gs, get_sc_gated,
                     node = "DNA",
                     channels = c("FSC PAR", "SSC", "Hoechst Red", "FITC", "APC", "Pe-TR", "BV650"))

# ----------------------------
# Save output
# ----------------------------
save_to <- file.path(data_path, "FCM_gated")
dir.create(save_to, recursive = TRUE, showWarnings = FALSE)
save(fl_data_ig, file = file.path(save_to, "fl_data_ig.save"), compress = TRUE)
cat("\n✅ DONE: Gated flow cytometry data saved.\n")
