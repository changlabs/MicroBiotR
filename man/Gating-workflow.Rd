\name{Gating_workflow}
\alias{Gating_workflow}
\alias{gate_debris}
\alias{plot_gate_debr}
\alias{gate_dna}
\alias{plot_gate_dna}
\alias{plot_gating}
\alias{get_sc_gated}
\title{Standalone Instrument-Specific Gating Workflow and Helpers}
\description{
Documents the root-level \file{Gating} R script and its six helper functions.
These functions are not installed or exported by MicroBiotR; they exist only
when their script definitions are evaluated in an analysis session.
Sourcing the complete script immediately runs the analysis and writes files.
}
\usage{
gate_debris(data)
plot_gate_debr(data, gate, path = tmp_path)
gate_dna(data, channels = c("FSC PAR", "Hoechst Red"))
plot_gate_dna(data, gate, path = tmp_path)
plot_gating(x, nodes = gs_get_pop_paths(ig_gs))
get_sc_gated(data, node, channels = NULL)
}
\arguments{
\item{data}{A compatible flowWorkspace gating object for gating/plotting,
or one gating hierarchy for extraction. Unused by gate_dna.}
\item{gate}{Population name used as the ggcyto subset, for example root
or non-debris.}
\item{path}{Unused plotting-helper argument. Its default refers to the
script's tmp_path object, but is normally not evaluated.}
\item{channels}{For gate_dna, two coordinate channel names in order.
For get_sc_gated, optional channel names to retain; NULL retains all.}
\item{x}{A gating hierarchy accepted by ggcyto autoplot.}
\item{nodes}{Unused argument; does not select populations. Its default
references the script's ig_gs object.}
\item{node}{Population name from which events are extracted, such as DNA.}
}
\details{
The script attaches flowCore, flowWorkspace, openCyto, ggcyto, ggplot2
and ggthemes. It reads lowercase .fcs files from the working directory,
with transformation disabled and original channel names retained.
At least six samples are required because the preview selects six
without replacement after setting seed 2025.

Required raw channels include FSC PAR, SSC, Hoechst Red, FITC, APC,
Pe-TR and BV650. Marker labels are DNA, hIgA2, hIgA1, hIgG and hIgM
for the last five channels. Gates are fixed in raw acquisition units;
validate them against the actual instrument and experiment.

The non-debris population is a rectangular boundary on FSC PAR and SSC,
both 2000 to 63000, under root. The DNA population is a rectangle on
FSC PAR (2000 to 65000) and Hoechst Red (18000 to 67000), under
non-debris. The parent relationship intersects these gates. Recompute
applies the hierarchy to all samples. Gates define operational event
populations rather than identifying organisms or testing hypotheses.

The debris helper calls the private openCyto .boundary API, which is
version dependent. The DNA helper ignores its data argument. Plotting
helpers return hexagonal density figures with descriptive gate statistics;
on caught errors the debris and DNA plotting helpers return an empty plot
titled Error. Their path arguments are unused. plot_gating ignores nodes.
The surrounding workflow prints lists of plots; inspect actual PDFs to
verify that expected pages have been rendered.

get_sc_gated extracts a population's event matrix, optionally selects
channels, and appends marker labels when selected marker names match.
It returns a sample name and events suitable for the SOM input list.
Subsetting lacks drop = FALSE, so very small or one-channel subsets can
lose matrix shape.
}
\value{
gate_debris returns the boundary result supplied by openCyto;
gate_dna returns a flowCore polygonGate. The plotting helpers return
plot objects or their documented error fallback. get_sc_gated returns
list(name = sample_name, data = event_matrix).
The complete script leaves analysis objects in its evaluation environment
and saves the object fl_data_ig rather than returning a package result.
}
\section{Files and session effects}{
Creates /R_out/data and /R_out/graphs, requiring filesystem permission.
Writes graphs/gating/Ig_gate_debris.pdf, Ig_gate_dna.pdf and gating_Ig.pdf,
and data/FCM_gated/fl_data_ig.save under /R_out. Existing outputs may be
replaced. The saved object uses save/load, not saveRDS/readRDS.
The script attaches packages, sets the RNG seed and creates many objects.
It should not be sourced merely to load helper functions.
}
\examples{
\dontrun{
# First evaluate ONLY the relevant helper definitions from Gating
# in your analysis copy, with their dependencies attached.
g <- gate_dna(NULL)
stopifnot(inherits(g, "polygonGate"))

# ig_gs must already contain validated gates and real experimental data.
x <- get_sc_gated(ig_gs[[1]], node = "DNA",
                  channels = c("FSC PAR", "SSC"))
stopifnot(is.matrix(x$data), ncol(x$data) == 2)

# Inspect existing saved output without rerunning the workflow:
e <- new.env(parent = emptyenv())
load("/R_out/data/FCM_gated/fl_data_ig.save", envir = e)
stopifnot(is.list(e$fl_data_ig))
}
}
\seealso{\code{\link{MBR_som}}, \code{\link{MBR_read}}.
The repository guide \file{docs/Gating.md} documents each stage and helper.}
