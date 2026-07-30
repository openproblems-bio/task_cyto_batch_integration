# Compares compute_mem_matrix() / compute_mem_rmsd() (../helper.R) against
# the reference cytoMEM package (https://github.com/cytolab/cytoMEM), using
# real test data.
#
# Not part of the component's viash tests - this needs a GitHub install and
# is meant to be run by hand from the repo root whenever helper.R's MEM logic
# changes:
#   Rscript src/metrics/abundance_preservation/tests/test_mem_against_reference.R

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("cytoMEM", quietly = TRUE)) {
  remotes::install_github("cytolab/cytoMEM", upgrade = "never")
}
if (!requireNamespace("gplots", quietly = TRUE)) {
  install.packages("gplots", repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(cytoMEM)
  library(anndata)
})
source("src/metrics/abundance_preservation/helper.R")

TOLERANCE <- 1e-6
N_CELLS_PER_POP <- 100
POPULATIONS <- c("B cells", "CD4 Naive", "CD8 Naive")
MARKERS <- c("CD19#PE-Cy5", "CD3e#PE", "CD4#BUV496", "CD8a#BUV615", "CD62L#BV421")

set.seed(1)

adata <- anndata::read_h5ad(
  "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/unintegrated.h5ad"
)

is_candidate <- adata$obs$is_control == 0 & adata$obs$cell_type %in% POPULATIONS
subsampled_indices <- unlist(lapply(POPULATIONS, function(pop) {
  sample(which(is_candidate & adata$obs$cell_type == pop), N_CELLS_PER_POP)
}))

expr_matrix <- get_layer_expression(adata, "preprocessed", MARKERS, cell_mask = subsampled_indices)
population_labels <- as.character(adata$obs$cell_type[subsampled_indices])

my_mem <- compute_mem_matrix(expr_matrix, population_labels, POPULATIONS, iqr_thresh = 0.5)

# cytoMEM::MEM() requires the population label to be a numeric "cluster"
# column: apply() converts the whole data.frame to a matrix first, so a
# character cluster column would silently coerce every marker column to
# character too (and median()/IQR() on strings return NA). Population names
# are mapped to integer codes here purely to satisfy that, then mapped back
# for comparison.
cluster_ids <- setNames(seq_along(POPULATIONS), POPULATIONS)
ref_df <- as.data.frame(expr_matrix)
ref_df$cluster <- cluster_ids[population_labels]

mem_result <- MEM(
  ref_df, transform = FALSE, choose.markers = FALSE, markers = "all",
  choose.ref = FALSE, zero.ref = FALSE, rename.markers = FALSE,
  file.is.clust = FALSE, add.fileID = FALSE, IQR.thresh = 0.5,
  output.prescaled.MEM = FALSE, scale.matrix = "linear"
)
their_mem <- mem_result$MEM_matrix[[1]]
rownames(their_mem) <- names(cluster_ids)[match(as.numeric(rownames(their_mem)), cluster_ids)]

# --- Check 1: MEM score for a single population ---
for (target_pop in POPULATIONS) {
  mem_diff <- max(abs(my_mem[target_pop, ] - their_mem[target_pop, colnames(my_mem)]))
  cat(sprintf("MEM score max abs diff for '%s' vs cytoMEM::MEM(): %.2e\n", target_pop, mem_diff))
  stopifnot(mem_diff < TOLERANCE)
}


# --- Check 2: MEM RMSD between a pair of populations ---
pairs <- combn(POPULATIONS, 2, simplify = FALSE)
for (pair in pairs) {

  my_similarity <- compute_mem_rmsd(my_mem[pair[1], , drop = FALSE], my_mem[pair[2], , drop = FALSE])

  pdf_path <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_path)
  their_similarity <- MEM_RMSD(mem_result, output.matrix = FALSE)
  grDevices::dev.off()
  unlink(pdf_path)

  their_pair_ids <- as.character(cluster_ids[pair])
  their_similarity_pair <- their_similarity[their_pair_ids[1], their_pair_ids[2]]

  rmsd_diff <- abs(my_similarity[1, 1] - their_similarity_pair)
  cat(sprintf("RMSD similarity max abs diff for '%s' vs '%s' vs cytoMEM::MEM_RMSD(): %.2e\n", pair[1], pair[2], rmsd_diff))
  stopifnot(rmsd_diff < TOLERANCE)
}

cat("\nAll MEM / MEM RMSD checks passed against the reference cytoMEM package.\n")
