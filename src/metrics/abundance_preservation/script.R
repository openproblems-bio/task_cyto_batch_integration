library(anndata, quietly = TRUE)
library(FlowSOM, quietly = TRUE)
library(flowCore, quietly = TRUE)

## VIASH START
par <- list(
  "input_unintegrated" = 'resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/unintegrated.h5ad',
  "input_integrated_split1" = 'resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/integrated_split1.h5ad',
  "input_integrated_split2" = 'resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/integrated_split2.h5ad',
  "output" = 'resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/score.h5ad'
)
meta <- list(
  "name" = 'abundance_preservation',
  "resources_dir" = "src/utils"
)
## VIASH END

source(paste0(meta$resources_dir, "/helper.R"))
source(paste0(meta$resources_dir, "/helper_functions.R"))

# With only a handful of samples per group per batch, the smallest two-sided
# Wilcoxon p-value that can ever be reached is around 0.1, so that is the
# threshold used to call a difference "significant" (same rationale as the
# functional_marker_preservation metric).
SIGNIFICANCE_THRESHOLD <- 0.1

print("Reading input files\n")
# Keep a raw, untouched copy of the unintegrated data: get_obs_var_for_integrated
# matches cells by name against it, so it must still contain every cell that
# is present in the freshly-read integrated splits (controls and unlabelled
# cells included).
unintegrated_raw <- anndata::read_h5ad(par[["input_unintegrated"]])

integrated_s1 <- anndata::read_h5ad(par[["input_integrated_split1"]]) |>
  get_obs_var_for_integrated(unintegrated_raw, split_id = 1) |>
  subset_markers_tocorrect() |>
  subset_nocontrols()

integrated_s2 <- anndata::read_h5ad(par[["input_integrated_split2"]]) |>
  get_obs_var_for_integrated(unintegrated_raw, split_id = 2) |>
  subset_markers_tocorrect() |>
  subset_nocontrols()

unintegrated <- unintegrated_raw |> subset_nocontrols()

print("Setting parameters\n")
lineage_markers <- rownames(unintegrated_raw$var)[unintegrated_raw$var$marker_type == "lineage"]
n_clusters <- unintegrated_raw$uns$parameter_num_clusters
grid_xdim <- unintegrated_raw$uns$parameter_som_xdim
grid_ydim <- unintegrated_raw$uns$parameter_som_ydim
real_cell_types <- unique(unintegrated$obs$cell_type)
real_cell_types <- real_cell_types[!tolower(real_cell_types) %in% c("unlabelled", "unlabeled")]
groups <- sort(as.character(unique(unintegrated$obs$group)))
group_a <- groups[1]
group_b <- groups[2]
batches <- sort(unique(unintegrated$obs$batch))

print("Testing baseline (unintegrated) abundance differences, per batch\n")
# Unlabelled cells count toward each sample's total (so real cell types'
# proportions reflect the true unannotated fraction) but are never tested
# themselves - run_wilcoxon_per_celltype skips them.
sig_per_batch <- lapply(batches, function(current_batch) {
  in_batch <- unintegrated$obs$batch == current_batch
  proportions <- compute_sample_proportions(
    cell_type = unintegrated$obs$cell_type[in_batch],
    sample = unintegrated$obs$sample[in_batch],
    group = unintegrated$obs$group[in_batch]
  )
  pvals <- run_wilcoxon_per_celltype(proportions, group_a, group_b)
  names(pvals)[pvals <= SIGNIFICANCE_THRESHOLD]
})
names(sig_per_batch) <- as.character(batches)
# A cell type only counts as a real baseline difference if it shows up in both batches.
sig_unintegrated <- Reduce(intersect, sig_per_batch)

print("Computing reference MEM scores for each cell type, averaged over batches\n")
mem_per_batch <- lapply(batches, function(current_batch) {
  in_batch <- unintegrated$obs$batch == current_batch
  expr <- get_layer_expression(unintegrated, "preprocessed", lineage_markers, cell_mask = in_batch)
  compute_mem_matrix(expr, unintegrated$obs$cell_type[in_batch], real_cell_types)
})
mem_celltype_reference <- Reduce(`+`, mem_per_batch) / length(mem_per_batch)

print("Clustering and testing each integrated split\n")
result_s1 <- process_integrated_split(
  integrated_s1, lineage_markers, n_clusters, grid_xdim, grid_ydim,
  mem_celltype_reference, group_a, group_b, SIGNIFICANCE_THRESHOLD
)
result_s2 <- process_integrated_split(
  integrated_s2, lineage_markers, n_clusters, grid_xdim, grid_ydim,
  mem_celltype_reference, group_a, group_b, SIGNIFICANCE_THRESHOLD
)

# A cell type only counts as "still significant post-integration" if it
# remains significant in both splits.
sig_integrated <- intersect(result_s1$sig, result_s2$sig)
sig_remaining <- intersect(sig_unintegrated, sig_integrated)

print("Computing final score\n")
score <- if (length(sig_unintegrated) == 0) {
  NA_real_
} else {
  length(sig_remaining) / length(sig_unintegrated)
}

print("Write output AnnData to file\n")
output <- anndata::AnnData(
  shape = c(0L, 0L),
  uns = list(
    dataset_id = integrated_s1$uns$dataset_id,
    method_id = integrated_s1$uns$method_id,
    metric_ids = list("abundance_preservation"),
    metric_values = list(score),
    group_a = group_a,
    group_b = group_b,
    sig_unintegrated = as.list(sig_unintegrated),
    sig_split1 = as.list(result_s1$sig),
    sig_split2 = as.list(result_s2$sig),
    sig_remaining = as.list(sig_remaining),
    metacluster_mapping_split1 = as.list(result_s1$mapping),
    metacluster_mapping_split2 = as.list(result_s2$mapping),
    metacluster_similarity_split1 = result_s1$similarity,
    metacluster_similarity_split2 = result_s2$similarity,
    fsom_parameters = list(
      "xdim" = grid_xdim,
      "ydim" = grid_ydim,
      "n_clusters" = n_clusters
    )
  )
)

output$write_h5ad(par[["output"]], compression = "gzip")
