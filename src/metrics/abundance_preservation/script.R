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
# is present in the freshly-read integrated splits (controls included).
unintegrated_raw <- anndata::read_h5ad(par[["input_unintegrated"]])

integrated_s1 <- anndata::read_h5ad(par[["input_integrated_split1"]]) |>
  get_obs_var_for_integrated(u_adata = unintegrated_raw, split_id = 1) |>
  subset_markers_tocorrect() |>
  subset_nocontrols()

integrated_s2 <- anndata::read_h5ad(par[["input_integrated_split2"]]) |>
  get_obs_var_for_integrated(u_adata = unintegrated_raw, split_id = 2) |>
  subset_markers_tocorrect() |>
  subset_nocontrols()

unintegrated <- unintegrated_raw |> subset_nocontrols()

print("Setting parameters\n")

# flowsom parameters
fsom_param <- get_flowsom_parameters(unintegrated_raw)

# dataset specific parameters
dataset_param <- get_dataset_parameters(unintegrated)

print("Cluster each unintegrated batch separately and compute MEM scores\n")
# Each batch is clustered on its own, so it is split off into its own AnnData
# first - a cluster's MEM score is relative to the other cells it was
# clustered alongside, which here should only ever be cells from the same batch.
unintegrated_b1 <- unintegrated[unintegrated$obs$batch == dataset_param$batch_1, ]$copy()
unintegrated_b2 <- unintegrated[unintegrated$obs$batch == dataset_param$batch_2, ]$copy()

clustering_b1 <- cluster_and_calculate_mem(
  adata = unintegrated_b1,
  layer_name = "preprocessed",
  lineage_markers = dataset_param$lineage_markers,
  fsom_param = fsom_param,
  label_prefix = paste0("batch", dataset_param$batch_1, "_c")
)
clustering_b2 <- cluster_and_calculate_mem(
  adata = unintegrated_b2,
  layer_name = "preprocessed",
  lineage_markers = dataset_param$lineage_markers,
  fsom_param = fsom_param,
  label_prefix = paste0("batch", dataset_param$batch_2, "_c")
)

print("Test which clusters are differentially abundant within each batch\n")
abundance_b1 <- test_differential_abundance(
  cluster_id_per_cell = clustering_b1$cluster_id_per_cell,
  sample_id_per_cell = unintegrated_b1$obs$sample,
  group_id_per_cell = unintegrated_b1$obs$group,
  group_a = dataset_param$group_a,
  group_b = dataset_param$group_b,
  significance_threshold = SIGNIFICANCE_THRESHOLD
)

abundance_b2 <- test_differential_abundance(
  cluster_id_per_cell = clustering_b2$cluster_id_per_cell,
  sample_id_per_cell = unintegrated_b2$obs$sample,
  group_id_per_cell = unintegrated_b2$obs$group,
  group_a = dataset_param$group_a,
  group_b = dataset_param$group_b,
  significance_threshold = SIGNIFICANCE_THRESHOLD
)

# keep just clusters which are DA
DA_clusters_b1 <- abundance_b1$significant
DA_clusters_b2 <- abundance_b2$significant


print("Finding best matches for batch 1 clusters\n")
batch_matches <- find_cluster_match(
  mem_a = clustering_b1$mem,
  mem_b = clustering_b2$mem
)
cluster_b1_best_matches <- batch_matches$cluster_a_matches
# keep only clusters which are we found "mutual" matches,
# i.e., if b1 cluster 1 best match in b2 is cluster 2,
# and b2 cluster 2 best match in b1 is cluster 1,
# then this is mutual matches. Otherwise, there were some discrepancies in 
# the matching that is hard to resolve (which one is the right one?)
# in that case, we just drop them.
cluster_b1_best_matches_filtered <- cluster_b1_best_matches[cluster_b1_best_matches$is_mutual, , drop = FALSE]

# now we check which of the mutually matched clusters are differentially 
# abundant in both batches
cluster_b1_best_matches_filtered <- cluster_b1_best_matches_filtered[
  cluster_b1_best_matches_filtered$cluster_a %in% DA_clusters_b1 &
    cluster_b1_best_matches_filtered$best_match_b %in% DA_clusters_b2,
  ,
  drop = FALSE
]

# Nothing survived both filters, so there is nothing to score against. This
# should not happen on a real dataset, so stop rather than report NaN.
if (nrow(cluster_b1_best_matches_filtered) == 0) {
  stop(
    "No batch 1 cluster was both mutually matched to a batch 2 cluster and ",
    "differentially abundant in both batches, so there is nothing to score."
  )
}

print("Clustering and running DA for each integrated split\n")
result_s1 <- process_integrated_split(
  adata = integrated_s1,
  dataset_param = dataset_param,
  fsom_param = fsom_param,
  b1_mem = clustering_b1$mem,
  retained_cluster_b1 = cluster_b1_best_matches_filtered$cluster_a,
  significance_threshold = SIGNIFICANCE_THRESHOLD,
  label_prefix = "split1_c"
)
result_s2 <- process_integrated_split(
  adata = integrated_s2,
  dataset_param = dataset_param,
  fsom_param = fsom_param,
  b1_mem = clustering_b1$mem,
  retained_cluster_b1 = cluster_b1_best_matches_filtered$cluster_a,
  significance_threshold = SIGNIFICANCE_THRESHOLD,
  label_prefix = "split2_c"
)

print("Computing final score\n")
# The retained b1 clusters are the ones carried into the splits, so they are
# what each split is scored against: a method that kept every one of them
# differentially abundant scores 1.
retained_clusters_b1 <- cluster_b1_best_matches_filtered$cluster_a

prop_still_da_s1 <- length(result_s1$da_results$significant) / length(retained_clusters_b1)
prop_still_da_s2 <- length(result_s2$da_results$significant) / length(retained_clusters_b1)
score <- mean(c(prop_still_da_s1, prop_still_da_s2))

print("Write output AnnData to file\n")
output <- anndata::AnnData(
  shape = c(0L, 0L),
  uns = list(
    dataset_id = integrated_s1$uns$dataset_id,
    method_id = integrated_s1$uns$method_id,
    metric_ids = list("abundance_preservation"),
    metric_values = list(score),
    group_a = dataset_param$group_a,
    group_b = dataset_param$group_b,
    da_clusters_batch1 = as.list(DA_clusters_b1),
    da_clusters_batch2 = as.list(DA_clusters_b2),
    pvals_batch1 = as.list(abundance_b1$pvals),
    pvals_batch2 = as.list(abundance_b2$pvals),
    retained_clusters_batch1 = as.list(retained_clusters_b1),
    retained_clusters_batch2 = as.list(cluster_b1_best_matches_filtered$best_match_b),
    batch_cluster_similarity = batch_matches$similarity,
    cluster_similarity_split1 = result_s1$similarity_to_batch1,
    cluster_similarity_split2 = result_s2$similarity_to_batch1,
    pvals_split1 = as.list(result_s1$da_results$pvals),
    pvals_split2 = as.list(result_s2$da_results$pvals),
    da_clusters_split1 = as.list(result_s1$da_results$significant),
    da_clusters_split2 = as.list(result_s2$da_results$significant),
    proportion_split1 = prop_still_da_s1,
    proportion_split2 = prop_still_da_s2,
    fsom_parameters = fsom_param
  )
)

output$write_h5ad(par[["output"]], compression = "gzip")
