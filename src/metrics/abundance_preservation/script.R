library(anndata, quietly = TRUE)
library(FlowSOM, quietly = TRUE)
library(flowCore, quietly = TRUE)

## VIASH START
home_path <- file.path("/Users/putri.g/Documents/GitHub/task_cyto_batch_integration/resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset")
par <- list(
  "input_unintegrated" = file.path(home_path, "unintegrated.h5ad"),
  "input_integrated_split1" = file.path(home_path, "integrated_split1.h5ad"),
  "input_integrated_split2" = file.path(home_path, "integrated_split2.h5ad"),
  "output" = file.path(home_path, "score.h5ad")
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
b1_cluster_matches <- find_cluster_match(
  mem_a = clustering_b1$mem,
  mem_b = clustering_b2$mem
)
cluster_b1_best_matches <- b1_cluster_matches$cluster_a_matches
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
# Clusters in integrated are mapped to retained clusters in b1.
# So it only makes sense for the final proportion to be calculated over the retained clusters in b1, not all clusters in b1.
retained_clusters_b1 <- cluster_b1_best_matches_filtered$cluster_a

prop_still_da_s1 <- length(result_s1$da_results$significant) / length(retained_clusters_b1)
prop_still_da_s2 <- length(result_s2$da_results$significant) / length(retained_clusters_b1)
score <- mean(c(prop_still_da_s1, prop_still_da_s2))

print("Write output AnnData to file\n")
# Everything below is grouped by the stage it came from, so a diagnosis never
# needs to stitch loose vectors back together. Similarity matrices are written
# as data frames because a plain matrix loses its dimnames in h5ad, which makes
# the cluster labels unrecoverable.

# One row per unintegrated cluster, saying whether it was differentially abundant.
da_status_for_clusters_in_unintegrated <- rbind(
  data.frame(
    cluster_id = abundance_b1$clusters_tested,
    batch = as.character(dataset_param$batch_1),
    is_da = abundance_b1$clusters_tested %in% DA_clusters_b1,
    stringsAsFactors = FALSE
  ),
  data.frame(
    cluster_id = abundance_b2$clusters_tested,
    batch = as.character(dataset_param$batch_2),
    is_da = abundance_b2$clusters_tested %in% DA_clusters_b2,
    stringsAsFactors = FALSE
  )
)

# What each batch 1 cluster matched to in batch 2, and whether that was mutual.
unintegrated_cluster_match_table <- data.frame(
  cluster_b1 = cluster_b1_best_matches$cluster_a,
  best_match_b2 = cluster_b1_best_matches$best_match_b,
  b2_best_match_back_in_b1 = cluster_b1_best_matches$best_match_b_matching_a,
  is_mutual = cluster_b1_best_matches$is_mutual,
  row.names = NULL,
  stringsAsFactors = FALSE
)

# One row per retained batch 1 cluster per split. A retained cluster that no
# split cluster mapped onto is never tested, and counts as not differentially
# abundant.
cluster_da_split1 <- data.frame(
  cluster_b1 = retained_clusters_b1,
  was_tested = retained_clusters_b1 %in% result_s1$da_results$clusters_tested,
  is_da = retained_clusters_b1 %in% result_s1$da_results$significant,
  stringsAsFactors = FALSE
)
cluster_da_split2 <- data.frame(
  cluster_b1 = retained_clusters_b1,
  was_tested = retained_clusters_b1 %in% result_s2$da_results$clusters_tested,
  is_da = retained_clusters_b1 %in% result_s2$da_results$significant,
  stringsAsFactors = FALSE
)

output <- anndata::AnnData(
  shape = c(0L, 0L),
  uns = list(
    dataset_id = integrated_s1$uns$dataset_id,
    method_id = integrated_s1$uns$method_id,
    metric_ids = list("abundance_preservation"),
    metric_values = list(score),

    parameters = list(
      group_a = dataset_param$group_a,
      group_b = dataset_param$group_b,
      batch_1 = as.character(dataset_param$batch_1),
      batch_2 = as.character(dataset_param$batch_2),
      significance_threshold = SIGNIFICANCE_THRESHOLD,
      som_xdim = fsom_param$xdim,
      som_ydim = fsom_param$ydim,
      n_clusters = fsom_param$n_clusters
    ),

    unintegrated = list(
      # cell_id and cluster_id indices are matching: cluster_id at index i belongs to cell_id at index i.
      batch1_cell_id = rownames(unintegrated_b1$obs),
      batch1_cluster_id = clustering_b1$cluster_id_per_cell,
      batch2_cell_id = rownames(unintegrated_b2$obs),
      batch2_cluster_id = clustering_b2$cluster_id_per_cell,
      # which clusters in b1 and b2 are DA? or not DA?
      which_clusters_are_da = da_status_for_clusters_in_unintegrated,
      # rows are batch 1 clusters, columns are batch 2 clusters
      b1_b2_cluster_similarity = as.data.frame(b1_cluster_matches$similarity),
      cluster_match_table = unintegrated_cluster_match_table,
      retained_cluster_b1 = retained_clusters_b1
    ),

    split1 = list(
      cell_id = rownames(integrated_s1$obs),
      cluster_id = result_s1$cluster_id_per_cell,
      mapped_b1_cluster_id = result_s1$mapped_b1_cluster,
      # rows are this split's clusters, columns are batch 1 clusters
      similarity_to_b1 = as.data.frame(result_s1$similarity_to_batch1),
      # note the cluster id is the mapped b1 clusters
      which_clusters_are_da = cluster_da_split1,
      proportion_still_da = prop_still_da_s1
    ),

    split2 = list(
      cell_id = rownames(integrated_s2$obs),
      cluster_id = result_s2$cluster_id_per_cell,
      mapped_b1_cluster_id = result_s2$mapped_b1_cluster,
      similarity_to_b1 = as.data.frame(result_s2$similarity_to_batch1),
      # note the cluster id is the mapped b1 clusters
      which_clusters_are_da = cluster_da_split2,
      proportion_still_da = prop_still_da_s2
    )
  )
)

output$write_h5ad(par[["output"]], compression = "gzip")
