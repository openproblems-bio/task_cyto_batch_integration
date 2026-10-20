library(anndata, quietly = TRUE)
library(FlowSOM, quietly = TRUE)
library(flowCore, quietly = TRUE)

## VIASH START
home_path <- file.path("resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset")
par <- list(
  "input_unintegrated" = file.path(home_path, "unintegrated.h5ad"),
  "output" = file.path(home_path, "precomputed_r_metric.h5ad")
)
meta <- list(
  "name" = "precompute_r_metric",
  "resources_dir" = "src/utils"
)
## VIASH END

source(paste0(meta$resources_dir, "/helper_functions.R"))
# helper.R from the abundance_preservation metric.
# Viash copies it into resources_dir.
source(paste0(meta$resources_dir, "/helper.R"))

cat("Reading input file\n")
unintegrated <- anndata::read_h5ad(par[["input_unintegrated"]]) |>
  subset_nocontrols()

# ── Abundance preservation ──────────────────────────────────────────────────
SIGNIFICANCE_THRESHOLD <- 0.1

cat("Setting parameters\n")
fsom_param <- get_flowsom_parameters(unintegrated)
dataset_param <- get_dataset_parameters(unintegrated)

cat("Cluster each unintegrated batch separately and compute MEM scores\n")
unintegrated_b1 <- unintegrated[unintegrated$obs$batch == dataset_param$batch_1, ]$copy()
unintegrated_b2 <- unintegrated[unintegrated$obs$batch == dataset_param$batch_2, ]$copy()

cat("Clustering and computing MEM scores for batch 1\n")
clustering_b1 <- cluster_and_calculate_mem(
  adata = unintegrated_b1,
  layer_name = "preprocessed",
  lineage_markers = dataset_param$lineage_markers,
  fsom_param = fsom_param,
  label_prefix = paste0("batch", dataset_param$batch_1, "_c")
)
cat("Clustering and computing MEM scores for batch 2\n")
clustering_b2 <- cluster_and_calculate_mem(
  adata = unintegrated_b2,
  layer_name = "preprocessed",
  lineage_markers = dataset_param$lineage_markers,
  fsom_param = fsom_param,
  label_prefix = paste0("batch", dataset_param$batch_2, "_c")
)

cat("Test which clusters are differentially abundant within each batch\n")
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

cat("Finding best matches for batch 1 clusters\n")
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
retained_clusters_b1 <- cluster_b1_best_matches_filtered$cluster_a

# If nothing is retained, the metric stops with an error. That check stays in
# the metric so a failure here doesn't take down the other metrics that run
# on this dataset.
print(sprintf("Retained %d batch 1 clusters\n", length(retained_clusters_b1)))

cat("Write output AnnData to file\n")
# Matrices are written as data frames because a plain matrix loses its
# dimnames in h5ad, which makes the cluster labels unrecoverable.

# One row per unintegrated cluster, saying whether it was differentially abundant.
da_status_for_clusters <- rbind(
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

# One row per unintegrated cell, with the batch and the cluster it was assigned to.
cluster_id_per_cell <- rbind(
  data.frame(
    cell_id = rownames(unintegrated_b1$obs),
    batch = as.character(dataset_param$batch_1),
    cluster_id = clustering_b1$cluster_id_per_cell,
    stringsAsFactors = FALSE
  ),
  data.frame(
    cell_id = rownames(unintegrated_b2$obs),
    batch = as.character(dataset_param$batch_2),
    cluster_id = clustering_b2$cluster_id_per_cell,
    stringsAsFactors = FALSE
  )
)

output <- anndata::AnnData(
  shape = c(0L, 0L),
  uns = list(
    # Results for each R metric are kept under that metric's name, so more
    # R metrics can add their own entries to this file.
    abundance_preservation = list(
      parameters = list(
        group_a = dataset_param$group_a,
        group_b = dataset_param$group_b,
        batch_1 = as.character(dataset_param$batch_1),
        batch_2 = as.character(dataset_param$batch_2),
        lineage_markers = dataset_param$lineage_markers,
        significance_threshold = SIGNIFICANCE_THRESHOLD,
        som_xdim = fsom_param$xdim,
        som_ydim = fsom_param$ydim,
        n_clusters = fsom_param$n_clusters
      ),

      # rows are batch 1 clusters, columns are markers.
      # The metric maps each integrated split cluster onto these.
      batch1_mem = as.data.frame(clustering_b1$mem),
      # batch 1 clusters that are mutually matched to a batch 2 cluster and
      # differentially abundant in both batches.
      retained_cluster_b1 = retained_clusters_b1,

      # Not used by the metric. Kept to check why a cluster was or wasn't retained.
      diagnostics = list(
        cluster_id_per_cell = cluster_id_per_cell,
        # which clusters in b1 and b2 are DA? or not DA?
        which_clusters_are_da = da_status_for_clusters,
        cluster_match_table = unintegrated_cluster_match_table,
        # rows are batch 1 clusters, columns are batch 2 clusters
        b1_b2_cluster_similarity = as.data.frame(b1_cluster_matches$similarity)
      )
    )
  )
)

output$write_h5ad(par[["output"]], compression = "gzip")
