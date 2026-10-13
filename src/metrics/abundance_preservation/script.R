library(anndata, quietly = TRUE)
library(FlowSOM, quietly = TRUE)
library(flowCore, quietly = TRUE)

## VIASH START
home_path <- file.path("/Users/putri.g/Documents/GitHub/task_cyto_batch_integration/resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset")
par <- list(
  "input_unintegrated" = file.path(home_path, "unintegrated.h5ad"),
  "input_integrated_split1" = file.path(home_path, "integrated_split1.h5ad"),
  "input_integrated_split2" = file.path(home_path, "integrated_split2.h5ad"),
  "input_precomputed_r_metric" = file.path(home_path, "precomputed_r_metric.h5ad"),
  "output" = file.path(home_path, "score.h5ad")
)
meta <- list(
  "name" = 'abundance_preservation',
  "resources_dir" = "src/utils"
)
## VIASH END

source(paste0(meta$resources_dir, "/helper_functions.R"))
source(paste0(meta$resources_dir, "/helper.R"))

print("Reading input files\n")

unintegrated <- anndata::read_h5ad(par[["input_unintegrated"]])

integrated_s1 <- anndata::read_h5ad(par[["input_integrated_split1"]]) |>
  get_obs_var_for_integrated(u_adata = unintegrated, split_id = 1) |>
  subset_markers_tocorrect() |>
  subset_nocontrols()

integrated_s2 <- anndata::read_h5ad(par[["input_integrated_split2"]]) |>
  get_obs_var_for_integrated(u_adata = unintegrated, split_id = 2) |>
  subset_markers_tocorrect() |>
  subset_nocontrols()

# The unintegrated half of the metric (per-batch clustering, DA tests and
# cluster matching) was pre-computed once per
# dataset by the precompute_r_metric component.
print("Reading precomputed unintegrated results\n")
precomputed <- anndata::read_h5ad(par[["input_precomputed_r_metric"]])$uns$abundance_preservation
precomputed_param <- precomputed$parameters

SIGNIFICANCE_THRESHOLD <- precomputed_param$significance_threshold

dataset_param <- list(
  lineage_markers = precomputed_param$lineage_markers,
  group_a = precomputed_param$group_a,
  group_b = precomputed_param$group_b,
  batch_1 = precomputed_param$batch_1,
  batch_2 = precomputed_param$batch_2
)
fsom_param <- list(
  n_clusters = precomputed_param$n_clusters,
  xdim = precomputed_param$som_xdim,
  ydim = precomputed_param$som_ydim
)

# The MEM matrix was stored as a data frame to keep its cluster and marker labels.
b1_mem <- as.matrix(precomputed$batch1_mem)
retained_clusters_b1 <- precomputed$retained_cluster_b1

# Nothing survived both filters, so there is nothing to score against. This
# should not happen on a real dataset, so stop rather than report NaN.
if (length(retained_clusters_b1) == 0) {
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
  b1_mem = b1_mem,
  retained_cluster_b1 = retained_clusters_b1,
  significance_threshold = SIGNIFICANCE_THRESHOLD,
  label_prefix = "split1_c"
)
result_s2 <- process_integrated_split(
  adata = integrated_s2,
  dataset_param = dataset_param,
  fsom_param = fsom_param,
  b1_mem = b1_mem,
  retained_cluster_b1 = retained_clusters_b1,
  significance_threshold = SIGNIFICANCE_THRESHOLD,
  label_prefix = "split2_c"
)

print("Computing final score\n")
# Clusters in integrated are mapped to retained clusters in b1.
# So it only makes sense for the final proportion to be calculated over the retained clusters in b1, not all clusters in b1.
prop_still_da_s1 <- length(result_s1$da_results$significant) / length(retained_clusters_b1)
prop_still_da_s2 <- length(result_s2$da_results$significant) / length(retained_clusters_b1)
score <- mean(c(prop_still_da_s1, prop_still_da_s2))

print("Write output AnnData to file\n")
# Everything below is grouped by the stage it came from, so a diagnosis never
# needs to stitch loose vectors back together. Similarity matrices are written
# as data frames because a plain matrix loses its dimnames in h5ad, which makes
# the cluster labels unrecoverable.
# The unintegrated results are in the precomputed file, not repeated here.

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
