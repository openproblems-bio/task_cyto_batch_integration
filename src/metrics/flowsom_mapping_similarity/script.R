library(anndata, quietly = TRUE)
library(FlowSOM, quietly = TRUE)

## VIASH START
par <- list(
  "input_unintegrated" = 'resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/unintegrated.h5ad',
  "input_integrated_split1" = 'resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/integrated_split1.h5ad',
  "input_integrated_split2" = 'resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/integrated_split2.h5ad',
  # if using perfect integration
  # input_integrated_split1 = "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/perfect_integrated_split1.h5ad",
  # input_integrated_split2 = "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/perfect_integrated_split2.h5ad",
  "output" = 'resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/score.h5ad'
)
meta <- list(
  "name" = 'flowsom_mapping_similarity',
  "resources_dir" = "src/utils"
)
## VIASH END

source(paste0(meta$resources_dir, "/helper.R"))
source(paste0(meta$resources_dir, "/helper_functions.R"))
library(anndata)

unintegrated <- anndata::read_h5ad(par[["input_unintegrated"]])

# read and filter split 1 data
integrated_s1 <- anndata::read_h5ad(par[["input_integrated_split1"]]) |>
  get_obs_var_for_integrated(unintegrated, split_id = 1) |>
  subset_markers_tocorrect() |>
  subset_nocontrols() |>
  remove_unlabelled()

# read and filter split 2 data
integrated_s2 <- anndata::read_h5ad(par[["input_integrated_split2"]]) |>
  get_obs_var_for_integrated(unintegrated, split_id = 2) |>
  subset_markers_tocorrect() |>
  subset_nocontrols() |>
  remove_unlabelled()

print("Setting parameters\n")
donor_list <- unique(integrated_s1$obs$donor)
lineage_markers <- integrated_s1$var_names[
  integrated_s1$var$marker_type == "lineage"
]
n_clusters <- unintegrated$uns$parameter_num_clusters
grid_xdim <- unintegrated$uns$parameter_som_xdim
grid_ydim <- unintegrated$uns$parameter_som_ydim

print("Computing mapping similarity\n")
fs_mapping_similarity_allres <- list()
fsom_absdiff_by_donor_refsplit <- list()
fs_mapping_similarity_all <- c()

for (donor in donor_list) {
  print(paste0("extracting information for donor: ", donor, " (split 1)"))
  split1_data <- get_ff_annotations(
    integrated_s1,
    donor_name = donor,
    layer_name = "integrated"
  )
  print(paste0("extracting information for donor: ", donor, " (split 2)"))
  split2_data <- get_ff_annotations(
    integrated_s2,
    donor_name = donor,
    layer_name = "integrated"
  )
  
  # Direction A: build on split 1, map split 2 into it
  print("Building FlowSOM tree with split 1 cells")
  fs_tree_s1_ref <- FlowSOM(
    split1_data$flowframe,
    colsToUse = lineage_markers,
    nClus = n_clusters,
    xdim = grid_xdim,
    ydim = grid_ydim,
    seed = 42
  )
  
  print("Mapping split 2 cells in the tree")
  fs_tree_s2_mapped <- NewData(
    fs_tree_s1_ref,
    split2_data$flowframe
  )
  
  print("Computing mapping similarity (tree split 1; mapped split 2)")
  FSOM_mapping_refsplit1 <- compute_fs_mapping_similarity(
    fs_tree_s1_ref,
    split1_data$ct_annotations,
    fs_tree_s2_mapped,
    split2_data$ct_annotations
  )

  # Direction B: build on split 2, map split 1 into it
  print("Building FlowSOM tree with split 2 cells")
  fs_tree_s2_ref <- FlowSOM(
    split2_data$flowframe,
    colsToUse = lineage_markers,
    nClus = n_clusters,
    xdim = grid_xdim,
    ydim = grid_ydim,
    seed = 42
  )

  print("Mapping split 1 cells in the tree")
  fs_tree_s1_mapped <- NewData(
    fs_tree_s2_ref,
    split1_data$flowframe
  )

  print("Computing mapping similarity (tree split 2; mapped split 1)")
  FSOM_mapping_refsplit2 <- compute_fs_mapping_similarity(
    fs_tree_s2_ref,
    split2_data$ct_annotations,
    fs_tree_s1_mapped,
    split1_data$ct_annotations
  )

  # add full results for a donor to a list (both directions)
  fs_mapping_similarity_allres[[donor]] <- list(
    split_1 = FSOM_mapping_refsplit1,
    split_2 = FSOM_mapping_refsplit2
  )

  # store per-donor absdiff matrices, keyed by donor and reference split
  fsom_absdiff_by_donor_refsplit[[paste(donor, "split_1", sep = "_")]] <-
    FSOM_mapping_refsplit1$absdiff_matrix
  fsom_absdiff_by_donor_refsplit[[paste(donor, "split_2", sep = "_")]] <-
    FSOM_mapping_refsplit2$absdiff_matrix

  # add both similarity scores
  fs_mapping_similarity_all <- c(
    fs_mapping_similarity_all,
    FSOM_mapping_refsplit1$similarity,
    FSOM_mapping_refsplit2$similarity
  )
}

# compute average mapping similarity for all donors
fs_mapping_similarity_avg <- mean(fs_mapping_similarity_all)

print("Write output AnnData to file\n")
output <- anndata::AnnData(
  shape = c(0L, 0L),
  uns = list(
    dataset_id = integrated_s1$uns$dataset_id,
    method_id = integrated_s1$uns$method_id,
    metric_ids = "flowsom_mean_mapping_similarity",
    metric_values = fs_mapping_similarity_avg,
    fsom_absdiff_by_donor_refsplit = fsom_absdiff_by_donor_refsplit,
    fsom_parameters = list(
      "xdim" = grid_xdim,
      "ydim" = grid_ydim,
      "n_clusters" = n_clusters
    )
  )
)

output$write_h5ad(par[["output"]], compression = "gzip")
