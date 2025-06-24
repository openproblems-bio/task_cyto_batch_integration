library(anndata, quietly = TRUE)
library(FlowSOM, quietly = TRUE)

## VIASH START
par <- list(
  input_validation = "resources_test/task_cyto_batch_integration/leomazzi_cyto_spleen_subset/validation.h5ad",
  input_unintegrated = "resources_test/task_cyto_batch_integration/leomazzi_cyto_spleen_subset/unintegrated.h5ad",
  input_integrated = "resources_test/task_cyto_batch_integration/leomazzi_cyto_spleen_subset/integrated.h5ad",
  output = "output.h5ad"
)
meta <- list(
  name = "flowsom_mapping_similarity",
  resources_dir = "src/metrics/flowsom_mapping_similarity",
  temp_dir = "/tmp"
)
## VIASH END

source(paste0(meta$resources_dir, "/helper.R"))
source(paste0(meta$resources_dir, "/helper_functions.R"))

print("Preprocess input files with python")
input_unintegrated <- anndata::read_h5ad(par[["input_unintegrated"]])

# read validation data
input_validation <- anndata::read_h5ad(par[["input_validation"]])

# read and filter integrated data
input_integrated <- anndata::read_h5ad(par[["input_integrated"]]) |>
  get_obs_var_for_integrated(
    input_validation,
    input_unintegrated
  ) |>
  subset_nocontrols() |>
  remove_unlabelled()

# filter validation data
input_validation <- input_validation |>
  remove_unlabelled()

print("Setting parameters\n")
donor_list <- unique(input_integrated$obs$donor)
lineage_markers <- input_validation$var_names[
  input_validation$var$marker_type == "lineage"
]
n_clusters <- input_integrated$uns$parameter_flowsom_nclus
grid_xdim <- input_integrated$uns$parameter_flowsom_xdim
grid_ydim <- input_integrated$uns$parameter_flowsom_ydim

print("Computing mapping similarity\n")
fs_mapping_similarity_allres <- list()
fs_mapping_similarity_all <- c()

for (donor in donor_list) {
  print(paste0("extracting information for donor: ", donor, " (integrated)"))
  integrated_data <- get_ff_annotations(
    input_integrated,
    donor_name = donor,
    layer_name = "integrated"
  )
  print(paste0("extracting information for donor: ", donor, " (validation)"))
  validation_data <- get_ff_annotations(
    input_validation,
    donor_name = donor,
    layer_name = "preprocessed"
  )

  print("Building FlowSOM tree with validation cells")
  fs_tree_validation <- FlowSOM(
    validation_data$flowframe,
    colsToUse = lineage_markers,
    nClus = n_clusters,
    xdim = grid_xdim,
    ydim = grid_ydim,
    seed = 42
  )

  print("Mapping new cells in the tree")
  fs_tree_integrated <- NewData(
    fs_tree_validation,
    integrated_data$flowframe
  )

  print("Computing mapping similarity")
  FSOM_mapping <- compute_fs_mapping_similarity(
    fs_tree_integrated,
    integrated_data$ct_annotations,
    fs_tree_validation,
    validation_data$ct_annotations
  )
  # add full results for a donor to a list
  fs_mapping_similarity_allres[[donor]] <- FSOM_mapping
  # add only similarity metric to a variable
  fs_mapping_similarity_all <- c(
    fs_mapping_similarity_all,
    FSOM_mapping$similarity
  )
}

# compute average mapping similarity for all donors
fs_mapping_similarity_avg <- mean(fs_mapping_similarity_all)

print("Write output AnnData to file\n")
output <- anndata::AnnData(
  shape = c(0L, 0L),
  uns = list(
    dataset_id = input_integrated$uns$dataset_id,
    method_id = input_integrated$uns$method_id,
    metric_ids = "flowsom_mean_mapping_similarity",
    metric_values = fs_mapping_similarity_avg,
    fsom_parameters = list(
      "xdim" = grid_xdim,
      "ydim" = grid_ydim,
      "n_clusters" = n_clusters
    )
  )
)

cat("Output AnnData object:\n")
print(output)

output$write_h5ad(par[["output"]], compression = "gzip")
