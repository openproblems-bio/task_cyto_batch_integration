library(anndata)
library(anndataR) #needed to write FSOM trees in the output without errors
library(FlowSOM)

## VIASH START
par <- list(
  input_validation = "resources_test/.../validation.h5ad",
  input_unintegrated = "resources_test/.../unintegrated.h5ad",
  input_integrated = "resources_test/.../integrated.h5ad",
  output = "output.h5ad"
)
meta <- list(
  name = "flowsom_mapping_similarity",
  temp_dir = "/tmp"
)
## VIASH END

source(paste0(meta$resources_dir, "/helper.R"))
preprocess_script <- paste0(meta$resources_dir, "/preprocess_with_utils.py")
preprocess_functions <- paste0(meta$resources_dir, "/helper_functions.py")

print("Preprocess input files with python")
int_file <- par[["input_integrated"]]
unint_file <- par[["input_unintegrated"]]
val_file <- par[["input_validation"]]
system(paste("python", preprocess_script, preprocess_functions, int_file, unint_file, val_file ))

cat("Reading preprocessed files\n")
input_integrated <- anndata::read_h5ad("tmp/integrated.h5ad")
input_validation <- anndata::read_h5ad("tmp/validation.h5ad")

print("Setting parameters\n")
donor_list <- unique(input_integrated$obs$donor)
lineage_markers <- input_validation$var_names[input_validation$var$marker_type=="lineage"] 

#TODO: make it dataset specific
n_clusters <- 25
grid_xdim <- 15
grid_ydim <- 15

print("Computing mapping similarity\n")
fs_mapping_similarity_allres <- list()
fs_mapping_similarity_all <- c()

for (donor in donor_list) {
  print(paste0("extracting information for donor: ",donor, " (integrated)"))
  integrated_data <- get_ff_annotations(input_integrated,
                                        donor_name = donor,
                                        layer_name = "integrated")
  print(paste0("extracting information for donor: ",donor, " (validation)"))
  validation_data <- get_ff_annotations(input_validation,
                                        donor_name = donor,
                                        layer_name = "preprocessed")
  
  ff_integrated <- integrated_data$flowframe
  ct_annotations_integrated <-integrated_data$ct_annotations
  ff_validation <- validation_data$flowframe
  ct_annotations_validation <-validation_data$ct_annotations
  
  print("Building FlowSOM tree with validation cells")
  fs_tree_validation <- FlowSOM(ff_validation,
                                colsToUse = lineage_markers,
                                nClus = n_clusters,
                                xdim = grid_xdim, ydim = grid_ydim,
                                seed = 42)
  
  print("Mapping new cells in the tree")
  fs_tree_integrated <- NewData(fs_tree_validation,ff_integrated)
  
  print("Computing mapping similarity")
  FSOM_mapping <- compute_fs_mapping_similarity(
    fs_tree_integrated,
    ct_annotations_integrated,
    fs_tree_validation,
    ct_annotations_validation
  )
  #add full results for a donor to a list
  fs_mapping_similarity_allres[[donor]] = FSOM_mapping 
  #add only similarity metric to a variable
  fs_mapping_similarity_all <- c(fs_mapping_similarity_all,FSOM_mapping$similarity)
}

fs_mapping_similarity_avg <- mean(fs_mapping_similarity_all)
print(fs_mapping_similarity_avg)
print("Write output AnnData to file\n")
output <- anndataR::AnnData(
  shape = c(0, 0),
  uns = list(
    dataset_id = input_integrated$uns$dataset_id,
    method_id = input_integrated$uns$method_id,
    metric_ids = "FlowSOM_mapping_similarity",
    metric_values = fs_mapping_similarity_avg,
    fsom_parameters = list("xdim" = grid_xdim,
                           "ydim" = grid_ydim,
                           "n_clusters" = n_clusters)
 )
)
output$write_h5ad(par[["output"]], compression = "gzip")