requireNamespace("flowCore", quietly = TRUE)
requireNamespace("anndata", quietly = TRUE)
requireNamespace("Biobase", quietly = TRUE)
requireNamespace("CytoNorm", quietly = TRUE)

## VIASH START
par <- list(
    input = "resources_test/task_cyto_batch_integration/cyto_spleen_subset/unintegrated_censored.h5ad",
    output = "resources_test/output.h5ad",
    som_grid_size = 10,
    num_metacluster = 10,
    n_quantiles = 99
)
meta <- list(
    name = "cytonorm_control",
    temp_dir = "resources_test/task_cyto_batch_integration/tmp",
    resources_dir = "src/utils"
)
## VIASH END

source(paste0(meta$resources_dir, "/anndata_to_fcs.R"))

tmp_path <- meta[["temp_dir"]]

cat("Reading input files\n")
adata <- anndata::read_h5ad(par[["input"]])

cat("Preparing training data\n")

# get the control samples to be used for training the model
fset_train <- anndata_to_fcs(adata[adata$obs$is_control != 0, ])
# every sample, including the controls, pretty much the entire unintegrated data
# will be corrected.
fset_all <- anndata_to_fcs(adata)

cat("Setting up some variables for training the model\n")

# get batch label for the training data
batch_lab_train <- sapply(sampleNames(fset_train), function(samp) {
    unique(adata[adata$obs$sample == samp]$obs$batch)[1]
})

markers_to_correct <- as.vector(adata$var$channel[adata$var$to_correct])

lineage_markers <- as.vector(adata$var$channel[adata$var$marker_type == "lineage"])

# get number of cells for clustering.
# we will define this as the minimum of the smallest sample and 1,000,000.
# and multiply this by how many samples we have - because internally,
# this number is divided by the number of files to determine the amount to select from
# each individual file.
n_cells_per_control_sample <- flowCore::fsApply(fset_train, function(ff) nrow(exprs(ff)))
n_cells_for_clustering <- min(n_cells_per_control_sample, 1000000) * length(n_cells_per_control_sample)

cat("Training Cytonorm model using all control samples\n")

# FlowSOM.params and normParams are the default parameters in cytonorm
model <- CytoNorm::CytoNorm.train(
    files = fset_train,
    labels = batch_lab_train,
    channels = markers_to_correct,
    outputDir = tmp_path,
    FlowSOM.params = list(
        nCells = n_cells_for_clustering,
        xdim = par[["som_grid_size"]],
        ydim = par[["som_grid_size"]],
        nClus = par[["num_metacluster"]],
        scale = FALSE,
        colsToUse = lineage_markers
    ),
    transformList = NULL,
    normParams = list(nQ = par[["n_quantiles"]], goal = "mean"),
    seed = 42,
    verbose = FALSE
)

# get batch label for the validation data
batch_labs <- sapply(sampleNames(fset_all), function(samp) {
    unique(adata[adata$obs$sample == samp]$obs$batch)[1]
})

cat("Normalising using Cytonorm model trained using all control samples\n")

norm_fset_all <- CytoNorm::CytoNorm.normalize(
    model = model,
    files = fset_all,
    labels = batch_labs,
    transformList = NULL,
    transformList.reverse = NULL,
    outputDir = tmp_path,
    prefix = "Norm_",
    clean = TRUE,
    write = FALSE,
    verbose = FALSE
)

cat("Preparing output anndata\n")
# cytonorm will return all markers corrected or not in the same order as the input data.
# so we can just directly replace the colnames with var_names
norm_mat <- flowCore::fsApply(norm_fset_all, exprs)
colnames(norm_mat) <- adata$var_names

norm_mat <- anndata::AnnData(
    obs = adata$obs[, integer(0)],
    var = adata$var[colnames(norm_mat), integer(0)],
    layers = list(integrated = norm_mat),
    uns = list(
        dataset_id = adata$uns$dataset_id,
        method_id = meta$name,
        parameters = list()
    )
)

cat("Write output AnnData to file\n")
norm_mat$write_h5ad(par[["output"]], compression = "gzip")
