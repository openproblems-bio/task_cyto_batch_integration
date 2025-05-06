library(flowCore)
library(anndata)
library(Biobase)
library(CytoNorm)

## VIASH START
par <- list(
    input = "resources_test/task_cyto_batch_integration/cyto_spleen_subset/unintegrated_censored.h5ad",
    output = "resources_test/output.h5ad"
)
meta <- list(
    name = "cytonorm_control",
    temp_dir = "resources_test/task_cyto_batch_integration/tmp"
)
## VIASH ENDs

source(paste0(meta$resources_dir, "/anndata_to_fcs.R"))

tmp_path <- meta[["temp_dir"]]

cat("Reading input files\n")
adata <- anndata::read_h5ad(par[["input"]])

# get the control samples to be used for training the model
fset_train <- anndata_to_fcs(adata[adata$obs$is_control != 0, ])
# every sample, including the controls, pretty much the entire unintegrated data 
# will be corrected.
fset_all <- anndata_to_fcs(adata)

# get batch label for the training data
batch_lab_train <- sapply(sampleNames(fset_train), function(samp) {
    unique(adata[adata$obs$sample == samp]$obs$batch)[1]
})

markers_to_correct <- as.vector(adata$var$channel[adata$var$to_correct])

# FlowSOM.params and normParams are the default parameters in cytonorm
model <- CytoNorm.train(
    files = fset_train,
    labels = batch_lab_train,
    channels = markers_to_correct,
    outputDir = tmp_path,
    FlowSOM.params = list(
        nCells = 1000000,
        xdim = 15,
        ydim = 15,
        nClus = 10,
        scale = FALSE
    ),
    transformList = NULL,
    normParams = list(nQ = 99, goal = "mean"),
    seed = 42,
    verbose = FALSE
)

# get batch label for the validation data
batch_labs <- sapply(sampleNames(fset_all), function(samp) {
    unique(adata[adata$obs$sample == samp]$obs$batch)[1]
})

norm_fset_all <- CytoNorm.normalize(
    model = model,
    files = fset_all,
    labels = batch_labs,
    transformList = NULL,
    transformList.reverse = NULL,
    outputDir = tmp_path,
    prefix = "Norm_",
    clean = FALSE,
    write = FALSE,
    verbose = FALSE
)

cat("Preparing output anndata\n")
# cytonorm will return all markers corrected or not in the same order as the input data.
# so we can just directly replace the colnames with var_names
norm_mat <- fsApply(norm_fset_all, exprs)
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
