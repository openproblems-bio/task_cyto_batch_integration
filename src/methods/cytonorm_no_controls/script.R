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
    name = "cytonorm_no_control",
    temp_dir = "resources_test/task_cyto_batch_integration/tmp",
    resources_dir = "src/utils"
)
## VIASH ENDs

source(paste0(meta$resources_dir, "/anndata_to_fcs.R"))

tmp_path <- meta[["temp_dir"]]

cat("Reading input files\n")
adata <- anndata::read_h5ad(par[["input"]])

cat("Creating aggregates per batch\n")

batches <- unique(adata$obs$batch)
fset_per_batch <- lapply(batches, function(bt) {
    anndata_to_fcs(adata[adata$obs$batch == bt, ])
})
names(fset_per_batch) <- batches

# create aggregate per batch
# TODO change to a parameter
n_cells_agg <- 1000
set.seed(42)
agg_per_batch <- lapply(batches, function(bt) {
    ff_obj <- FlowSOM::AggregateFlowFrames(
        fileNames = fset_per_batch[[bt]],
        cTotal = length(fset_per_batch[[bt]]) * n_cells_agg
    )
    # change the object name as otherwise it seems to just refer to the 1st sample
    # in the batch.. and it is unintuitive.
    description(ff_obj)$GUID <- paste0("batch", bt)
    return(ff_obj)
})

markers_to_correct <- as.vector(adata$var$channel[adata$var$to_correct])

lineage_markers <- as.vector(adata$var$channel[adata$var$marker_type == "lineage"])

cat("Training Cytonorm model using aggregates\n")

# FlowSOM.params and normParams are the default parameters in cytonorm
model <- CytoNorm::CytoNorm.train(
    files = flowSet(agg_per_batch),
    labels = batches,
    channels = markers_to_correct,
    outputDir = tmp_path,
    FlowSOM.params = list(
        nCells = 1000000,
        xdim = 15,
        ydim = 15,
        nClus = 10,
        scale = FALSE,
        colsToUse = lineage_markers
    ),
    transformList = NULL,
    normParams = list(nQ = 99, goal = "mean"),
    seed = 42,
    verbose = FALSE
)

cat("Normalising using Cytonorm model trained using aggregates\n")

# every sample, including the controls, pretty much the entire unintegrated data
# will be corrected.
fset_all <- anndata_to_fcs(adata)
batch_labs_per_samp <- sapply(sampleNames(fset_all), function(samp) {
    unique(adata[adata$obs$sample == samp]$obs$batch)[1]
})

norm_fset_all <- CytoNorm::CytoNorm.normalize(
    model = model,
    files = fset_all,
    labels = batch_labs_per_samp,
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
