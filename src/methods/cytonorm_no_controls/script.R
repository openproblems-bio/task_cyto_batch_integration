requireNamespace("flowCore", quietly = TRUE)
requireNamespace("anndata", quietly = TRUE)
requireNamespace("Biobase", quietly = TRUE)
requireNamespace("CytoNorm", quietly = TRUE)
requireNamespace("FlowSOM", quietly = TRUE)

## VIASH START
par <- list(
    input = "resources_test/task_cyto_batch_integration/cyto_spleen_subset/unintegrated_censored.h5ad",
    output = "resources_test/output.h5ad",
    som_grid_size = 10,
    num_metacluster = 10,
    n_quantiles = 99
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
set.seed(42)
agg_per_batch <- lapply(batches, function(bt) {

    n_cells_per_sample <- flowCore::fsApply(fset_per_batch[[bt]], function(ff) nrow(exprs(ff)))
    n_cells_to_aggregate <- min(n_cells_per_sample, 1000000) * length(n_cells_per_sample)

    ff_obj <- FlowSOM::AggregateFlowFrames(
        fileNames = fset_per_batch[[bt]],
        cTotal = n_cells_to_aggregate
    )
    # change the object name as otherwise it seems to just refer to the 1st sample
    # in the batch.. and it is unintuitive.
    # TODO remove the suppress warning for keyword when flowCore fixed it.
    suppressWarnings(description(ff_obj)$GUID <- paste0("batch", bt))
    return(ff_obj)
})
# convert to flowSet
agg_per_batch <- flowCore::flowSet(agg_per_batch)

cat("Setting up some variables for training the model\n")

markers_to_correct <- as.vector(adata$var$channel[adata$var$to_correct])

lineage_markers <- as.vector(adata$var$channel[adata$var$marker_type == "lineage"])

# get number of cells for clustering.
# we will define this as the minimum of the smallest aggregate and 1,000,000.
# and multiply this by how many batches we have.
n_cells_per_batch <- flowCore::fsApply(agg_per_batch, function(ff) nrow(exprs(ff)))
n_cells_for_clustering <- min(n_cells_per_batch, 1000000) * length(n_cells_per_batch)

cat("Training Cytonorm model using aggregates\n")

# FlowSOM.params and normParams are the default parameters in cytonorm
model <- CytoNorm::CytoNorm.train(
    files = agg_per_batch,
    labels = batches,
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
