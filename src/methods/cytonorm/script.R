requireNamespace("flowCore", quietly = TRUE)
requireNamespace("anndata", quietly = TRUE)
requireNamespace("Biobase", quietly = TRUE)
requireNamespace("CytoNorm", quietly = TRUE)
requireNamespace("FlowSOM", quietly = TRUE)

## VIASH START
par <- list(
    input = "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/censored_split1.h5ad",
    output = "resources_test/output.h5ad",
    controls = "all",
    target = "mid",
    som_grid_size = 10,
    num_metacluster = 10,
    n_quantiles = 99
)
meta <- list(
    name = "cytonorm",
    temp_dir = "resources_test/task_cyto_batch_integration/tmp",
    resources_dir = "src/utils"
)
## VIASH END

source(paste0(meta$resources_dir, "/anndata_to_fcs.R"))
source(paste0(meta$resources_dir, "/helper_functions.R"))

tmp_path <- get_temp_dir(meta)
print(paste0("Using temp dir: ", tmp_path))
on.exit(clean_temp_dir(tmp_path))

cat("Reading input files\n")
adata <- anndata::read_h5ad(par[["input"]])

# subset the input depending on which controls are used
if (par[["controls"]] == "none") {
    adata <- subset_nocontrols(adata)
} else if (par[["controls"]] == "one") {
    adata <- subset_onecontrol(adata)
}

cat("Preparing training data\n")

if (par[["controls"]] == "none") {
    # without controls, train on one aggregate per batch as pseudo-reference
    cat("Creating aggregates per batch\n")

    batches <- as.character(unique(adata$obs$batch))
    fset_per_batch <- lapply(batches, function(bt) {
        anndata_to_fcs(adata[adata$obs$batch == as.numeric(bt), ])
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
    fset_train <- flowCore::flowSet(agg_per_batch)
    batch_lab_train <- batches
} else {
    # get the control samples to be used for training the model
    fset_train <- anndata_to_fcs(adata[adata$obs$is_control != 0, ])

    # get batch label for the training data
    batch_lab_train <- vapply(sampleNames(fset_train), function(samp) {
        as.character(
            unique(
                adata[adata$obs$sample == samp]$obs$batch
            )[1]
        )
    }, FUN.VALUE = character(1))
}

# every sample, including the controls, pretty much the entire unintegrated data
# will be corrected.
fset_all <- anndata_to_fcs(adata)

# get batch label for the all data
batch_labs <- vapply(sampleNames(fset_all), function(samp) {
    as.character(
        unique(
            adata[adata$obs$sample == samp]$obs$batch
        )[1]
    )
}, FUN.VALUE = character(1))

cat("Setting up some variables for training the model\n")

markers_to_correct <- as.vector(adata$var$channel[adata$var$to_correct])

lineage_markers <- as.vector(adata$var$channel[adata$var$marker_type == "lineage"])

# get number of cells for clustering.
# we will define this as the minimum of the smallest training sample and 1,000,000.
# and multiply this by how many training samples we have - because internally,
# this number is divided by the number of files to determine the amount to select from
# each individual file.
n_cells_per_train_sample <- flowCore::fsApply(fset_train, function(ff) nrow(exprs(ff)))
n_cells_for_clustering <- min(n_cells_per_train_sample, 1000000) * length(n_cells_per_train_sample)

# align distributions to batch 1 (goal) or to the mean across batches (mid)
norm_goal <- if (par[["target"]] == "goal") "1" else "mean"

cat("Training Cytonorm model\n")

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
    normParams = list(
        nQ = par[["n_quantiles"]],
        goal = norm_goal
    ),
    seed = 42,
    verbose = FALSE,
    recompute = TRUE
)

cat("Normalising using trained Cytonorm model\n")

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
        parameters = list(
            controls = par[["controls"]],
            target = par[["target"]]
        )
    )
)

cat("Write output AnnData to file\n")
norm_mat$write_h5ad(par[["output"]], compression = "gzip")

cat("Written anndata of shape ", dim(norm_mat), " to file: ", par[["output"]], "\n")
