requireNamespace("anndata", quietly = TRUE)
requireNamespace("Seurat", quietly = TRUE)

## VIASH START
par <- list(
    input = "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/unintegrated_censored.h5ad",
    output = "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/output.h5ad",
    npcs = 10,
    n_neighbours = 50
)
meta <- list(
    name = "rpca_to_goal"
)
## VIASH END

cat("Reading input files\n")
input_adata <- anndata::read_h5ad(par[["input"]])

cat("Preparing input Anndata\n")
input_adata$obs$batch <- as.factor(input_adata$obs$batch)

adata_to_correct <- input_adata[, input_adata$var$to_correct]
markers_to_correct <- input_adata$var_names[input_adata$var$to_correct]

cat("Creating Seurat object and preprocess\n")

# create one seurat object per batch
batches <- unique(input_adata$obs$batch)

seurat_objs <- lapply(batches, function(batch) {

    cat(paste("Processing batch", batch))

    adata_batch <- input_adata[
        input_adata$obs$batch == batch,
        input_adata$var$to_correct
    ]
    # batch <- batches[1]
    mat <- Matrix::as.matrix(adata_batch$layers["preprocessed"])

    # have to transpose so cells are columns..
    mat <- Matrix::t(mat)

    # convert to sparse matrix
    mat <- Matrix::Matrix(mat, sparse = TRUE)

    seurat_obj <- Seurat::CreateSeuratObject(
        counts = mat,
        # data = mat,
        assay = "cyto",
        meta.data = adata_batch$obs
    )

    SeuratObject::SetAssayData(
        object = seurat_obj,
        slot = "data",
        new.data = mat,
        assay = "cyto"
    )

    # save RAM
    rm(mat)

    # scale all features/markers
    seurat_obj <- Seurat::ScaleData(
        object = seurat_obj,
        features = markers_to_correct,
        assay = "cyto",
        verbose = FALSE
    )

    # run pca. mandatory
    # if num pcs is more than number of markers, it'll be capped at
    # the number of markers
    # not using approximate pca as we don't have many markers
    seurat_obj <- Seurat::RunPCA(
        object = seurat_obj,
        features = markers_to_correct,
        assay = "cyto",
        npcs = par[["npcs"]],
        approx = FALSE,
        verbose = FALSE
    )

    return(seurat_obj)
})

names(seurat_objs) <- batches

cat("Finding anchors\n")

# get how many PCs we have calculated
# if the number of PCs is more than how many markers
# seurat set that to how many markers.
# hence we can't just use par[["npcs"]] below.
npcs_computed <- dim(seurat_objs[[1]][["pca"]])[2]

anchors <- Seurat::FindIntegrationAnchors(
    object.list = seurat_objs,
    anchor.features = markers_to_correct,
    dims = seq(npcs_computed),
    k.anchor = par[["n_neighbours"]],
    reduction = "rpca",
    verbose = FALSE,
    reference = which(names(seurat_objs) == "1")
)

cat("Batch correct\n")

# Warning will say Layer counts isn't present in the assay object; returning NULL
# Even though the original assay has counts layer.
# Not sure why. But the object has data layer.
batch_corrected_seurat_obj <- Seurat::IntegrateData(
    anchorset = anchors,
    features = markers_to_correct,
    features.to.integrate = markers_to_correct,
    dims = seq(npcs_computed),
    verbose = FALSE
)
# just to be sure!
Seurat::DefaultAssay(batch_corrected_seurat_obj) <- "integrated"

cat("Creating output AnnData\n")

batch_corrected_mat <- Matrix::t(
    Matrix::as.matrix(Seurat::GetAssayData(batch_corrected_seurat_obj))
)
# cbind corrected matrix to matrix containing markers not corrected
batch_corrected_mat <- cbind(
    batch_corrected_mat,
    input_adata[, !input_adata$var$to_correct]$layers[["preprocessed"]]
)

# make sure the row and column orders are matching
# between input adata and the batch corrected matrix
batch_corrected_mat <- batch_corrected_mat[
    input_adata$obs_names, input_adata$var_names
]

cat("Write output AnnData to file\n")
output <- anndata::AnnData(
    obs = input_adata$obs[, integer(0)],
    var = input_adata$var[colnames(batch_corrected_mat), integer(0)],
    layers = list(integrated = batch_corrected_mat),
    uns = list(
        dataset_id = input_adata$uns$dataset_id,
        method_id = meta$name,
        parameters = list(
            "npcs" = par[["npcs"]],
            "n_neighbours" = par[["n_neighbours"]]
        )
    )
)
output$write_h5ad(par[["output"]], compression = "gzip")