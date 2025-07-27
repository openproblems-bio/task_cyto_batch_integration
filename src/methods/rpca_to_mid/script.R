requireNamespace("anndata", quietly = TRUE)
requireNamespace("Seurat", quietly = TRUE)

## VIASH START
par <- list(
  input = "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/unintegrated_censored.h5ad",
  output = "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/output.h5ad",
  npcs = 50,
  n_neighbours = 5
)
meta <- list(
  name = "rpca_to_mid"
)
## VIASH END

cat("Reading input files\n")
input_adata <- anndata::read_h5ad(par[["input"]])

cat("Preparing input Anndata and df\n")
adata_to_correct <- input_adata[, input_adata$var$to_correct]
markers_to_correct <- input_adata$var_names[input_adata$var$to_correct]

cat("Creating Seurat object and preprocess\n")
# create one seurat object per batch
batches <- unique(input_adata$obs$batch)

seurat_objs <- lapply(batches, function(batch) {

  cat(paste("Processing batch", batch))
  # batch <- batches[1]
  mat <- Matrix::as.matrix(input_adata[
    input_adata$obs$batch == batch,
    input_adata$var$to_correct
  ]$layers["preprocessed"])

  # have to transpose so cells are columns..
  mat <- t(mat)

  # convert to sparse matrix
  mat <- Matrix::Matrix(mat, sparse = TRUE)

  seurat_obj <- Seurat::CreateSeuratObject(
      counts = mat,
      data = mat,
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

n_pcs_computed <- dim(
  Seurat::Embeddings(
    seurat_objs[[1]], reduction = "pca"
  )
)[2]
# setting scale to FALSE as we have scaled the features.
anchors <- Seurat::FindIntegrationAnchors(
    object.list = seurat_objs,
    anchor.features = markers_to_correct,
    dims = seq(n_pcs_computed),
    k.anchor = par[["n_neighbours"]],
    reduction = "rpca",
    scale = FALSE
)

cat("Write output AnnData to file\n")
output <- anndata::AnnData(
  
)
output$write_h5ad(par[["output"]], compression = "gzip")
