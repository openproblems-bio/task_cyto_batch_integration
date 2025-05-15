library(anndata)
requireNamespace("cyCombine", quietly = TRUE)

## VIASH START
par <- list(
    input = "resources_test/task_cyto_batch_integration/leomazzi_cyto_spleen_subset/unintegrated_censored.h5ad",
    output = "resources_test/task_cyto_batch_integration/leomazzi_cyto_spleen_subset/output.h5ad"
)
meta <- list(name = "cycombine")
## VIASH END

cat("Reading input files\n")
input_adata <- anndata::read_h5ad(par[["input"]])

cat("Preparing input Anndata and df\n")

adata_to_correct <- input_adata[, input_adata$var$to_correct]

markers_to_correct <- input_adata$var_names[input_adata$var$to_correct]

lineage_markers <- as.vector(input_adata$var_names[
    input_adata$var$marker_type == "lineage"
])

df_to_correct <- as.data.frame(
    adata_to_correct$layers[["preprocessed"]],
    check.names = FALSE
)
df_to_correct$batch <- adata_to_correct$obs$batch
df_to_correct$sample <- adata_to_correct$obs$sample

cat("Run cyCombine\n")

# use the default parameters in normalize
# do the normalisation on all markers to be corrected
df_to_correct_norm <- cyCombine::normalize(
    df = df_to_correct,
    markers = markers_to_correct,
    norm_method = "scale",
    ties.method = "average"
)

# again, using default parameter values
cluster_labels <- cyCombine::create_som(
    df = df_to_correct_norm,
    markers = lineage_markers,
    rlen = 10,
    seed = 42,
    xdim = 8,
    ydim = 8
)

# Batch correct using default parameter values
df_corrected <- cyCombine::correct_data(
    df = df_to_correct_norm,
    label = cluster_labels,
    markers = markers_to_correct,
    method = "ComBat",
    covar = NULL,
    anchor = NULL,
    ref.batch = NULL,
    parametric = TRUE
)

cat("Preparing output Anndata\n")
df_not_corrected <- as.data.frame(
    input_adata[, !input_adata$var$to_correct]$layers[["preprocessed"]]
)

df_output <- cbind(
    df_corrected[, markers_to_correct],
    df_not_corrected
)

# reorder column to match input_adata
df_output <- df_output[, input_adata$var_names]

output <- anndata::AnnData(
    obs = input_adata$obs[, integer(0)],
    var = input_adata$var[colnames(df_output), integer(0)],
    layers = list(integrated = as.matrix(df_output)),
    uns = list(
        dataset_id = input_adata$uns$dataset_id,
        method_id = meta$name,
        parameters = list()
    )
)

cat("Write output AnnData to file\n")

output$write_h5ad(par[["output"]], compression = "gzip")
