requireNamespace("cyCombine", quietly = TRUE)
requireNamespace("anndata", quietly = TRUE)

## VIASH START
par <- list(
    input = "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/censored_split1.h5ad",
    output = "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/output.h5ad",
    controls = "all",
    target = "mid",
    som_grid_size = 8,
    rlen = 10
)
meta <- list(
    name = "cycombine",
    resources_dir = "src/utils"
)
## VIASH END

source(paste0(meta$resources_dir, "/helper_functions.R"))

cat("Reading input files\n")
input_adata <- anndata::read_h5ad(par[["input"]])

# subset the input depending on which controls are used
if (par[["controls"]] == "none") {
    input_adata <- subset_nocontrols(input_adata)
} else if (par[["controls"]] == "one") {
    input_adata <- subset_onecontrol(input_adata)
}

cat("Preparing input Anndata and df\n")

markers_to_correct <- input_adata$var_names[input_adata$var$to_correct]

adata_to_correct <- input_adata[, markers_to_correct]

# convert adata to data.frame

df_to_correct <- as.data.frame(
    adata_to_correct$layers[["preprocessed"]]
)
df_to_correct$batch <- as.factor(adata_to_correct$obs$batch)
df_to_correct$sample <- as.factor(adata_to_correct$obs$sample)

# add an "anchor" column which specify which samples are the technical replicates
# this is a bit weird as the anchor information should be a column name or character vector
# giving all anchors the same label and every other sample a unique label.
# so for our controls (replicate samples), i have to give it the same name,
# while non-replicate sample should be given a unique identifier.
# hence for non-replicate samples, we will just use the values in the "sample" column.
# for the replicate samples, we will use the values in the "is_control" column.
# that way, samples from the same donor will have the same label, and
# cycombine can identify which samples are technical replicates of each other.
# the as.character function is needed as otherwise we will get NAs for those controls..

anchor <- NULL
if (par[["controls"]] != "none") {
    df_to_correct$anchor <- ifelse(
        adata_to_correct$obs$is_control == 0,
        as.character(adata_to_correct$obs$sample),
        paste0("control_", adata_to_correct$obs$is_control)
    )
    df_to_correct$anchor <- as.factor(df_to_correct$anchor)
    anchor <- "anchor"
}

# correct towards batch 1 (goal) or towards a midpoint across batches (mid)
ref_batch <- if (par[["target"]] == "goal") "1" else NULL

lineage_markers <- as.vector(input_adata$var_names[
    input_adata$var$marker_type == "lineage"
])


cat("Run cyCombine\n")

# use the default parameters in normalize
# do the normalisation on all markers to be corrected
# use z-score normalisation because we are merging batches
# from a single study.
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
    rlen = par[["rlen"]],
    seed = 42,
    xdim = par[["som_grid_size"]],
    ydim = par[["som_grid_size"]]
)

# Batch correct using default parameter values
df_corrected <- cyCombine::correct_data(
    df = df_to_correct,
    label = cluster_labels,
    markers = markers_to_correct,
    method = "ComBat",
    covar = NULL,
    anchor = anchor,
    ref.batch = ref_batch,
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
        parameters = list(
            "controls" = par[["controls"]],
            "target" = par[["target"]],
            "normalize" = list(
                "markers" = markers_to_correct,
                "norm_method" = "scale",
                "ties.method" = "average"
            ),
            "create_som" = list(
                "markers" = lineage_markers,
                "seed" = 42,
                "rlen" = par[["rlen"]],
                "xdim" = par[["som_grid_size"]],
                "ydim" = par[["som_grid_size"]]
            ),
            "correct_data" = list(
                "markers" = markers_to_correct,
                "method" = "ComBat",
                "covar" = NULL,
                "anchor" = anchor,
                "ref.batch" = ref_batch,
                "parametric" = TRUE
            )
        )
    )
)

cat("Write output AnnData to file\n")

output$write_h5ad(par[["output"]], compression = "gzip")

cat("Written anndata of shape ", dim(output), " to file: ", par[["output"]], "\n")
