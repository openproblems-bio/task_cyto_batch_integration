library(anndata)
requireNamespace("sva", quietly = TRUE)

## VIASH START
par <- list(
  input = "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/censored_split1.h5ad",
  output = "output.h5ad",
  target = "mid"
)
meta <- list(
  name = "combat",
  resources_dir = "src/utils"
)
## VIASH END

source(paste0(meta$resources_dir, "/helper_functions.R"))

cat("Reading input files\n")
input <- anndata::read_h5ad(par[["input"]]) |>
  subset_nocontrols()

cat("Subset data\n")
data_not_correct <- input[, !input$var$to_correct]
data_to_correct <- input[, input$var$to_correct]

# correct towards batch 1 (goal) or towards the pooled mean across batches (mid)
ref_batch <- if (par[["target"]] == "goal") "1" else NULL

cat("Run ComBat\n")
output <- sva::ComBat(
  dat = Matrix::t(data_to_correct$layers[["preprocessed"]]),
  batch = data_to_correct$obs$batch,
  ref.batch = ref_batch
)

cat("Preparing output Anndata\n")
output <- Matrix::t(output)
output <- cbind(output, data_not_correct$layers[["preprocessed"]])

cat("Write output AnnData to file\n")
output <- anndata::AnnData(
  obs = input$obs[, integer(0)],
  var = input$var[colnames(output), integer(0)],
  layers = list(integrated = output),
  uns = list(
    dataset_id = input$uns$dataset_id,
    method_id = meta$name,
    parameters = list(target = par[["target"]])
  )
)

# reorder the var in output
output <- output[, input$var_names]

output$write_h5ad(par[["output"]], compression = "gzip")

cat("Written anndata of shape ", dim(output), " to file: ", par[["output"]], "\n")
