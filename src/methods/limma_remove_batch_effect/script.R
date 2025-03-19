library(anndata)
requireNamespace("limma", quietly = TRUE)

## VIASH START
par <- list(
  input = "resources_test/task_cyto_batch_integration/cyto_spleen_subset/unintegrated_censored.h5ad",
  output = "output.h5ad"
)
meta <- list(
  name = "limma_remove_batch_effect"
)
## VIASH END

cat("Reading input files\n")
input <- anndata::read_h5ad(par[["input"]])

cat("Subset data\n")
data_not_correct <- input[, !input$var$to_correct]
data_to_correct <- input[, input$var$to_correct]

cat("Run limma\n")
output <- limma::removeBatchEffect(
  Matrix::t(data_to_correct$layers[["preprocessed"]]),
  data_to_correct$obs$batch
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
    parameters = list()
  )
)

# reorder the var in output
output <- output[, input$var_names]

output$write_h5ad(par[["output"]], compression = "gzip")
