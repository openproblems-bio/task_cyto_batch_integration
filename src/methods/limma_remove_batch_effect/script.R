library(anndata)
requireNamespace("limma", quietly = TRUE)

## VIASH START
par <- list(
  input = "resources_test/task_cyto_batch_integration/starter_file/unintegrated.h5ad",
  output = "output.h5ad"
)
meta <- list(
  name = "limma_remove_batch_effect"
)
## VIASH END

cat("Reading input files\n")
input <- anndata::read_h5ad(par[["input"]])

cat("Subset data\n")
input <- input[, input$var$to_correct]

cat("Run limma\n")
output <- limma::removeBatchEffect(
  Matrix::t(input$layers[["preprocessed"]]),
  input$obs$batch
)

cat("Write output AnnData to file\n")
output <- anndata::AnnData(
  obs = input$obs[, integer(0)],
  var = input$var[, integer(0)],
  layers = list(
    integrated = Matrix::t(output)
  ),
  uns = list(
    dataset_id = input$uns$dataset_id,
    method_id = meta$name,
    parameters = list()
  )
)

output$write_h5ad(par[["output"]], compression = "gzip")
