library(anndata)
library(batchelor)

## VIASH START
par <- list(
  input = "resources_test/.../input.h5ad",
  output = "output.h5ad"
)
meta <- list(
  name = "mnn"
)
## VIASH END

cat("Reading input files\n")
input <- anndata::read_h5ad(par[["input"]])

cat("Subset data\n")
data_not_correct <- input[, !input$var$to_correct]
data_to_correct <- input[, input$var$to_correct]

cat("Run MNN\n")

# If prop_num_nn is set to 0, it is set to NULL
if (par[["prop_num_nn"]] == 0.0) {
  print("'prop_num_nn' = 0, setting it to NULL")
  par[["prop_num_nn"]] <- NULL
}

corrected_data <- mnnCorrect(Matrix::t(data_to_correct$layers[["preprocessed"]]),
              batch = data_to_correct$obs$batch,
              k= par[["num_nn"]],
              prop.k = par[["prop_num_nn"]],
              sigma = par[["sigma_value"]],
              cos.norm.in = FALSE,
              cos.norm.out = FALSE)

cat("Preparing output Anndata\n")
corrected_data <- Matrix::t(assay(corrected_data))
corrected_data <- cbind(corrected_data, data_not_correct$layers[["preprocessed"]])

cat("Write output AnnData to file\n")
output <- anndata::AnnData(
  obs = input$obs[, integer(0)],
  var = input$var[colnames(corrected_data), integer(0)],
  layers = list(integrated = corrected_data),
  uns = list(
    dataset_id = input$uns$dataset_id,
    method_id = meta$name,
    parameters = list(
      "k" = par[["num_nn"]],
      "prop.k" = par["prop_num_nn"],
      "sigma" = par[["sigma_value"]]
    )
  )
)

output$write_h5ad(par[["output"]], compression = "gzip")
