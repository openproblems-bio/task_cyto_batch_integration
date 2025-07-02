library(anndata)

## VIASH START
par <- list(
  input_validation = "resources_test/.../validation.h5ad",
  input_unintegrated = "resources_test/.../unintegrated.h5ad",
  input_integrated = "resources_test/.../integrated.h5ad",
  output = "output.h5ad"
)
meta <- list(
  name = "cms"
)
## VIASH END

cat("Reading input files\n")
input_validation <- anndata::read_h5ad(par[["input_validation"]])
input_unintegrated <- anndata::read_h5ad(par[["input_unintegrated"]])
input_integrated <- anndata::read_h5ad(par[["input_integrated"]])

cat("Compute metrics\n")
# metric_ids and metric_values can have length > 1
# but should be of equal length
uns_metric_ids <- c("cms")
uns_metric_values <- c(0.5)

cat("Write output AnnData to file\n")
output <- anndata::AnnData(
  
)
output$write_h5ad(par[["output"]], compression = "gzip")
